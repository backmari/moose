/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "WeibullStressAtCrackFrontEdge.h"
#include "SymmTensor.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "PlaneTracing.h"
#include <cmath>
#include <algorithm>
#include <set>

template<>
InputParameters validParams<WeibullStressAtCrackFrontEdge>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params += validParams<MaterialTensorCalculator>();

  params.addRequiredParam<UserObjectName>("crack_front_definition", "The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index", "The index of the point on the crack front");
  params.addRequiredParam<Real>("yield_stress", "Yield stress of the material");
  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addParam<Real>("weibull_r_max","Max radius for Weibull stress calculation");
  params.addParam<unsigned int>("symmetry_plane", "Account for a symmetry plane passing through the plane of the crack, normal to the specified axis (0=x, 1=y, 2=z)");
  params.set<MooseEnum>("quantity") = "MaxPrincipal";
  params.set<MultiMooseEnum>("execute_on") = "timestep_end";

  return params;
}

WeibullStressAtCrackFrontEdge :: WeibullStressAtCrackFrontEdge(const InputParameters & parameters) :
  ElementIntegralPostprocessor(parameters),
  MaterialTensorCalculator(parameters),
  _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
  _has_crack_front_point_index(isParamValid("crack_front_point_index")),
  _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
  _stress_tensor(getMaterialProperty<SymmTensor>("stress")),
  _m(getParam<Real>("m")),
  _lambda(getParam<Real>("lambda")),
  _yield_stress(getParam<Real>("yield_stress")),
  _r_max(getParam<Real>("weibull_r_max")),
  _has_symmetry_plane(isParamValid("symmetry_plane"))
{
  _cutoff = _lambda * _yield_stress;
}

void
WeibullStressAtCrackFrontEdge::initialSetup()
{
  _treat_as_2d = _crack_front_definition->treatAs2D();

  if (_treat_as_2d)
  {
    //Loop over all elements and pick the ones within radius _r_max from the crack front node
    MeshBase::const_element_iterator el = _mesh.getMesh().elements_begin();
    const MeshBase::const_element_iterator end_el = _mesh.getMesh().elements_end();
    for ( ; el != end_el ; ++el)
    {
      const Elem * elem = *el;
      Point el_pos = elem->centroid();
      Real r;
      Real theta;
      _crack_front_definition->calculateRThetaToCrackFront(el_pos, _crack_front_point_index, r, theta);

      if (r < _r_max)
        _intersected_elems.push_back(const_cast<Elem *>(elem));
    }
  }
  else
  {  
    //Get the midpoint between this crack front point and the next and the tangent at the midpoint
    const Point * p0 = _crack_front_definition->getCrackFrontPoint(_crack_front_point_index);
    RealVectorValue normal0 = _crack_front_definition->getCrackFrontTangent(_crack_front_point_index);

    const Point * p1 = _crack_front_definition->getCrackFrontPoint(_crack_front_point_index+1);
    RealVectorValue normal1 = _crack_front_definition->getCrackFrontTangent(_crack_front_point_index+1);

    const Point pmid = (*p0 + *p1)/2;
    RealVectorValue normalmid = (normal0 + normal1)/2;

    Moose::elementsIntersectedByPlane(pmid, normalmid, _mesh.getMesh(), _intersected_elems);
  }

  _max_weibull_stress = 0.0;
}

void
WeibullStressAtCrackFrontEdge::execute()
{
  if (std::find(_intersected_elems.begin(), _intersected_elems.end(), _current_elem) != _intersected_elems.end())
    _integral_value += computeIntegral();
}

Real
WeibullStressAtCrackFrontEdge::computeQpIntegral()
{
  Real value = 0;
  Real r;
  Real theta;
  _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], _crack_front_point_index, r, theta);

  if (r < _r_max)
  {
    const SymmTensor & tensor(_stress_tensor[_qp]);
    RealVectorValue direction;
    Real max_principal_stress = getTensorQuantity(tensor, _q_point[_qp], direction);
    if (max_principal_stress > _cutoff)
      value = std::pow(max_principal_stress, _m);
  }

  Real edge_length = 1.0;
  if (!_treat_as_2d)
    edge_length = _crack_front_definition->getCrackFrontForwardSegmentLength(_crack_front_point_index);

  return value/edge_length;
}

Real
WeibullStressAtCrackFrontEdge::getValue()
{
  gatherSum(_integral_value);

  Real value = std::pow(_integral_value, 1.0 / _m);
  if (_has_symmetry_plane)
    value *= std::pow(2.0, 1.0 / _m);

  //Weibull stress cannot decrease
  if (value < _max_weibull_stress)
    value = _max_weibull_stress;
  else if (value > _max_weibull_stress)
    _max_weibull_stress = value;

  return value;
}
