/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "WeibullStressAtCrackFrontEdge.h"
#include "SymmTensor.h"
#include "FEProblem.h"
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
  params.addParam<Real>("r_max","Max radius for Weibull stress calculation");
  params.addParam<unsigned int>("symmetry_plane", "Account for a symmetry plane passing through the plane of the crack, normal to the specified axis (0=x, 1=y, 2=z)");
  params.set<MooseEnum>("quantity") = "MaxPrincipal";
  params.set<MultiMooseEnum>("execute_on") = "timestep";

  return params;
}

WeibullStressAtCrackFrontEdge :: WeibullStressAtCrackFrontEdge(const std::string & name, InputParameters parameters) :
  ElementIntegralPostprocessor(name, parameters),
  MaterialTensorCalculator(name, parameters),
  _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
  _has_crack_front_point_index(isParamValid("crack_front_point_index")),
  _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
  _stress_tensor(getMaterialProperty<SymmTensor>("stress")),
  _m(getParam<Real>("m")),
  _lambda(getParam<Real>("lambda")),
  _yield_stress(getParam<Real>("yield_stress")),
  _r_max(getParam<Real>("r_max")),
  _has_symmetry_plane(isParamValid("symmetry_plane"))
{
  _cutoff = _lambda * _yield_stress;
}

void
WeibullStressAtCrackFrontEdge::initialize()
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

void
WeibullStressAtCrackFrontEdge::execute()
{
  if (std::find(_intersected_elems.begin(), _intersected_elems.end(), _current_elem) != _intersected_elems.end())
    _integral_value += computeIntegral();
}

Real
WeibullStressAtCrackFrontEdge::computeQpIntegral()
{
  Real sum = 0;
  Real r;
  Real theta;
  _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp],_crack_front_point_index,r,theta);

  if (r < _r_max)
  {
    const SymmTensor & tensor(_stress_tensor[_qp]);
    RealVectorValue direction;
    Real principal_stress = getTensorQuantity(tensor, &_q_point[_qp], direction);
    if (principal_stress > _cutoff)
      sum = std::pow(principal_stress,_m);
  }

  Real edge_length = _crack_front_definition->getCrackFrontForwardSegmentLength(_crack_front_point_index);
  return sum/edge_length;
}

Real
WeibullStressAtCrackFrontEdge::getValue()
{
  gatherSum(_integral_value);
  if (_has_symmetry_plane)
    _integral_value *= 2.0;

  return std::pow(_integral_value,1.0/_m);
}
