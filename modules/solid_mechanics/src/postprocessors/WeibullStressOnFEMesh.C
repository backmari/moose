/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "WeibullStressOnFEMesh.h"
#include "SymmTensor.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "PlaneTracing.h"
#include <cmath>
#include <algorithm>
#include <set>

template<>
InputParameters validParams<WeibullStressOnFEMesh>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params += validParams<MaterialTensorCalculator>();

  params.addRequiredParam<UserObjectName>("crack_front_definition", "The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index", "The index of the point on the crack front");
  params.addRequiredParam<PostprocessorName>("KI_name","The name of the KI postprocessor");
  params.addRequiredParam<PostprocessorName>("KII_name","The name of the KII postprocessor");
  params.addRequiredParam<PostprocessorName>("KIII_name","The name of the KIII postprocessor");
  params.addRequiredParam<Real>("poissons_ratio","Poisson's ratio for the material.");
  params.addRequiredParam<Real>("yield_stress", "Yield stress of the material");
  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addParam<Real>("weibull_r_max","Max radius for Weibull stress calculation");
  params.addParam<unsigned int>("symmetry_plane", "Account for a symmetry plane passing through the plane of the crack, normal to the specified axis (0=x, 1=y, 2=z)");
  params.addParam<Real>("crack_tip_radius","Radius of the blunt crack tip.");
  params.set<MooseEnum>("quantity") = "MaxPrincipal";
  params.set<MultiMooseEnum>("execute_on") = "timestep_end";

  return params;
}

WeibullStressOnFEMesh :: WeibullStressOnFEMesh(const InputParameters & parameters) :
    ElementIntegralPostprocessor(parameters),
    MaterialTensorCalculator(parameters),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _has_crack_front_point_index(isParamValid("crack_front_point_index")),
    _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    _ki_value(getPostprocessorValue("KI_name")),
    _kii_value(getPostprocessorValue("KII_name")),
    _kiii_value(getPostprocessorValue("KIII_name")),
    _stress_tensor(getMaterialProperty<SymmTensor> ("stress")),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _m(getParam<Real>("m")),
    _lambda(getParam<Real>("lambda")),
    _yield_stress(getParam<Real>("yield_stress")),
    _r_max(getParam<Real>("weibull_r_max")),
    _has_symmetry_plane(isParamValid("symmetry_plane"))
{
  _cutoff = _lambda * _yield_stress;

  if (isParamValid("crack_tip_radius"))
    _rho = 2 * getParam<Real>("crack_tip_radius");
}

void
WeibullStressOnFEMesh::initialSetup()
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
}

void
WeibullStressOnFEMesh::execute()
{
  if (std::find(_intersected_elems.begin(), _intersected_elems.end(), _current_elem) != _intersected_elems.end())
    _integral_value += computeIntegral();
}

Real
WeibullStressOnFEMesh::computeQpIntegral()
{
  Real value = 0;
  Real r;
  Real theta;

  Point p(_q_point[_qp]);
  p(0) = p(0) + 0.5*_rho;
  _crack_front_definition->calculateRThetaToCrackFront(p, _crack_front_point_index, r, theta);

  if (r < _r_max)
  {
    const SymmTensor & tensor(_stress_tensor[_qp]);
    RealVectorValue direction;
    Real principal_stress1 = getTensorQuantity(tensor, &_q_point[_qp], direction);
    
    Real principal_stress2 = computePrincipalStress(r, theta);

//    if (principal_stress1 > _cutoff || principal_stress2 > _cutoff)
//      Moose::out<<"DEBUG2 "<<(principal_stress1 - principal_stress2)/principal_stress1<<std::endl;
    
    if (principal_stress2 > _cutoff && principal_stress2 < 3*_yield_stress)
      value = std::pow(principal_stress2,_m);
  }

  Real edge_length = 1.0;
  if (!_treat_as_2d)
    _crack_front_definition->getCrackFrontForwardSegmentLength(_crack_front_point_index);

  return value/edge_length;
}

Real
WeibullStressOnFEMesh::getValue()
{
  gatherSum(_integral_value);
  return std::pow(_integral_value,1.0/_m);
}


Real
WeibullStressOnFEMesh::computePrincipalStress(Real r, Real theta)
{
  Real st2 = std::sin(theta / 2);
  Real stt2 = std::sin(3 * theta / 2);
  Real ct2 = std::cos(theta / 2);
  Real ctt2 = std::cos(3 * theta / 2);
  Real ki2PiR =   _ki_value   / std::sqrt(2 * libMesh::pi * r);
  Real kii2PiR =  _kii_value  / std::sqrt(2 * libMesh::pi * r);
  Real kiii2PiR = _kiii_value / std::sqrt(2 * libMesh::pi * r);  
  Real rho2R = _rho / (2 * r);

  ColumnMajorMatrix stress(3,3);
  stress(0,0) = ki2PiR * ct2 * (1 - st2 * stt2) - ki2PiR * rho2R * ctt2 - kii2PiR * st2 * (2 + ct2 * ctt2) + kii2PiR * rho2R * stt2;
  stress(1,1) = ki2PiR * ct2 * (1 + st2 * stt2) + ki2PiR * rho2R * ctt2 + kii2PiR * st2 * ct2 * ctt2 - kii2PiR * rho2R * stt2;
  stress(0,1) = ki2PiR * st2 * ct2 * ctt2 - ki2PiR * rho2R * stt2 + kii2PiR * ct2 * (1 - st2 * stt2) - kii2PiR * rho2R * ct2;
  stress(1,0) = stress(0,1);
  stress(0,2) = -kiii2PiR * st2;
  stress(2,0) = stress(0,2);
  stress(1,2) = kiii2PiR * ct2;
  stress(2,1) = stress(1,2);
  stress(2,2) = _poissons_ratio * (2 * ki2PiR * ct2 - 2 * kii2PiR * st2);

  ColumnMajorMatrix eval(3,1);
  ColumnMajorMatrix evec(3,3);
  stress.eigen(eval,evec);

  return eval(2);
}
