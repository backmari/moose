/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This post processor calculates the Weibull stress
//
#include "WeibullStressFromSFIs.h"

template<>
InputParameters validParams<WeibullStressFromSFIs>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();

  params.addCoupledVar("q", "The q function, aux variable");
  params.addRequiredParam<UserObjectName>("crack_front_definition","The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index","The index of the point on the crack front corresponding to this q function");
  params.addRequiredParam<PostprocessorName>("KI_name","The name of the KI postprocessor");
  params.addRequiredParam<PostprocessorName>("KII_name","The name of the KII postprocessor");
  params.addRequiredParam<PostprocessorName>("KIII_name","The name of the KIII postprocessor");
  params.addRequiredParam<Real>("poissons_ratio","Poisson's ratio for the material.");
  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addRequiredParam<Real>("yield_stress","Yield stress of the material");
  params.addParam<unsigned int>("symmetry_plane", "Account for a symmetry plane passing through the plane of the crack, normal to the specified axis (0=x, 1=y, 2=z)");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

WeibullStressFromSFIs::WeibullStressFromSFIs(const InputParameters & parameters):
    ElementIntegralPostprocessor(parameters),
    _scalar_q(coupledValue("q")),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _has_crack_front_point_index(isParamValid("crack_front_point_index")),
    _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    _ki_value(getPostprocessorValue("KI_name")),
    _kii_value(getPostprocessorValue("KII_name")),
    _kiii_value(getPostprocessorValue("KIII_name")),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _m(getParam<Real>("m")),
    _lambda(getParam<Real>("lambda")),
    _yield_stress(getParam<Real>("yield_stress")),
    _has_symmetry_plane(isParamValid("symmetry_plane"))
{
  _cutoff = _lambda * _yield_stress;
}

void
WeibullStressFromSFIs::initialSetup()
{
}

Real
WeibullStressFromSFIs::computeQpIntegral()
{
  Moose::out<<"now in WeibullStressFromSFIs::computeQpIntegral()"<<std::endl;
  
  Point p(_q_point[_qp]);
  _crack_front_definition->calculateRThetaToCrackFront(p,_crack_front_point_index,_r,_theta);

  Real principal_stress = computePrincipalStress(_ki_value,_kii_value,_kiii_value);
  //  Moose::out<<"node "<<_crack_front_point_index<<" element "<<_current_elem->id()<<" "<<_qp<<" stress "<<principal_stress<<std::endl;
//  Moose::out<<"node "<<_crack_front_point_index<<" element "<<_current_elem->id()<<" "<<_qp<<" r "<<_r<<" theta "<<_theta<<std::endl;
  Real value(0.0);
  if (principal_stress > _cutoff && _scalar_q[_qp] > 0.0)
  {
    value = std::pow(principal_stress,_m);
  }

  return value;
}

Real
WeibullStressFromSFIs::getValue()
{
  Moose::out<<"now in WeibullStressFromSFIs::getValue()"<<std::endl;

  gatherSum(_integral_value);
  if (_has_symmetry_plane)
    _integral_value *= 2.0;

  return std::pow(_integral_value,1.0/_m);
}

Real
WeibullStressFromSFIs::computePrincipalStress(const PostprocessorValue _ki, const PostprocessorValue _kii, const PostprocessorValue _kiii)
{
  Real sqrt2PiR = std::sqrt(2*libMesh::pi*_r);
  Real t = _theta;
  Real t2 = _theta/2;
  Real tt2 = 3*_theta/2;
  Real st2 = std::sin(t2);
  Real ct2 = std::cos(t2);
  Real stt2 = std::sin(tt2);
  Real ctt2 = std::cos(tt2);

  ColumnMajorMatrix stress(3,3);
  stress(0,0) = 1 / sqrt2PiR * ( _ki * ct2 * (1 - st2 * stt2)
                                  - _kii * st2 * (2 + ct2 * ctt2));
  stress(1,1) = 1 / sqrt2PiR * ( _ki * ct2 * (1 + st2 * stt2)
                                  + _kii * st2 * ct2 * ctt2);
  stress(0,1) = 1 / sqrt2PiR * ( _ki * ct2 * st2 * ctt2
                                  + _kii * ct2 * (1 - st2*stt2));
  stress(1,0) = stress(0,1);
  stress(0,2) = -1 / sqrt2PiR * _kiii * st2;
  stress(2,0) = stress(0,2);
  stress(1,2) =  1 / sqrt2PiR * _kiii * ct2;
  stress(2,1) = stress(1,2);
  stress(2,2) = _poissons_ratio * (stress(0,0) + stress(1,1));

//  stress.print();
  
  ColumnMajorMatrix eval(3,1);
  ColumnMajorMatrix evec(3,3);
  stress.eigen(eval,evec);

  return eval(0);
}
