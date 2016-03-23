/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "WeibullPrincipalStressDifference.h"
#include "PlaneTracing.h"
#include "SymmTensor.h"
#include "FEProblem.h"

#include <algorithm>

template<>
InputParameters validParams<WeibullPrincipalStressDifference>()
{
  InputParameters params = validParams<AuxKernel>();
  params += validParams<MaterialTensorCalculator>();

  params.addRequiredParam<UserObjectName>("crack_front_definition", "The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index", "The index of the point on the crack front");
  params.addRequiredParam<PostprocessorName>("KI_name","The name of the KI postprocessor");
  params.addRequiredParam<PostprocessorName>("KII_name","The name of the KII postprocessor");
  params.addRequiredParam<PostprocessorName>("KIII_name","The name of the KIII postprocessor");
  MooseEnum tip_shape("Sharp Blunt","Sharp");
  params.addParam<MooseEnum>("crack_tip_shape", tip_shape, "Shape of the crack tip. Choices are: " + tip_shape.getRawNames());
  params.addRequiredParam<Real>("poissons_ratio","Poisson's ratio for the material.");
  params.addRequiredParam<Real>("yield_stress", "Yield stress of the material");
  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addParam<Real>("crack_tip_radius", 0.0, "Radius of the blunt crack tip.");
  params.addParam<Real>("weibull_r_max","Max radius for Weibull stress calculation");
  params.set<MooseEnum>("quantity") = "MaxPrincipal";
  params.set<MultiMooseEnum>("execute_on") = "timestep_end";

  return params;
}

WeibullPrincipalStressDifference::WeibullPrincipalStressDifference(const InputParameters & parameters) :
    AuxKernel(parameters),
    MaterialTensorCalculator(parameters),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _has_crack_front_point_index(isParamValid("crack_front_point_index")),
    _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    _ki_value(getPostprocessorValue("KI_name")),
    _kii_value(getPostprocessorValue("KII_name")),
    _kiii_value(getPostprocessorValue("KIII_name")),
    _crack_tip_shape(getParam<MooseEnum>("crack_tip_shape")),
    _stress_tensor(getMaterialProperty<SymmTensor> ("stress")),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _m(getParam<Real>("m")),
    _lambda(getParam<Real>("lambda")),
    _yield_stress(getParam<Real>("yield_stress")),
    _r_max(getParam<Real>("weibull_r_max"))
{
  if ( _nodal ) mooseError("From WeibullPrincipalStressDifference: element on plane id must be an element property (use CONSTANT MONOMIAL)");

  if (isParamValid("crack_tip_radius"))
    _rho = 2 * getParam<Real>("crack_tip_radius");

  if (_rho > 0 && _crack_tip_shape == "Sharp")
    mooseError("A sharp crack tip cannot have a non-zero crack tip radius");

  _cutoff = _lambda * _yield_stress;
}

void
WeibullPrincipalStressDifference::initialSetup()
{
}

Real
WeibullPrincipalStressDifference::computeValue()
{
  Real value = 0;
  Real r;
  Real theta;

  Point p(_q_point[_qp]);
  p(0) = p(0) + 0.5*_rho;
  _crack_front_definition->calculateRThetaToCrackFront(p, _crack_front_point_index, r, theta);

  const SymmTensor & tensor(_stress_tensor[_qp]);
  RealVectorValue direction;
    
  Real principal_stress1 = getTensorQuantity(tensor, &p, direction);
  Real principal_stress2 = computePrincipalStress(r, theta);

  value = (principal_stress2 - principal_stress1)/ principal_stress2;
  
  return value;
}

Real
WeibullPrincipalStressDifference::computePrincipalStress(Real r, Real theta)
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


Real
WeibullPrincipalStressDifference::computePrincipalStress2(Real r, Real theta)
{
  Real ki2PiR =   _ki_value   / std::sqrt(2 * libMesh::pi * r);
  Real principal_stress = ki2PiR * (std::cos(theta / 2) + 0.5 * std::sqrt(std::pow(_rho/r,2)+std::pow(std::sin(theta),2)));
  return principal_stress;
}
