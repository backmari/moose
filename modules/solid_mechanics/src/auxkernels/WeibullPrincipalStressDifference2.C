/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "WeibullPrincipalStressDifference2.h"
#include "PlaneTracing.h"
#include "SymmTensor.h"
#include "FEProblem.h"

#include <algorithm>

template<>
InputParameters validParams<WeibullPrincipalStressDifference2>()
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

WeibullPrincipalStressDifference2::WeibullPrincipalStressDifference2(const InputParameters & parameters) :
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
  if ( _nodal ) mooseError("From WeibullPrincipalStressDifference2: element on plane id must be an element property (use CONSTANT MONOMIAL)");

  if (isParamValid("crack_tip_radius"))
    _rho = 2 * getParam<Real>("crack_tip_radius");

  if (_rho > 0 && _crack_tip_shape == "Sharp")
    mooseError("A sharp crack tip cannot have a non-zero crack tip radius");

  _cutoff = _lambda * _yield_stress;
}

void
WeibullPrincipalStressDifference2::initialSetup()
{
}

Real
WeibullPrincipalStressDifference2::computeValue()
{
  Real value = 0;
  Real r;
  Real theta;

  Point p(_q_point[_qp]);
  p(0) = p(0) + 0.5*_rho;
  _crack_front_definition->calculateRThetaToCrackFront(p, _crack_front_point_index, r, theta);

  const SymmTensor & tensor(_stress_tensor[_qp]);
  RealVectorValue direction;
    
  Real principal_stress1 = getTensorQuantity(tensor, p, direction);
  Real principal_stress2 = computePrincipalStress(r, theta);

  value = (principal_stress2 - principal_stress1)/ principal_stress2;
  
  return value;
}

Real
WeibullPrincipalStressDifference2::computePrincipalStress(Real r, Real theta)
{
  Real st2 = std::sin(theta / 2);
  Real st32 = std::sin(3 * theta / 2);
  Real st52 = std::sin(5 * theta / 2);
  Real st72 = std::sin(7 * theta / 2);
  Real ct = std::cos(theta);
  Real ct2 = std::cos(theta / 2);
  Real ct32 = std::cos(3 * theta / 2);
  Real ct52 = std::cos(5 * theta / 2);
  Real ct72 = std::cos(7 * theta / 2);
  Real ki2PiR =   _ki_value   / std::sqrt(2 * libMesh::pi * r);
  Real kii2PiR =  _kii_value  / std::sqrt(2 * libMesh::pi * r);
  Real kiii2PiR = _kiii_value / std::sqrt(2 * libMesh::pi * r);
  Real rhoR = _rho / r;
  Real rhoRsq = std::pow(_rho / r, 2);
  Real rhoRcu = std::pow(_rho / r, 3);

  ColumnMajorMatrix stress(3,3);

  stress(0,0) = ki2PiR * (0.75 * ct2 + (0.25 - 0.75 * rhoRsq) * ct52 + (3.0 / 8.0 * rhoR - 5.0 / 8.0 * rhoRcu) * ct72)
    + kii2PiR * (-7.0 / 4.0 * st2 - 2 * rhoR * st32 - (0.25 - 0.75 * rhoRsq) * st52 + (-9.0 / 8.0 * rhoR + 15.0 / 8.0 * rhoRcu) * st72);

  stress(1,1) = ki2PiR * (1.25 * ct2 + rhoR * ct32 + (-0.25 + 0.75 * rhoRsq) * ct52 + (-3.0 / 8.0 * rhoR + 5.0 / 8.0 * rhoRcu) * ct72)
    + kii2PiR * (-0.25 * st2 - rhoR * st32 + (0.25 - 0.75 * rhoRsq) * st52 + (9.0 / 8.0 * rhoR - 15.0 / 8.0 * rhoRcu) * st72);
  
  stress(0,1) = ki2PiR * (-0.25 * st2 - 0.5 * rhoR * st32 + (0.25 - 0.75 * rhoRsq) * st52 + (3.0 / 8.0 * rhoR - 5.0 / 8.0 * rhoRcu) * st72)
    + kii2PiR * (0.75 * ct2 + 0.5 * rhoR * ct32 + (0.25 - 0.75 * rhoRsq) * ct52 + (9.0 / 8.0 * rhoR - 15.0 / 8.0 * rhoRcu) * ct72);

  stress(1,0) = stress(0,1);

  stress(1,2) = kiii2PiR * (ct2 - rhoR * ct2 + 2 * rhoR * ct2 * ct);

  stress(2,1) = stress(1,2);

  stress(0,2) = -kiii2PiR * (st2 + rhoR * st2 + 2 * rhoR * st2 * ct);

  stress(2,0) = stress(0,2);

  stress(2,2) = _poissons_ratio * (ki2PiR * (2 * ct2 + rhoR * ct32) - kii2PiR * (2 * st2 + 3 * rhoR * st32));

  ColumnMajorMatrix eval(3,1);
  ColumnMajorMatrix evec(3,3);
  stress.eigen(eval,evec);

  return eval(2);
}
