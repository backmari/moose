/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ElementsContribToWeibull.h"
#include "PlaneTracing.h"
#include "MooseMesh.h"
#include <algorithm>

template<>
InputParameters validParams<ElementsContribToWeibull>()
{
  InputParameters params = validParams<AuxKernel>();
  params += validParams<MaterialTensorCalculator>();

  params.addRequiredParam<UserObjectName>("crack_front_definition", "The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index", "The index of the point on the crack front");
  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addParam<Real>("weibull_r_max","Max radius for Weibull stress calculation");
  params.addRequiredParam<Real>("yield_stress","Yield stress of the material");
  params.set<MooseEnum>("quantity") = "MaxPrincipal";
  
  params.set<MultiMooseEnum>("execute_on") = "timestep_end";

  return params;
}

ElementsContribToWeibull::ElementsContribToWeibull(const InputParameters & parameters):
    AuxKernel(parameters),
    MaterialTensorCalculator(parameters),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _has_crack_front_point_index(isParamValid("crack_front_point_index")),
    _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    _stress_tensor(getMaterialProperty<SymmTensor>("stress")),
    _m(getParam<Real>("m")),
    _lambda(getParam<Real>("lambda")),
    _yield_stress(getParam<Real>("yield_stress")),
    _r_max(getParam<Real>("weibull_r_max"))
{
  _cutoff = _lambda * _yield_stress;
}

void
ElementsContribToWeibull::initialSetup()
{
  _treat_as_2d = _crack_front_definition->treatAs2D();
  _crack_front_length = 1.0;
  if (!_treat_as_2d)
    _crack_front_length = _crack_front_definition->getCrackFrontLength(); 
}

Real
ElementsContribToWeibull::computeValue()
{
  Real value = 0;
  Real r;
  Real theta;
  _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], _crack_front_point_index, r, theta);

  if (r < _r_max)
  {
    const SymmTensor & tensor(_stress_tensor[_qp]);
    RealVectorValue direction;
    Real principal_stress = getTensorQuantity(tensor, &_q_point[_qp], direction);
    if (principal_stress > _cutoff && principal_stress < 3*_yield_stress)
      value = std::pow(principal_stress,_m);
  }

  return std::pow(value/_crack_front_length, 1.0/_m);
}
