/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This post processor calculates the J-Integral
//
#include "WeibullStress.h"

template<>
InputParameters validParams<WeibullStress>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params += validParams<MaterialTensorCalculator>();

  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addRequiredParam<Real>("yield_stress","Yield stress of the material");
  params.set<MooseEnum>("quantity") = "MaxPrincipal";
  params.addParam<unsigned int>("symmetry_plane", "Account for a symmetry plane passing through the plane of the crack, normal to the specified axis (0=x, 1=y, 2=z)");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

WeibullStress::WeibullStress(const std::string & name, InputParameters parameters):
    ElementIntegralPostprocessor(name, parameters),
    MaterialTensorCalculator(name, parameters),
    _stress_tensor(getMaterialProperty<SymmTensor>("stress")),
    _m(getParam<Real>("m")),
    _lambda(getParam<Real>("lambda")),
    _yield_stress(getParam<Real>("yield_stress")),
    _has_symmetry_plane(isParamValid("symmetry_plane"))
{
  _cutoff = _lambda * _yield_stress;
}

void
WeibullStress::initialSetup()
{
}

Real
WeibullStress::computeQpIntegral()
{
  const SymmTensor & tensor(_stress_tensor[_qp]);
  RealVectorValue direction;

  Real principal_stress = getTensorQuantity(tensor,&_q_point[_qp],direction);
  Real value;
  if (principal_stress > _cutoff)
    value = std::pow(principal_stress,_m);

  return value;
}

Real
WeibullStress::getValue()
{
  gatherSum(_integral_value);
  if (_has_symmetry_plane)
    _integral_value *= 2.0;

  return std::pow(_integral_value,1.0/_m);
}
