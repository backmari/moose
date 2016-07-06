/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This post processor calculates the Weibull stress
//
#include "WeibullStressAtCrackFrontNode.h"

template<>
InputParameters validParams<WeibullStressAtCrackFrontNode>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params += validParams<MaterialTensorCalculator>();

  params.addCoupledVar("q", "The q function, aux variable");
  params.addRequiredParam<UserObjectName>("crack_front_definition","The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index","The index of the point on the crack front corresponding to this q function");
  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addRequiredParam<Real>("yield_stress","Yield stress of the material");
  params.set<MooseEnum>("quantity") = "MaxPrincipal";
  params.addParam<unsigned int>("symmetry_plane", "Account for a symmetry plane passing through the plane of the crack, normal to the specified axis (0=x, 1=y, 2=z)");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

WeibullStressAtCrackFrontNode::WeibullStressAtCrackFrontNode(const InputParameters & parameters):
    ElementIntegralPostprocessor(parameters),
    MaterialTensorCalculator(parameters),
    _scalar_q(coupledValue("q")),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _has_crack_front_point_index(isParamValid("crack_front_point_index")),
    _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    _stress_tensor(getMaterialProperty<SymmTensor>("stress")),
    _m(getParam<Real>("m")),
    _lambda(getParam<Real>("lambda")),
    _yield_stress(getParam<Real>("yield_stress")),
    _has_symmetry_plane(isParamValid("symmetry_plane"))
{
  _cutoff = _lambda * _yield_stress;
}

void
WeibullStressAtCrackFrontNode::initialSetup()
{
}

Real
WeibullStressAtCrackFrontNode::computeQpIntegral()
{
  const SymmTensor & tensor(_stress_tensor[_qp]);
  RealVectorValue direction;

  Real principal_stress = getTensorQuantity(tensor, _q_point[_qp], direction);
  Real value(0.0);
  if (principal_stress > _cutoff && _scalar_q[_qp] > 0.0)
  {
    value = std::pow(principal_stress,_m);
//    Moose::out<<"node "<<_crack_front_point_index<<" element "<<_current_elem->id()<<" "<<_qp<<" stress "<<principal_stress<<std::endl;
  }

  return value;
}

Real
WeibullStressAtCrackFrontNode::getValue()
{
  gatherSum(_integral_value);
  if (_has_symmetry_plane)
    _integral_value *= 2.0;

  return std::pow(_integral_value,1.0/_m);
}
