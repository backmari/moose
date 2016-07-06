/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This post processor calculates the J-Integral
//
#include "WeibullStress.h"
#include "PiecewiseLinear.h"
#include "MaterialTensorCalculatorTools.h"

// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<WeibullStress>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params += validParams<MaterialTensorCalculator>();

  params.addRequiredParam<UserObjectName>("crack_front_definition", "The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index", "The index of the point on the crack front");
  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addParam<Real>("weibull_r_max","Max radius for Weibull stress calculation");
  params.addRequiredParam<Real>("yield_stress","Yield stress of the material");
  params.set<MooseEnum>("quantity") = "MaxPrincipal";
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

WeibullStress::WeibullStress(const InputParameters & parameters):
    ElementIntegralPostprocessor(parameters),
    MaterialTensorCalculator(parameters),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _has_crack_front_point_index(isParamValid("crack_front_point_index")),
    _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    //    _von_mises_calculator(NULL),
    //    _max_princ_calculator(NULL),
    _stress_tensor(getMaterialProperty<SymmTensor>("stress")),
    _m(getParam<Real>("m")),
    _lambda(getParam<Real>("lambda")),
    _yield_stress(getParam<Real>("yield_stress")),
    _r_max(getParam<Real>("weibull_r_max")),
    _has_symmetry_plane(isParamValid("symmetry_plane"))
{
  _cutoff = _lambda * _yield_stress;
//  _max_cutoff = 3 * _yield_stress;

  //  InputParameters von_mises_params = emptyInputParameters();
  //  von_mises_params += this->parameters();
  //  von_mises_params.set<MooseEnum>("quantity") = "VonMises";
  //  _von_mises_calculator = new MaterialTensorCalculator::MaterialTensorCalculator(von_mises_params);

  //  InputParameters max_princ_params = emptyInputParameters();
  //  max_princ_params += this->parameters();
  //  max_princ_params.set<MooseEnum>("quantity") = "MaxPrincipal";
  //  _max_princ_calculator = new MaterialTensorCalculator::MaterialTensorCalculator(max_princ_params);
}

//WeibullStress::~WeibullStress()
//{
//  delete _von_mises_calculator;
//  delete _max_princ_calculator;
//}

void
WeibullStress::initialSetup()
{
  _treat_as_2d = _crack_front_definition->treatAs2D();
}

Real
WeibullStress::computeQpIntegral()
{
  Real value = 0;
  Real r;
  Real theta;
  _crack_front_definition->calculateRThetaToCrackFront(_q_point[_qp], _crack_front_point_index, r, theta);

  if (r < _r_max)
  {
    const SymmTensor & tensor(_stress_tensor[_qp]);
    RealVectorValue direction;
    //    Real von_mises = MaterialTensorCalculatorTools::vonMisesStress(tensor);
    Real max_principal_stress = MaterialTensorCalculatorTools::maxPrinciple(tensor, direction);
    if (max_principal_stress > _cutoff)
      value = std::pow(max_principal_stress, _m);
  }

  return value;
};

// Real
// WeibullStress::computeQpIntegral()
// {
//   mooseError("From WeibullStress: computeQpIntegral() is not defined");
//   return 0;
// }

// Real
// WeibullStress::computeIntegral()
// {
//   Real value = 0.0;
//   Real principal_stress = 0.0;
//   RealVectorValue direction;

//   Moose::out<<"qp ";
//   for (_qp=0; _qp<_qrule->n_points(); _qp++)
//   {
//     const SymmTensor & tensor(_stress_tensor[_qp]);
//     principal_stress += _JxW[_qp]*_coord[_qp]*getTensorQuantity(tensor, &_q_point[_qp], direction);
//     Moose::out<<getTensorQuantity(tensor, &_q_point[_qp], direction)<<" "<<_JxW[_qp]*_coord[_qp]<<" qp ";
//   }
//   Moose::out<<std::endl;
//   //Is this correct for all the possible cases of element types?
//   principal_stress /= _current_elem_volume;
//   Moose::out<<"elem "<<principal_stress<<std::endl;

//   if (principal_stress > _max_cutoff  && principal_stress < 3 * _yield_stress)
//     value = std::pow(principal_stress, _m);

// //  Real r;
// //  Real theta;
// //  _crack_front_definition->calculateRThetaToCrackFront(_current_elem->centroid(), _crack_front_point_index, r, theta);
// //  if (r < 0.5)
// //    value = std::pow(principal_stress, _m);

//   return value;
// }

Real
WeibullStress::getValue()
{
  gatherSum(_integral_value);

  Real crack_front_length = 1.0;
  if (!_treat_as_2d)
    crack_front_length = _crack_front_definition->getCrackFrontLength();    

  Moose::out<<"result "<<std::pow(_integral_value/crack_front_length, 1.0/_m)<<" "<<_integral_value<<" "<<crack_front_length<<" "<<_m<<" "<<_treat_as_2d<<std::endl;

  if (std::isnan(_integral_value) == true)
    Moose::out<<"NAN "<<_integral_value<<std::endl;

  return std::pow(_integral_value/crack_front_length, 1.0/_m);
}
