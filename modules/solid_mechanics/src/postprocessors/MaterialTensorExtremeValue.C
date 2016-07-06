/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "MaterialTensorExtremeValue.h"

#include "SymmTensor.h"
#include "MaterialTensorCalculator.h"

#include <algorithm>
#include <limits>

template<>
InputParameters validParams<MaterialTensorExtremeValue>()
{
  // Define the min/max enumeration
  MooseEnum type_options("max=0 min=1", "max");

  InputParameters params = validParams<ElementExtremeValue>();
  params += validParams<MaterialTensorCalculator>();
  params.addRequiredParam<std::string>("tensor", "The material tensor name.");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

MaterialTensorExtremeValue::MaterialTensorExtremeValue(const InputParameters & parameters) :
    ElementExtremeValue(parameters),
//    _type((ExtremeType)(int)parameters.get<MooseEnum>("value_type")),
//    _value(_type == 0 ? -std::numeric_limits<Real>::max() : std::numeric_limits<Real>::max()),
    _material_tensor_calculator(parameters),
    _tensor(getMaterialProperty<SymmTensor>(getParam<std::string>("tensor")))
{}

void
MaterialTensorExtremeValue::computeQpValue()
{
  RealVectorValue direction;
  Real prop_value = _material_tensor_calculator.getTensorQuantity(_tensor[_qp],
                                                       _q_point[_qp],
                                                       direction);
  switch (_type)
  {
    case MAX:
      _value = std::max(_value, prop_value);
      break;

    case MIN:
      _value = std::min(_value, prop_value);
      break;
  }
}
