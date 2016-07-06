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

#ifndef MATERIALTENSOREXTREMEVALUE_H
#define MATERIALTENSOREXTREMEVALUE_H

#include "ElementExtremeValue.h"
#include "MaterialTensorCalculator.h"

class MaterialTensorExtremeValue;
class SymmTensor;

template<>
InputParameters validParams<MaterialTensorExtremeValue>();

class MaterialTensorExtremeValue : public ElementExtremeValue
{
public:
  MaterialTensorExtremeValue(const InputParameters & parameters);

protected:
  virtual void computeQpValue();

  MaterialTensorCalculator _material_tensor_calculator;
  const MaterialProperty<SymmTensor> & _tensor;
};

#endif