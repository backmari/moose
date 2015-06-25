/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef WEIBULLSTRESS_H
#define WEIBULLSTRESS_H

#include "ElementIntegralPostprocessor.h"
#include "MaterialTensorCalculator.h"

//Forward Declarations
class WeibullStress;

template<>
InputParameters validParams<WeibullStress>();

/**
 * This postprocessor computes the Weibull stress
 *
 */
class WeibullStress:
  public ElementIntegralPostprocessor,
  public MaterialTensorCalculator
{
public:
  WeibullStress(const std::string & name, InputParameters parameters);
  virtual Real getValue();

protected:
  virtual void initialSetup();
  virtual Real computeQpIntegral();
  const MaterialProperty<SymmTensor> & _stress_tensor;
  Real _m;
  Real _lambda;
  Real _yield_stress;
  bool _has_symmetry_plane;
  Real _cutoff;
};

#endif //WEIBULLSTRESS_H
