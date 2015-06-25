/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef OLDWEIBULLSTRESS_H
#define OLDWEIBULLSTRESS_H

#include "ElementIntegralPostprocessor.h"
#include "MaterialTensorCalculator.h"

//Forward Declarations
class OldWeibullStress;

template<>
InputParameters validParams<OldWeibullStress>();

/**
 * This postprocessor computes the OldWeibull stress
 *
 */
class OldWeibullStress:
  public ElementIntegralPostprocessor,
  public MaterialTensorCalculator
{
public:
  OldWeibullStress(const InputParameters & parameters);
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

#endif //OLDWEIBULLSTRESS_H
