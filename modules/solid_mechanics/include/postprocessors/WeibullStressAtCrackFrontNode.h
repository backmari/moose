/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef WEIBULLSTRESSATCRACKFRONTNODE_H
#define WEIBULLSTRESSATCRACKFRONTNODE_H

#include "ElementIntegralPostprocessor.h"
#include "MaterialTensorCalculator.h"
#include "CrackFrontDefinition.h"

//Forward Declarations
class WeibullStressAtCrackFrontNode;

template<>
InputParameters validParams<WeibullStressAtCrackFrontNode>();

/**
 * This postprocessor computes the Weibull stress
 *
 */
class WeibullStressAtCrackFrontNode:
  public ElementIntegralPostprocessor,
  public MaterialTensorCalculator
{
public:
  WeibullStressAtCrackFrontNode(const InputParameters & parameters);
  virtual Real getValue();

protected:
  virtual void initialSetup();
  virtual Real computeQpIntegral();
  const VariableValue & _scalar_q;
  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  const MaterialProperty<SymmTensor> & _stress_tensor;
  Real _m;
  Real _lambda;
  Real _yield_stress;
  bool _has_symmetry_plane;
  Real _cutoff;
};

#endif //WEIBULLSTRESSATCRACKFRONTNODE_H
