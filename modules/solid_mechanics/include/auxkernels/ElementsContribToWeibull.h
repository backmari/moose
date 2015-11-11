/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELEMENTSCONTRIBTOWEIBULL_H
#define ELEMENTSCONTRIBTOWEIBULL_H

#include "AuxKernel.h"
#include "MaterialTensorCalculator.h"
#include "CrackFrontDefinition.h"

class ElementsContribToWeibull :
  public AuxKernel,
  public MaterialTensorCalculator
{

public:

  ElementsContribToWeibull(const InputParameters & parameters);

  virtual ~ElementsContribToWeibull() {}

protected:

  virtual void initialSetup();
//  virtual void compute();
  virtual Real computeValue();
  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  const MaterialProperty<SymmTensor> & _stress_tensor;
  Real _m;
  Real _lambda;
  Real _yield_stress;
  Real _r_max;
  bool _has_symmetry_plane;
  Real _cutoff;
  Real _crack_front_length;

private:
  bool _treat_as_2d;
};

template<>
InputParameters validParams<ElementsContribToWeibull>();

#endif // ELEMENTSCONTRIBTOWEIBULL_H
