/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef WEIBULLPRINCIPALSTRESSDIFFERENCE_H
#define WEIBULLPRINCIPALSTRESSDIFFERENCE_H

#include "AuxKernel.h"
#include "MaterialTensorCalculator.h"
#include "CrackFrontDefinition.h"

class WeibullPrincipalStressDifference :
  public AuxKernel,
  public MaterialTensorCalculator
{

public:

  WeibullPrincipalStressDifference(const InputParameters & parameters);

  virtual ~WeibullPrincipalStressDifference() {}

protected:

  virtual void initialSetup();
  virtual Real computeValue();
  Real computePrincipalStress(Real r, Real theta);
  Real computePrincipalStress2(Real r, Real theta);

  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  const PostprocessorValue & _ki_value;
  const PostprocessorValue & _kii_value;
  const PostprocessorValue & _kiii_value;
  MooseEnum _crack_tip_shape;
  const MaterialProperty<SymmTensor> & _stress_tensor;
  Real _poissons_ratio;
  Real _m;
  Real _lambda;
  Real _yield_stress;
  Real _r_max;
  Real _cutoff;
  Real _rho;

private:

  std::vector<Elem *> _intersected_elems;
};

template<>
InputParameters validParams<WeibullPrincipalStressDifference>();

#endif // ELEMENTSINTERSECTEDBYPLANE_H
