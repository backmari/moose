/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef WEIBULLSTRESSFROMSFIS_H
#define WEIBULLSTRESSFROMSFIS_H

#include "ElementIntegralPostprocessor.h"
#include "CrackFrontDefinition.h"

//Forward Declarations
class WeibullStressFromSFIs;

template<>
InputParameters validParams<WeibullStressFromSFIs>();

/**
 * This postprocessor computes the Weibull stress
 *
 */
class WeibullStressFromSFIs:
  public ElementIntegralPostprocessor
{
public:
  WeibullStressFromSFIs(const InputParameters & parameters);
  virtual Real getValue();

protected:
  virtual void initialSetup();
  virtual Real computeQpIntegral();
  Real computePrincipalStress(PostprocessorValue ki, PostprocessorValue kii, PostprocessorValue kiii);
  VariableValue & _scalar_q;
  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  const PostprocessorValue & _ki_value;
  const PostprocessorValue & _kii_value;
  const PostprocessorValue & _kiii_value;
  Real _poissons_ratio;
  Real _m;
  Real _lambda;
  Real _yield_stress;
  bool _has_symmetry_plane;
  Real _cutoff;
  Real _r;
  Real _theta;
};

#endif //WEIBULLSTRESSFROMSFIS_H
