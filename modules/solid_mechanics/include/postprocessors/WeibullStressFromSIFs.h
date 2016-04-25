/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef WEIBULLSTRESSFROMSIFS_H
#define WEIBULLSTRESSFROMSIFS_H

#include "GeneralPostprocessor.h"
#include "CrackFrontDefinition.h"

//Forward Declarations
class WeibullStressFromSIFs;

template<>
InputParameters validParams<WeibullStressFromSIFs>();

/**
 * This postprocessor computes the Weibull stress
 *
 */
class WeibullStressFromSIFs:
  public GeneralPostprocessor
{
public:
  WeibullStressFromSIFs(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual Real getValue();

protected:

  Real computeWeibullStress(PostprocessorValue ki, PostprocessorValue kii, PostprocessorValue kiii);
  Real computePrincipalStress(Real r, Real theta);
  void generateRegularMesh();
  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  const PostprocessorValue & _ki_value;
  const PostprocessorValue & _kii_value;
  const PostprocessorValue & _kiii_value;
  MooseEnum _crack_tip_shape;
  Real _poissons_ratio;
  Real _m;
  Real _lambda;
  Real _yield_stress;
  Real _r_max;
  Real _cutoff;
  Real _rho;
  std::vector<Real> _r;
  std::vector<Real> _theta;
  unsigned int _npoints;
  Real _dr;
  Real _dtheta;
};

#endif //WEIBULLSTRESSFROMSIFS_H
