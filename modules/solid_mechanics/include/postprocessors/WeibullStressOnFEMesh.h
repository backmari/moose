/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef WEIBULLSTRESSATTEST_H
#define WEIBULLSTRESSATTEST_H

#include "ElementIntegralPostprocessor.h"
#include "MaterialTensorCalculator.h"
#include "CrackFrontDefinition.h"

class WeibullStressOnFEMesh;
class SymmTensor;

template<>
InputParameters validParams<WeibullStressOnFEMesh>();

class WeibullStressOnFEMesh :
  public ElementIntegralPostprocessor,
  public MaterialTensorCalculator
{
public:
  WeibullStressOnFEMesh(const InputParameters & parameters);

  virtual void initialSetup();
  virtual void execute();
//  virtual void threadJoin(const UserObject & u );
  virtual Real getValue();

protected:

  virtual Real computeQpIntegral();
  Real computePrincipalStress(Real r, Real theta);

  std::map<std::pair<unsigned int,unsigned int>,Real> _dist;
  std::map<std::pair<unsigned int,unsigned int>,Real> _value;

  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  const PostprocessorValue & _ki_value;
  const PostprocessorValue & _kii_value;
  const PostprocessorValue & _kiii_value;
  const MaterialProperty<SymmTensor> & _stress_tensor;
  Real _poissons_ratio;
  Real _m;
  Real _lambda;
  Real _yield_stress;
  Real _r_max;
  bool _has_symmetry_plane;
  Real _cutoff;
  Real _rho;

private:
  std::vector<Elem *> _intersected_elems;
  bool _treat_as_2d;
};

#endif
