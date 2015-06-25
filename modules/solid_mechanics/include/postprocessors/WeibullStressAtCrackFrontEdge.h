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

class WeibullStressAtCrackFrontEdge;
class SymmTensor;

template<>
InputParameters validParams<WeibullStressAtCrackFrontEdge>();

class WeibullStressAtCrackFrontEdge :
  public ElementIntegralPostprocessor,
  public MaterialTensorCalculator
{
public:
  WeibullStressAtCrackFrontEdge(const std::string & name, InputParameters parameters);

  virtual void initialize();
  virtual void execute();
//  virtual void threadJoin(const UserObject & u );
  virtual Real getValue();

protected:

  virtual Real computeQpIntegral();

  std::map<std::pair<unsigned int,unsigned int>,Real> _dist;
  std::map<std::pair<unsigned int,unsigned int>,Real> _value;

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

private:
  std::vector<Elem *> _intersected_elems;
};

#endif
