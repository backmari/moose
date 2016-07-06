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
#include "CrackFrontDefinition.h"

//Forward Declarations
class WeibullStress;
class PiecewiseLinear;

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
  WeibullStress(const InputParameters & parameters);
  //  virtual ~WeibullStress();
  virtual Real getValue();

protected:
  virtual void initialSetup();
  virtual Real computeQpIntegral();
//  virtual Real computeIntegral();
  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  //  MaterialTensorCalculator * _von_mises_calculator;
  //  MaterialTensorCalculator * _max_princ_calculator;
  const MaterialProperty<SymmTensor> & _stress_tensor;
  Real _m;
  Real _lambda;
  Real _yield_stress;
  Real _r_max;
  bool _has_symmetry_plane;
  Real _cutoff;
  //  Function * _cutoff_function;

private:
  bool _treat_as_2d;
};

#endif //WEIBULLSTRESS_H
