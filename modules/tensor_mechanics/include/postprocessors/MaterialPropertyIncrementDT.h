/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef MATERIALPROPERTYINCREMENTDT_H
#define MATERIALPROPERTYINCREMENTDT_H

#include "ElementPostprocessor.h"

//Forward Declarations
class MaterialPropertyIncrementDT;

template<>
InputParameters validParams<MaterialPropertyIncrementDT>();

/**
 * This postprocessor computes the maximum time step based on the
 * permissible change in a material property over a time step
 */
class MaterialPropertyIncrementDT: public ElementPostprocessor
{
public:
  MaterialPropertyIncrementDT(const InputParameters & parameters);

  virtual void initialize() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void execute() override;

protected:
  virtual void computeQpValue();

  const MaterialProperty<Real> & _scalar;
  const MaterialProperty<Real> & _scalar_old;
  const Real _max_increment;
  const unsigned int _growth_delay_steps;
  unsigned int _grow_dt_count;

  Real _value;
  unsigned int _qp;
};

#endif //MATERIALPROPERTYINCREMENTDT_H
