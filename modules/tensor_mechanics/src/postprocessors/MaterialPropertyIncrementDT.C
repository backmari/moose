/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MaterialPropertyIncrementDT.h"
#include "MooseTypes.h"

// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<MaterialPropertyIncrementDT>()
{
  InputParameters params = validParams<ElementPostprocessor>();
  params.addRequiredParam<MaterialPropertyName>("mat_prop", "The name of the scalar material property.");
  params.addRequiredParam<Real>("max_increment", "The maximum permissible increment in the material property");
  params.addParam<unsigned int>("growth_delay_steps", 5, "The number of steps to waith after cutting back the time step to allow the time step to grow");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

MaterialPropertyIncrementDT::MaterialPropertyIncrementDT(const InputParameters & parameters) :
    ElementPostprocessor(parameters),
    _scalar(getMaterialProperty<Real>(getParam<MaterialPropertyName>("mat_prop"))),
    _scalar_old(getMaterialPropertyOld<Real>(getParam<MaterialPropertyName>("mat_prop"))),
    _max_increment(parameters.get<Real>("max_increment")),
    _growth_delay_steps(parameters.get<unsigned int>("growth_delay_steps")),
    _grow_dt_count(0),
    _qp(0)
{
}

void
MaterialPropertyIncrementDT::initialize()
{
  _value = -std::numeric_limits<Real>::max();
}

void
MaterialPropertyIncrementDT::execute()
{
  for (_qp=0; _qp<_qrule->n_points(); _qp++)
    computeQpValue();
}

void
MaterialPropertyIncrementDT::computeQpValue()
{
  _value = std::max(_value, _scalar[_qp] - _scalar_old[_qp]);
}

Real
MaterialPropertyIncrementDT::getValue()
{
  //Only limit the max time step if the increment is too big
  const Real large = 1.0e3;
  Real dt_scale_factor = large;

  gatherMax(_value);
  if (_value > _max_increment)
  {
    dt_scale_factor = 0.5;
    _grow_dt_count = _growth_delay_steps;
  }
  else if (_grow_dt_count > 0)
  {
    _grow_dt_count -= 1;
    dt_scale_factor = 1.0;
  }

  return dt_scale_factor * _dt;
}

void
MaterialPropertyIncrementDT::threadJoin(const UserObject & y)
{
  const MaterialPropertyIncrementDT & pps = static_cast<const MaterialPropertyIncrementDT &>(y);

  _value = std::max(_value, pps._value);
}
