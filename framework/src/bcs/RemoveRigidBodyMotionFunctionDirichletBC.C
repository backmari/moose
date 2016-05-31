/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "RemoveRigidBodyMotionFunctionDirichletBC.h"
#include "Function.h"

template<>
InputParameters validParams<RemoveRigidBodyMotionFunctionDirichletBC>()
{
  InputParameters params = validParams<FunctionDirichletBC>();
  params.addRequiredParam<RealVectorValue>("point", "The point at which to compare the function.");
  return params;
}

RemoveRigidBodyMotionFunctionDirichletBC::RemoveRigidBodyMotionFunctionDirichletBC(const InputParameters & parameters) :
    FunctionDirichletBC(parameters),
    _point(getParam<RealVectorValue>("point"))
{
}

//TODO: Look up the value at the reference point one time instead of at every node

Real
RemoveRigidBodyMotionFunctionDirichletBC::f()
{
  return _func.value(_t, *_current_node) - _func.value(_t, _point);
}
