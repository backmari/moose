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

#ifndef REMOVERIGIDBODYMOTIONFUNCTIONDIRICHLETBC_H
#define REMOVERIGIDBODYMOTIONFUNCTIONDIRICHLETBC_H

#include "FunctionDirichletBC.h"

//Forward Declarations
class RemoveRigidBodyMotionFunctionDirichletBC;

template<>
InputParameters validParams<RemoveRigidBodyMotionFunctionDirichletBC>();

/**
 * Defines a boundary condition that forces the value to be a user specified
 * function at the boundary.
 */
class RemoveRigidBodyMotionFunctionDirichletBC : public FunctionDirichletBC
{
public:
  RemoveRigidBodyMotionFunctionDirichletBC(const InputParameters & parameters);

protected:
  /**
   * Evaluate the function at the current quadrature point and timestep.
   */
  virtual Real f();

  RealVectorValue _point;
};

#endif //REMOVERIGIDBODYMOTIONFUNCTIONDIRICHLETBC_H
