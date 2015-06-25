/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ElementsIntersectedByPlane.h"
#include "PlaneTracing.h"

#include <algorithm>

template<>
InputParameters validParams<ElementsIntersectedByPlane>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredParam<RealVectorValue>("point", "Point on plane from which to pick elements");
  params.addRequiredParam<RealVectorValue>("normal", "Normal vector of plane from which to pick elements");
  params.addParam<int>("plane_id", 1, "ID of the plane from which to pick elements");

  params.set<MultiMooseEnum>("execute_on") = "initial";

  return params;
}

ElementsIntersectedByPlane::ElementsIntersectedByPlane(const std::string & name, InputParameters parameters) : AuxKernel(name, parameters),

  _p0( getParam<RealVectorValue>("point") ),
  _normal( getParam<RealVectorValue>("normal") ),
  _plane_id( getParam<int>("plane_id"))

{
  if ( _nodal ) mooseError("From ElementsIntersectedByPlane: element on plane id must be an element property (use CONSTANT MONOMIAL)");
}

void
ElementsIntersectedByPlane::initialSetup()
{
  Moose::elementsIntersectedByPlane(_p0,_normal,_mesh.getMesh(),_intersected_elems);
}

void
ElementsIntersectedByPlane::compute()
{
  if (std::find(_intersected_elems.begin(),_intersected_elems.end(),_current_elem) != _intersected_elems.end())
  {
    /**
     * Update the variable data refernced by other kernels.
     * Note that this will update the values at the quadrature points too
     * (because this is an Elemental variable)
     */
    _var.setNodalValue(_plane_id);
  }
}

Real
ElementsIntersectedByPlane::computeValue()
{
  mooseError("From ElementsIntersectedByPlane: computeValue() is not defined");
  return 0;
}
