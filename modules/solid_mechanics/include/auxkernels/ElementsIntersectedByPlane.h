/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELEMENTSINTERSECTEDBYPLANE_H
#define ELEMENTSINTERSECTEDBYPLANE_H

#include "AuxKernel.h"

class ElementsIntersectedByPlane : public AuxKernel
{

public:

  ElementsIntersectedByPlane(const std::string & name, InputParameters parameters);

  virtual ~ElementsIntersectedByPlane() {}

protected:

  virtual void initialSetup();
  virtual void compute();
  virtual Real computeValue();

private:

  const Point _p0;
  const Point _normal;
  int _plane_id;
  std::vector<Elem *> _intersected_elems;

};

template<>
InputParameters validParams<ElementsIntersectedByPlane>();

#endif // ELEMENTSINTERSECTEDBYPLANE_H
