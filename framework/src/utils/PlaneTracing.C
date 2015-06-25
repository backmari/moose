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

// Moose includes
#include "PlaneTracing.h"

// libMesh includes
#include "libmesh/plane.h"
#include "libmesh/point.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

#include <algorithm>

namespace Moose
{

void findElementsIntersectedByPlane(const Plane & plane, const MeshBase & mesh, std::vector<Elem *> & intersected_elems)
{
  // Loop over all elements to find elements intersected by the plane
  MeshBase::const_element_iterator el = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();
  for ( ; el != end_el ; ++el)
  {
    const Elem * elem = *el;
    std::vector<bool> above_plane;

    // Loop over all nodes belonging to this element and determine
    // whether they are on the same or different sides of the plane
    for (unsigned int i = 0; i < elem->n_nodes(); ++i)
    {
      const Node * node = elem->get_node(i);

      Point p;
      int dim = 3;
      for (unsigned int k=0; k<dim; ++k)
        p(k) = (*node)(k);

      bool node_above = plane.above_surface(p);
      above_plane.push_back(node_above);
    }

    if (std::find(above_plane.begin(),above_plane.end(),true)  != above_plane.end() &&
        std::find(above_plane.begin(),above_plane.end(),false) != above_plane.end())
      intersected_elems.push_back(const_cast<Elem *>(elem));
  }
}

void elementsIntersectedByPlane(const Point & p0, const Point & normal, const MeshBase & mesh, std::vector<Elem *> & intersected_elems)
{
  // Make sure our list is clear
  intersected_elems.clear();

  // Create plane from point and normal:
  Plane plane(p0,normal);

  // Find 'em!
  findElementsIntersectedByPlane(plane, mesh, intersected_elems);
}

void elementsIntersectedByPlane(const Point & p0, const Point & p1, const Point & p2, const MeshBase & mesh, std::vector<Elem *> & intersected_elems)
{
  // Make sure our list is clear
  intersected_elems.clear();

  // Create plane from three points:
  Plane plane(p0,p1,p2);

  // Find 'em!
  findElementsIntersectedByPlane(plane, mesh, intersected_elems);
}
}
