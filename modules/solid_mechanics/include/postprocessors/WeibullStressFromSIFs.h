/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef WEIBULLSTRESSFROMSIFS_H
#define WEIBULLSTRESSFROMSIFS_H

#include "GeneralPostprocessor.h"
#include "CrackFrontDefinition.h"

//Forward Declarations
class WeibullStressFromSIFs;

template<>
InputParameters validParams<WeibullStressFromSIFs>();

/**
 * This postprocessor computes the Weibull stress
 *
 */
class WeibullStressFromSIFs:
  public GeneralPostprocessor
{
public:
  WeibullStressFromSIFs(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual Real getValue();

protected:
  enum MESH_TYPE
  {
    RANDOM,
    REGULAR
  };

  Real computeWeibullStress(PostprocessorValue ki, PostprocessorValue kii, PostprocessorValue kiii);
  void generateRegularMesh();
  void generateRandomMesh();
  void generateRThetaMesh();
  const CrackFrontDefinition * const _crack_front_definition;
  bool _has_crack_front_point_index;
  const unsigned int _crack_front_point_index;
  const PostprocessorValue & _ki_value;
  const PostprocessorValue & _kii_value;
  const PostprocessorValue & _kiii_value;
  Real _poissons_ratio;
  Real _m;
  Real _lambda;
  Real _yield_stress;
  bool _has_symmetry_plane;
  Real _cutoff;
  MooseEnum _mesh_type;
  Real _rho;
  Real _mesh_area_radius;
  std::vector<Real> _r;
  std::vector<Real> _theta;
};

#endif //WEIBULLSTRESSFROMSIFS_H
