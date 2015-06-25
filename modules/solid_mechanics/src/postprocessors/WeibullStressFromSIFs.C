/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  This post processor calculates the Weibull stress
//
#include "WeibullStressFromSIFs.h"

template<>
InputParameters validParams<WeibullStressFromSIFs>()
{
  InputParameters params = validParams<GeneralPostprocessor>();

  params.addRequiredParam<UserObjectName>("crack_front_definition","The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index","The index of the point on the crack front corresponding to this q function");
  params.addRequiredParam<PostprocessorName>("KI_name","The name of the KI postprocessor");
  params.addRequiredParam<PostprocessorName>("KII_name","The name of the KII postprocessor");
  params.addRequiredParam<PostprocessorName>("KIII_name","The name of the KIII postprocessor");
  params.addRequiredParam<Real>("poissons_ratio","Poisson's ratio for the material.");
  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addRequiredParam<Real>("yield_stress","Yield stress of the material");
  params.addParam<unsigned int>("symmetry_plane", "Account for a symmetry plane passing through the plane of the crack, normal to the specified axis (0=x, 1=y, 2=z)");
  MooseEnum mesh_type("Random Regular","Random");
  params.addParam<MooseEnum>("weibull_stress_mesh_type",mesh_type,"The method used to generate the mesh points used to calculate the Weibull stress, options are: "+mesh_type.getRawNames());
  params.addParam<Real>("crack_tip_radius","Radius of the blunt crack tip.");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

WeibullStressFromSIFs::WeibullStressFromSIFs(const InputParameters & parameters):
    GeneralPostprocessor(parameters),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _has_crack_front_point_index(isParamValid("crack_front_point_index")),
    _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    _ki_value(getPostprocessorValue("KI_name")),
    _kii_value(getPostprocessorValue("KII_name")),
    _kiii_value(getPostprocessorValue("KIII_name")),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _m(getParam<Real>("m")),
    _lambda(getParam<Real>("lambda")),
    _yield_stress(getParam<Real>("yield_stress")),
    _has_symmetry_plane(isParamValid("symmetry_plane")),
    _mesh_type(getParam<MooseEnum>("weibull_stress_mesh_type"))
{
  if (isParamValid("crack_tip_radius"))
    _rho = 2 * getParam<Real>("crack_tip_radius");

  _mesh_area_radius = 50.0 * _rho;
  Moose::out<<"mesh area radius "<<_mesh_area_radius<<std::endl;
  _cutoff = _lambda * _yield_stress;
}

void
WeibullStressFromSIFs::initialize()
{
  generateRThetaMesh();
}

void
WeibullStressFromSIFs::execute()
{
  std::vector<Real> values;
  Real principal_stress = computeWeibullStress(_ki_value,_kii_value,_kiii_value);
  //  Moose::out<<"node "<<_crack_front_point_index<<" element "<<_current_elem->id()<<" "<<_qp<<" stress "<<principal_stress<<std::endl;
  //  Moose::out<<"node "<<_crack_front_point_index<<" element "<<_current_elem->id()<<" "<<_qp<<" r "<<_r<<" theta "<<_theta<<std::endl;
  Real value(0.0);
  if (principal_stress > _cutoff)
  {
    value = std::pow(principal_stress,_m);
  }
}

Real
WeibullStressFromSIFs::getValue()
{
  Real weibull_stress = computeWeibullStress(_ki_value,_kii_value,_kiii_value);
  
  if (_has_symmetry_plane)
    weibull_stress *= 2.0;

  return std::pow(weibull_stress,1.0/_m);
}

Real
WeibullStressFromSIFs::computeWeibullStress(const PostprocessorValue _ki, const PostprocessorValue _kii, const PostprocessorValue _kiii)
{
  Real weibull_stress;
  mooseAssert(_r.size() == _theta.size(), "mesh vector lengths are not identical");

  unsigned int i;
  unsigned int j;
  for (i=0, j=0; i<_r.size(), j<_theta.size(); ++i, ++j)
  {
    Real st2 = std::sin(_theta[j]/2);
    Real stt2 = std::sin(3*_theta[j]/2);
    Real ct2 = std::cos(_theta[j]/2);
    Real ctt2 = std::cos(3*_theta[j]/2);
    Real ki2PiR = _ki / std::sqrt(2*libMesh::pi*_r[i]);
    Real kii2PiR = _kii / std::sqrt(2*libMesh::pi*_r[i]);
    Real kiii2PiR = _kiii / std::sqrt(2*libMesh::pi*_r[i]);
    Real rho2R = _rho/(2*_r[i]);

    ColumnMajorMatrix stress(3,3);
    stress(0,0) = ki2PiR * ct2 * (1 - st2 * stt2) - ki2PiR * rho2R * ctt2 - kii2PiR * st2 * (2 + ct2 * ctt2) + kii2PiR * rho2R * stt2;
    stress(1,1) = ki2PiR * ct2 * (1 + st2 * stt2) + ki2PiR * rho2R * ctt2 + kii2PiR * st2 * ct2 * ctt2 - kii2PiR * rho2R * stt2;
    stress(0,1) = ki2PiR * st2 * ct2 * ctt2 - ki2PiR * rho2R * stt2 + kii2PiR * ct2 * (1 - st2 * stt2) - kii2PiR * rho2R * ct2;
    stress(1,0) = stress(0,1);
    stress(0,2) = -kiii2PiR * st2;
    stress(2,0) = stress(0,2);
    stress(1,2) = kiii2PiR * ct2;
    stress(2,1) = stress(1,2);
    stress(2,2) = _poissons_ratio * (2 * ki2PiR * ct2 - 2 * kii2PiR * st2);

//    Moose::out<<"ki2PiR "<<ki2PiR<<" kii2PiR "<<kii2PiR<<" kiii2PiR "<<kiii2PiR<<std::endl;
    ColumnMajorMatrix eval(3,1);
    ColumnMajorMatrix evec(3,3);
    stress.eigen(eval,evec);

    Real principal_stress = eval(2);
//    Moose::out<<"r "<<_r[i]<<" theta "<<_theta[j]<<std::endl;
//    stress.print();
//    Moose::out<<std::endl;
//    Moose::out<<"Stress "<<stress(0,0)<<std::endl;
    if (principal_stress > _cutoff)
    {
      weibull_stress += std::pow(principal_stress,_m);
      Moose::out<<"Principal stress "<<eval(2)<<" r "<<_r[i]<<" theta "<<_theta[j]<<" "<<_ki<<" "<<_kii<<" "<<_kiii<<std::endl;
      stress.print();
    }
  }
  
  return weibull_stress;
}

void
WeibullStressFromSIFs::generateRThetaMesh()
{
  if (_mesh_type == RANDOM)
    generateRandomMesh();

  else if (_mesh_type == REGULAR)
    generateRegularMesh();
}

void
WeibullStressFromSIFs::generateRegularMesh()
{
  int imax = 30;
  int jmax = 30;

  for (int i=0; i<imax; ++i)
  {
    for (int j=0; j<jmax; ++j)
    {
      Real r = _rho + sqrt((double) i/ (double) imax) *(_mesh_area_radius-_rho);
      _r.push_back(_rho + sqrt((double) i/ (double) imax) *(_mesh_area_radius-_rho));
      Real theta = 2*pi* (double) j / (double) jmax;
      _theta.push_back(2*pi* (double) j / (double) jmax);
//      Moose::out<<"r "<<r<<" theta "<<theta<<std::endl;
      //      double x = r * cos(theta);
//      double y = r * sin(theta);
//      std::cout<<x<<" "<<y<<std::endl;
//      std::cout<<r<<" "<<theta*180.0/pi<<" "<<i<<" "<<(double) j / (double) jmax<<std::endl;
//      double pstress = principalStress(_rho, r, theta, ki, kii, kiii);
    }
  }
}

void
WeibullStressFromSIFs::generateRandomMesh()
{
  srand(time(NULL));

  for (int i=0; i<=100000; ++i)
  {
    double u = ((double) rand() / ((double) RAND_MAX + 1.0));
    double v = ((double) rand() / ((double) RAND_MAX + 1.0));
    _r.push_back(_mesh_area_radius * sqrt( u ));
    _theta.push_back(2*pi*v - pi);
//    if (r <= _rho) continue; //FIX!
//    double x = r * cos(theta);
//    double y = r * sin(theta);
//    std::cout<<x<<" "<<y<<std::endl;
//    double pstress = principalStress(_rho, r, theta, ki, kii, kiii);
  }
}
