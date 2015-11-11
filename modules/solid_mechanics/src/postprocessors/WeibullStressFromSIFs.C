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
  params += validParams<RandomInterface>();

  params.addRequiredParam<UserObjectName>("crack_front_definition","The CrackFrontDefinition user object name");
  params.addParam<unsigned int>("crack_front_point_index","The index of the point on the crack front corresponding to this q function");
  params.addRequiredParam<PostprocessorName>("KI_name","The name of the KI postprocessor");
  params.addRequiredParam<PostprocessorName>("KII_name","The name of the KII postprocessor");
  params.addRequiredParam<PostprocessorName>("KIII_name","The name of the KIII postprocessor");
  MooseEnum tip_shape("Sharp Blunt","Sharp");
  params.addParam<MooseEnum>("crack_tip_shape", tip_shape, "Shape of the crack tip. Choices are: " + tip_shape.getRawNames());
  params.addRequiredParam<Real>("poissons_ratio","Poisson's ratio for the material.");
  params.addParam<Real>("m", "Weibull modulus");
  params.addParam<Real>("lambda", 2.0, "Stress cut-off scaling factor");
  params.addRequiredParam<Real>("yield_stress","Yield stress of the material");
  params.addParam<Real>("crack_tip_radius", 0.0, "Radius of the crack tip.");
  params.addParam<Real>("weibull_r_max","Max radius for Weibull stress calculation");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

WeibullStressFromSIFs::WeibullStressFromSIFs(const InputParameters & parameters):
    GeneralPostprocessor(parameters),
    RandomInterface(parameters, _fe_problem, _tid, true),
    _crack_front_definition(&getUserObject<CrackFrontDefinition>("crack_front_definition")),
    _has_crack_front_point_index(isParamValid("crack_front_point_index")),
    _crack_front_point_index(_has_crack_front_point_index ? getParam<unsigned int>("crack_front_point_index") : 0),
    _ki_value(getPostprocessorValue("KI_name")),
    _kii_value(getPostprocessorValue("KII_name")),
    _kiii_value(getPostprocessorValue("KIII_name")),
    _crack_tip_shape(getParam<MooseEnum>("crack_tip_shape")),
    _poissons_ratio(getParam<Real>("poissons_ratio")),
    _m(getParam<Real>("m")),
    _lambda(getParam<Real>("lambda")),
    _yield_stress(getParam<Real>("yield_stress")),
    _r_max(getParam<Real>("weibull_r_max"))
{
  if (isParamValid("crack_tip_radius"))
    _rho = 2 * getParam<Real>("crack_tip_radius");

  if (_rho < 0)
    mooseError("Crack tip radius cannot be negative");

  if (_rho > 0 && _crack_tip_shape == "Sharp")
    mooseError("A sharp crack tip cannot have a non-zero crack tip radius SIF");

  _cutoff = _lambda * _yield_stress;

  setRandomResetFrequency(EXEC_TIMESTEP_END);

  _npoints = 100000;

  generateRegularMesh();
}

void
WeibullStressFromSIFs::initialize()
{
}

void
WeibullStressFromSIFs::execute()
{
}

Real
WeibullStressFromSIFs::getValue()
{
  Real weibull_stress = computeWeibullStress(_ki_value,_kii_value,_kiii_value);

  return std::pow(weibull_stress, 1.0 / _m);
}

// Real
// WeibullStressFromSIFs::computeWeibullStress(const PostprocessorValue _ki, const PostprocessorValue _kii, const PostprocessorValue _kiii)
// {
//   unsigned int i = 0;
//   unsigned int n_fracture_zone = 0;
//   Real weibull_stress = 0;

//   Real x;
//   Real y;
//   Real r;
//   Real theta;
  
//   while (i < _npoints)
//   {
//     x = -_r_max + 2 * _r_max * getRandomReal();
//     y = -_r_max + 2 * _r_max * getRandomReal();
//     r = std::sqrt(x*x + y*y);
//     theta = std::atan2(y,x);
//     if (r < _rho/2 || r > _r_max)
//       continue;

//     // Real r = _r_max * std::sqrt(getRandomReal());
//     // if (r < _rho/2)
//     //   continue;
//     // Real theta = 2 * libMesh::pi * getRandomReal();

//     i++;
//     Real max_principal_stress = computePrincipalStress(r,theta);
// //    if (max_principal_stress > _cutoff)
//     if (max_principal_stress > _cutoff && max_principal_stress < 3*_yield_stress)
//     {
//       weibull_stress = weibull_stress + std::pow(max_principal_stress, _m);
//       n_fracture_zone++;
//     }
//   }

// //  Real edge_length = _crack_front_definition->getCrackFrontForwardSegmentLength(_crack_front_point_index);
// //  return edge_length * libMesh::pi * std::pow(_r_max, 2) * (weibull_stress / _npoints);
// //  Real volume = (n_fracture_zone / _npoints) * libMesh::pi * std::pow(_r_max, 2);
//   Real volume = libMesh::pi * (std::pow(_r_max, 2) - std::pow(0.5*_rho, 2));
//   Moose::out<<"n_fracture_zone "<<n_fracture_zone<<" _npoints "<<_npoints<<" volume "<<volume<<std::endl;
//   return (volume / _npoints) * weibull_stress;
// //  return weibull_stress / _npoints;
// }

Real
WeibullStressFromSIFs::computePrincipalStress(Real r, Real theta)
{
  Real st2 = std::sin(theta / 2);
  Real stt2 = std::sin(3 * theta / 2);
  Real ct2 = std::cos(theta / 2);
  Real ctt2 = std::cos(3 * theta / 2);
  Real ki2PiR =   _ki_value   / std::sqrt(2 * libMesh::pi * r);
  Real kii2PiR =  _kii_value  / std::sqrt(2 * libMesh::pi * r);
  Real kiii2PiR = _kiii_value / std::sqrt(2 * libMesh::pi * r);
  // Real ki2PiR =   7.185655e+01   / std::sqrt(2 * libMesh::pi * r);
//  Real ki2PiR =   0.0 / std::sqrt(2 * libMesh::pi * r);
//  Real kii2PiR =  0.0 / std::sqrt(2 * libMesh::pi * r);
//  Real kiii2PiR = 0.0 / std::sqrt(2 * libMesh::pi * r);
  Real rho2R = _rho / (2 * r);

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

  ColumnMajorMatrix eval(3,1);
  ColumnMajorMatrix evec(3,3);
  stress.eigen(eval,evec);

  return eval(2);
}

Real
WeibullStressFromSIFs::computeWeibullStress(const PostprocessorValue _ki, const PostprocessorValue _kii, const PostprocessorValue _kiii)
{
  Real weibull_stress = 0;

  unsigned int i;
  unsigned int j;
  unsigned int n = 0;

  for (i=0; i<_r.size(); ++i)
  {
    for(j=0; j<_theta.size(); ++j)
    {
      Real x = _r[i] * cos(_theta[j]);
      Real y = _r[i] * sin(_theta[j]);
      if (x < 0 && (y < _rho/2 && y > -_rho/2))
          continue;
//      else
//        Moose::out<<"theta "<<_theta[j]<<std::endl;
//        Moose::out<<"x y "<<x<<" "<<y<<std::endl;
      
      Real principal_stress = computePrincipalStress(_r[i],_theta[j]);

      if (principal_stress > _cutoff && principal_stress < 3*_yield_stress)
      {
        Real r_inner = _r[i] - 0.5 * _dr;
        Real area(0.0);
        if (i == 0 && _crack_tip_shape == "Sharp")
        {
          Real r_outer = _r[i] + 0.5 * _dr;
          area = (_dtheta/2)*std::pow(r_outer,2);
        }
        else
          area = _dr * r_inner * _dtheta;

        weibull_stress += area * std::pow(principal_stress,_m);
        n++;
      }
    }
  }
  return weibull_stress;
}

void
WeibullStressFromSIFs::generateRegularMesh()
{
  int imax = 1000;
  int jmax = 1000;
  _dr = (_r_max-0.5*_rho) / (Real) imax;
  _dtheta = 2*libMesh::pi / (Real) jmax;

  for (int i=0; i<imax; ++i)
  {
      Real r = 0.5*_rho + 0.5 * _dr + (Real) i * _dr;
      _r.push_back(r);
  }
  
  for (int j=0; j<jmax; ++j)
  {
    Real theta = (Real) j * _dtheta;
    _theta.push_back(theta);
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
    _r.push_back(_r_max * sqrt( u ));
    _theta.push_back(2*pi*v - pi);
//    if (r <= _rho) continue; //FIX!
//    double x = r * cos(theta);
//    double y = r * sin(theta);
//    std::cout<<x<<" "<<y<<std::endl;
//    double pstress = principalStress(_rho, r, theta, ki, kii, kiii);
  }
}
