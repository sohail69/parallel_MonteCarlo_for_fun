#pragma once
#include <functional>
#include "../globalMacros.hpp"
#include "RNG.hpp"

/**************************************\
!
! Random unit vector that points from 
! the centre to a point on the surface
! of a unit sphere
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
template<typename real, size_t sdim>
FORCE_INLINE VecND<real,sdim> sampleUnitSphereUniform( VecND<real,sdim> randVecUnit)
{
  //Zero out the vector
  VecND<real,sdim> uniformSphereUnitRandVec;
  for(size_t i=0; i<sdim; i++) uniformSphereUnitRandVec[i] = real(0.00);

  //Use random numbers to generate
  //vectors on the unit N-Sphere
  //for dims 1:3 ()
  if(sdim == 1) uniformSphereUnitRandVec[0] = ((randVecUnit[0] < 0.5) ? -1.00: 1.00);
  if(sdim == 2){
    real phi = 2.0*M_PI*randVecUnit[0];
    uniformSphereUnitRandVec = {std::cos(phi), std::sin(phi)};
  }
  if(sdim == 3){
    real z = 1.0 - 2.0*randVecUnit[0];
    real r = std::sqrt(std::max(0.00, 1.00 - z*z));
    real phi = 2.0*M_PI*randVecUnit[1];
    uniformSphereUnitRandVec = {r*std::cos(phi), r*std::sin(phi), z};
  }
  return uniformSphereUnitRandVec;
};

/**************************************\
!
! Random unit vector that points to a 
! point on the surface of a concentric
! iunit disk
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
/*
template<typename real>
FORCE_INLINE VecND<real,sdim>  sampleUnitDiskConcentric(VecND<real,sdim> randVecUnit)
{
  //Zero out the vector
  VecND<real,sdim> uniformHemiSphereUnitRandVec;
  for(size_t i=0; i<sdim; i++) uniformHemiSphereUnitRandVec[i] = real(0.00);
  if(sdim == 2){
    real u1 = 2.00*randVecUnit[0] - 1.00;
    real z = std::sqrt(std::max(0.00, 1.00 - u1*u1));
    uniformHemiSphereUnitRandVec = {u1, z};
  }
//---------------------- Repair this, start here
  if(sdim == 3){
    real d2=0.00;
    VecND<real,2> d = sampleUnitDiskConcentric(randVecUnit);
    for(size_t i=0; i<sdim; i++) d2 += d[i];
    real z = std::sqrt(std::max(0.00, 1.00 - d2));
    uniformHemiSphereUnitRandVec = {d[0], d[1], z};
  }
//----------------------


  return uniformHemiSphereUnitRandVec;
};*/

/**************************************\
!
! Random unit vector that points to a 
! point on the surface of a unit hemi
! sphere
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
template<typename real, size_t sdim>
FORCE_INLINE VecND<real,sdim>  sampleUnitHemisphereCosine(VecND<real,sdim> randVecUnit)
{
  //Zero out the vector
  VecND<real,sdim> uniformHemiSphereUnitRandVec;
  for(size_t i=0; i<sdim; i++) uniformHemiSphereUnitRandVec[i] = real(0.00);
  if(sdim == 2){
    real u1 = 2.00*randVecUnit[0] - 1.00;
    real z = std::sqrt(std::max(0.00, 1.00 - u1*u1));
    uniformHemiSphereUnitRandVec = {u1, z};
  }

/***
  if(sdim == 3){
    real d2=0.00;
    VecND<real,2> d = sampleUnitDiskConcentric(randVecUnit);
    for(size_t i=0; i<sdim; i++) d2 += d[i];
    real z = std::sqrt(std::max(0.00, 1.00 - d2));
    uniformHemiSphereUnitRandVec = {d[0], d[1], z};
  }
***/


  return uniformHemiSphereUnitRandVec;
};


/**************************************\
!
! Random unit vector that points from 
! the centre to a point inside
! of a sphere
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
template<typename real, size_t sdim>
FORCE_INLINE VecND<real,sdim> sampleUnitBallUniform(VecND<real,sdim> randVecUnit)
{
  VecND<real,sdim> uniformBallUnitRandVec = sampleUnitSphereUniform<real,sdim>(randVecUnit);
  real r = std::pow( randVecUnit[sdim-1], 1.00/real(sdim) );
  for(size_t i=0; i<sdim; i++) uniformBallUnitRandVec[i] *= r;
  return uniformBallUnitRandVec;
};

/**************************************\
!
! The probability density function on
! a sample uniform ball
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
template<typename real, size_t sdim>
FORCE_INLINE real pdfSampleBallUniform(real r)
{
  real pdfUN = 0.00;
  real pdf2D = 1.0/(M_PI*r*r);
  real pdf3D = 3.0/(4.0f*M_PI*r*r*r);
  return ((sdim==2) ? pdf2D : ((sdim==3) ? pdf3D : pdfUN) );
};














