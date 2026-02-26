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
  if constexpr(sdim == 1) uniformSphereUnitRandVec[0] = ((randVecUnit[0] < 0.5) ? -1.00: 1.00);
  if constexpr(sdim == 2){
    real phi = 2.0*M_PI*randVecUnit[0];
    uniformSphereUnitRandVec = {std::cos(phi), std::sin(phi)};
  }
  if constexpr(sdim == 3){
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
template<typename real>
FORCE_INLINE VecND<real,2> sampleConcentricUnitDisk(VecND<real,3> randVecUnit)
{
  //Zero out the vector
  VecND<real,2> uniformConcentricDiskUnitRandVec = {0.00, 0.00};

  // map uniform random numbers to [-1,1]^2
  real u1 = 2.00*randVecUnit[0] - 1.00;
  real u2 = 2.00*randVecUnit[1] - 1.00;
  real theta, r;

  // handle degeneracy at the origin
  if (u1 != 0 && u2 != 0) {
    // apply concentric mapping to point
    if(std::abs(u1) > std::abs(u2)){
      r = u1;
      theta = 0.25*M_PI*(u2/u1);
    }else{
      r = u2;
      theta = 0.50*M_PI*(1.00 - 0.50*(u1/u2));
    }
    uniformConcentricDiskUnitRandVec = {r*std::cos(theta), r*std::sin(theta)};
  }
  return uniformConcentricDiskUnitRandVec;
};

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
  if constexpr(sdim == 2){
    real u1 = 2.00*randVecUnit[0] - 1.00;
    real z = std::sqrt(std::max(0.00, 1.00 - u1*u1));
    uniformHemiSphereUnitRandVec = {u1, z};
  }
  if constexpr(sdim == 3){
    real d2=0.00;
    VecND<real,2> d = sampleConcentricUnitDisk(randVecUnit);
    for(size_t i=0; i<sdim; i++) d2 += d[i];
    real z = std::sqrt(std::max(0.00, 1.00 - d2));
    uniformHemiSphereUnitRandVec = {d[0], d[1], z};
  }
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
! The probability density function for
! a sample on the uniform sphere
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
template<typename real, size_t sdim>
FORCE_INLINE real pdfSampleSphereUniform(real r)
{
  real pdfUN = 0.00;
  real pdf2D = 1.00/(2.00*M_PI*r);
  real pdf3D = 1.00/(4.00*M_PI*r*r);
  return ( (sdim==2) ? pdf2D : ((sdim==3) ? pdf3D : pdfUN) );
};

/**************************************\
!
! The probability density function for
! a sample on the uniform ball
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
template<typename real, size_t sdim>
FORCE_INLINE real pdfSampleUnitHemisphereCosine(real cosTheta)
{
  real pdfUN = 0.00;
  real pdf2D = cosTheta/2.00;
  real pdf3D = cosTheta/M_PI;
  return ( (sdim==2) ? pdf2D : ((sdim==3) ? pdf3D : pdfUN) );
};

/**************************************\
!
! The probability density function for
! a sample in the uniform unit ball
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
  return ( (sdim==2) ? pdf2D : ((sdim==3) ? pdf3D : pdfUN) );
};

/**************************************\
!
! Perform a coordinate transform on
! a given vector
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
template<typename real, size_t sdim>
FORCE_INLINE void transformCoordinates(const VecND<real,sdim> & n, VecND<real,sdim> & d)
{
  VecND<real,sdim> dtmp = d;
  if constexpr(sdim == 2){
    // compute orthonormal basis
    VecND<real,sdim> s = {n[1], -n[0]};
    //Perform Transform
    d[0] = dtmp[0]*s[0] + dtmp[1]*n[0];
    d[1] = dtmp[1]*s[1] + dtmp[1]*n[1];
  }
  if constexpr(sdim == 3){
    // compute orthonormal basis; source: https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    real sign = std::copysignf(1.00, n[2]);
    const real a = -1.0f/(sign + n[2]);
    const real b = n[0]*n[1]*a;
    VecND<real,sdim> b1 = {1.00 + sign*n[0]*n[0]*a, sign*b, -sign*n[0]};
    VecND<real,sdim> b2 = {b, sign + n[1]*n[1]*a, -n[1]};

    //Perform Transform
    d[0] = dtmp[0]*b1[0] + dtmp[1]*b2[0] + dtmp[2]*n[0];
    d[1] = dtmp[0]*b1[1] + dtmp[1]*b2[1] + dtmp[2]*n[1];
    d[2] = dtmp[0]*b1[2] + dtmp[1]*b2[2] + dtmp[2]*n[2];
  }
};








