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
template<typename real, typename RNGData, size_t sdim>
FORCE_INLINE VecND<real,sdim> sampleUnitSphereUniform(std::function<void(RNGData&)> rngUpdate // RN-update
                                                    , RNGData & seedData)                     // RN-data
{
  constexpr size_t hdim = ((sdim-1)>0) ? (sdim-1) : 1;
  VecND<real,hdim> randVec;
  VecND<real,sdim> uniformUnitRandVec;

  //Generate a series of random values
  //on the unit line (general case)
  for(size_t i=0; i<hdim; i++){
    rngUpdate(seedData);
    randVec[i] = RNG_reNormalise<real,RNGData>(seedData, 0.00, 1.00);
  };

  //Zero out the vector
  for(size_t i=0; i<sdim; i++) uniformUnitRandVec[i] = real(0.00);

  //Use random numbers to generate
  //vectors on the unit N-Sphere
  //for dims 1:3
  if(sdim == 1) uniformUnitRandVec[0] = ((randVec[0] < 0.5) ? -1.00: 1.00);
  if(sdim == 2){
    real phi = 2.0*M_PI*randVec[0];
    uniformUnitRandVec = {std::cos(phi), std::sin(phi)};
  }
  if(sdim == 3){
    real z = 1.0 - 2.0*randVec[0];
    real r = std::sqrt(std::max(0.00, 1.00 - z*z));
    real phi = 2.0*M_PI*randVec[1];
    uniformUnitRandVec = {r*std::cos(phi), r*std::sin(phi), z};
  }
  return uniformUnitRandVec;
};


/**************************************\
!
! Random unit vector that points to a 
! point inside a unit sphere
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
template<typename real, typename RNGData, size_t sdim>
FORCE_INLINE VecND<real,sdim>  sampleUnitHemisphereCosine(std::function<void(RNGData&)> rngUpdate // RN-update
                                                        , RNGData & seedData)                     // RN-data
{
  return VecND<real,sdim>();
};
