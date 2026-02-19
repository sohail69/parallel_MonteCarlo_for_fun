#pragma once
#include <random>
#include <functional>
#include "../globalMacros.hpp"


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

  return VecND<real,sdim>();
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
