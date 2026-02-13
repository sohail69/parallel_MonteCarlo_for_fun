#pragma once
#include <random>
#include <functional>
#include "../globalMacros.hpp"


/**************************************\
!
! General interface for RNG's
! used in the model
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/

//Gets the random number then normalises
//it on the range of interest
template<typename real, typename RNGData>
FORCE_INLINE real RNG_reNormalise(std::function<REAL(const RNGData&)> rngNormalised
                                , const RNGData & seedData
                                , const real    & rMin
                                , const real    & rMax)
{
  return  rngNormalised(seedData)*(rMax-rMin) + rMin;
};


/**************************************\
!
! Original RNG used by the WoStr
! simplified example
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
//Data for Original RNG
struct OG_randomData
{
  const size_t rand_max = RAND_MAX;
  size_t seedval=0;
};

//Original RNG update
template<typename real>
FORCE_INLINE void OG_randomUpdate(OG_randomData & rqd_seed){
  const real rRandMax = 1.0/real(rqd_seed.rand_max);
  OG_randomData.seedval = rand();
}


//Original RNG renormalisation/re-Range
template<typename real>
FORCE_INLINE real OG_randomNormalised(const OG_randomData & rqd_seed){
  const real rRandMax = 1.0/real(rqd_seed.rand_max);
  return real(OG_randomData.seedval)*rRandMax;
}

/**************************************\
!
! Generate random number
! using a linear congruential
! unsigned int 32-bit generator
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
struct LCG32_rngData
{
  const size_t rand_max = (1L << 32);
  uint32_t seedval=0UL;
};

template<typename real>
FORCE_INLINE void LCG32_rngUpdate(LCG32_rngData & rqd_seed){
  rqd_seed.seedval  = (uint32_t) (1664525UL * rqd_seed.seedval + 1013904223UL);
};


// returns a random value in the range [rMin,rMax]
template<typename REAL>
FORCE_INLINE REAL LCG32_rngNormalised(const LCG32_rngData & rqd_seed){
  const real rRandMax = 1.0/real(rqd_seed.rand_max);
  return real(OG_randomData.seedval)*rRandMax;
};
