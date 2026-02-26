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
! Date  : 21/02/2025
!
\**************************************/

//Gets the random number then normalises
//it on the range of interest
template<typename real, typename RNGData>
FORCE_INLINE real RNG_reNormalise(const RNGData & seedData
                                , const real    & rMin
                                , const real    & rMax)
{
  return  ( real( seedData.randomNumber*(rMax-rMin) )/real(seedData.rand_max) )  + rMin;
};

/**************************************\
!
! Original RNG used by the WoStr
! simplified example (just uses 
! std::random)
!
! Author: Sohail Rathore
! Date  : 21/02/2025
!
\**************************************/
//Data for Original RNG
struct OG_randomData
{
  const size_t rand_max = RAND_MAX;
  size_t randomNumber=0;
};

//Original RNG update
template<typename real>
FORCE_INLINE void OG_randomUpdate(OG_randomData & rqd_seed){
  const real rRandMax = 1.0/real(rqd_seed.rand_max);
  rqd_seed.randomNumber = rand();
}

/**************************************\
!
! Generate random number
! using a linear congruential
! unsigned int 32-bit generator
!
! Author: Sohail Rathore
! Date  : 21/02/2025
!
\**************************************/
struct LCG32_rngData
{
  const size_t rand_max = std::numeric_limits<uint32_t>::max();
  uint32_t randomNumber=0UL;
};

FORCE_INLINE void LCG32_rngUpdate(LCG32_rngData & rqd_seed){
  rqd_seed.randomNumber  = (uint32_t) (1664525UL * rqd_seed.randomNumber + 1013904223UL);
};

/**************************************\
!
! Generate random number
! using  xor-shift256++
!
! Author: Sohail Rathore
! Date  : 21/02/2025
!
\**************************************/
struct XORSHIFT256_rngData
{
  const size_t rand_max = std::numeric_limits<uint64_t>::max();
  uint64_t s[4];
  uint64_t randomNumber=1UL;
};

FORCE_INLINE uint64_t XORSHIFT256_rol64(const uint64_t & x, const int & k){
  return (x << k) | (x >> (64 - k));
};

FORCE_INLINE void XORSHIFT256_rngUpdate(XORSHIFT256_rngData & rqd_seed){
  rqd_seed.randomNumber = XORSHIFT256_rol64(rqd_seed.s[0] + rqd_seed.s[3],23) + rqd_seed.s[0];
  uint64_t t = rqd_seed.s[1] << 17;
  rqd_seed.s[2] ^= rqd_seed.s[0];
  rqd_seed.s[3] ^= rqd_seed.s[1];
  rqd_seed.s[1] ^= rqd_seed.s[2];
  rqd_seed.s[0] ^= rqd_seed.s[3];

  rqd_seed.s[2] ^= t;
  rqd_seed.s[3] ^= XORSHIFT256_rol64(rqd_seed.s[3], 45);
};
