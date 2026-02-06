#pragma once
#include <random>
#include "../globalMacros.hpp"

// returns a random value in the range [rMin,rMax]
template<typename REAL>
FORCE_INLINE REAL random(REAL rMin, REAL rMax){
   const REAL rRandMax = 1.0/REAL(RAND_MAX);
   REAL u = rRandMax*REAL(rand());
   return u*(rMax-rMin) + rMin;
}


//
// Generate random number
// using a linear congruential
// unsigned int 32-bit generator
//
FORCE_INLINE uint32_t randqd_uint32(uint32_t rqd_seed){
    rqd_seed = (uint32_t) (1664525UL * rqd_seed + 1013904223UL);
    return rqd_seed;
};


 // returns a random value in the range [rMin,rMax]
template<typename REAL>
FORCE_INLINE REAL my_random(REAL rMin, REAL rMax, uint32_t & rqd_seed){
   const REAL rRandMax = 1.0/REAL(1L << 32);
   REAL u = rRandMax*REAL(randqd_uint32(rqd_seed));
   return u*(rMax-rMin) + rMin;
}

