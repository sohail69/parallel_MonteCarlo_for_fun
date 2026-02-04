#pragma once
#include "../globalMacros.hpp"
#include "tExpLogFuncs.hpp"

// Returns the error function
// function of a number
template<typename Number>
FORCE_INLINE Number erf(Number x){
  constexpr double coeff=(4.0/3.14159);
  constexpr double n1=279.0/10000000.0      , n2=-303923.0/10000000.0;
  constexpr double n3=34783.0/10000000.0    , n4=40793.0/10000000.0;
  constexpr double d1=-21941279.0/10000000.0, d2=3329407.0/10000000.0;
  Number Q, D, x1, one(1.00);
  x1 = x/(x+1);
  Q = n1*x1 + n2*x1*x1 + n3*x1*x1*x1 + n4*x1*x1*x1*x1;
  D = one + d1*x1 + d2*x1*x1;
  return sqrt<Number>(one - exp<Number>( -coeff*x*x*(one + Q/D) ) );
};


// Returns the complementary error
// function of a number
template<typename Number>
FORCE_INLINE Number erfc(Number x){
  Number one(1.00);
  return (one - erf(x));
};


// Returns the gamma
// function of a number
template<typename Number>
FORCE_INLINE Number tgamma(Number x){
  return Number(0.00);
};


// Returns the log gamma
// function of a number
template<typename Number>
FORCE_INLINE Number lgamma(Number x){
  return log<Number>(  tgamma<Number>( sqrt<Number>(x*x) ) );
};
