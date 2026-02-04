#pragma once
#include "../globalMacros.hpp"
#include "tPowFuncs.hpp"


// Returns the absolute value of
// a function, doesn't work for
// complex numbers
template<typename Number>
FORCE_INLINE Number abs(const Number z){
  return  sqrt<Number>(z*z);
};


// Returns the absolute value of
// a function, doesn't work for
// complex numbers
template<typename Number>
FORCE_INLINE Number fabs(const Number z){
  return  sqrt<Number>(z*z);
};


// Returns the multiply add function
// for single number types
template<typename Number>
FORCE_INLINE Number fma(const Number x, const Number y, const Number z){
  return x*y + z;
};

