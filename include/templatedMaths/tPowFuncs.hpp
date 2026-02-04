#pragma once
#include "../globalMacros.hpp"
#include "tExpLogFuncs.hpp"

// Returns the pow(x,y) function
// using log and exponent relations
// x^y = exp(y.log(x))
template<typename Number>
FORCE_INLINE Number pow(const Number x, const Number y)
{
  Number YlogX  = y*log<Number>(x);
  Number XpowY = exp<Number>(YlogX);
  return XpowY;
};


// Returns the sqrt of a number
// by using 5 iterations of 
// halley's method
template<typename Number>
FORCE_INLINE Number sqrt(const Number x)
{
  double two=2.0, three=3.0;
  Number sqrtX(0.6*x), Q(0.0), D(0.0);

  #pragma unroll
  for(unsigned I=0; I<3; I++){
    Q = two*sqrtX*(sqrtX*sqrtX - x);
    D = three*sqrtX*sqrtX + x;
    sqrtX = sqrtX - Q/D;
  }
  return sqrtX;
};


// Returns the cbrt of a number
// by using 6 iterations of 
// halley's method
template<typename Number>
FORCE_INLINE Number cbrt(const Number x)
{
  double three=3.00, six=6.00;
  Number cbrtX(0.1*x), Q(0.0), D(0.0);

  #pragma unroll
  for(unsigned I=0; I<6; I++){
    Q = three*cbrtX*(cbrtX*cbrtX*cbrtX - x);
    D = six*cbrtX*cbrtX*cbrtX + three*x;
    cbrtX = cbrtX - Q/D;
  }
  return cbrtX;
};


// Gets the hypotenuse of
// a right-angled triangle
template<typename Number>
FORCE_INLINE Number hypot(const Number x, const Number y)
{
  return sqrt<Number>(x*x + y*y);
};
