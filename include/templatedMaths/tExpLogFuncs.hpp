#pragma once
#include "../globalMacros.hpp"

// Returns the exponential
// of a number
template<typename Number>
FORCE_INLINE Number exp(const Number x)
{
  Number bigNum(33554432.0), one(1.0), xCopy(x);
  Number expX = one + (xCopy/bigNum);

  #pragma unroll
  for(unsigned I=0; I<25; I++) expX = expX*expX; // m=2^24
  return expX;
};


// Returns the natural log of a 
// number using 10 iterations
// of Halley's method
template<typename Number>
FORCE_INLINE Number log(const Number x)
{
  //B = log(x) => x = exp(B)
  Number B(0.00), expB(0.00), const2(2.0);
  Number xPb(0.00), xMb(0.00);

  //Use iterative formula x-exp(b) = 0
  //and solve with 5 iterations of
  //Halley's Method
  #pragma unroll
  for(unsigned I=0; I<10; I++){
    expB =  exp<Number>(B);
    xPb  = (expB + x);
    xMb  = (expB - x);
    B = B - const2*(xMb/xPb);
  }
  return B;
};


// Returns the common log
// of a number (base 10)
template<typename Number>
FORCE_INLINE Number log10(const Number x)
{
  Number ln10 = log<Number>(Number(10.0));
  Number lnX = log<Number>(x);
  return lnX/ln10;
};


// Returns the exponential function
// base 2
template<typename Number>
FORCE_INLINE Number exp2(const Number x)
{
  Number ln2 = log<Number>(Number(2.0));
  return exp<Number>(x*ln2);
};


// Returns the exponential function
// minus 1
template<typename Number>
FORCE_INLINE Number expm1(const Number x)
{
  Number one(1.0);
  return exp<Number>(x) - one;
};


// Returns the binary log
// of a number (base 2)
template<typename Number>
FORCE_INLINE Number log2(const Number x)
{
  Number ln2 = log<Number>(Number(2.0));
  Number lnX = log<Number>(x);
  return lnX/ln2;
};


// Returns the binary log
// of a number (base 2)
template<typename Number>
FORCE_INLINE Number logb(const Number x)
{
  return  log2<Number>(x);
};


// Returns the value from significand
// and exponent
template<typename Number>
FORCE_INLINE Number ldexp(const Number x, const Number y)
{
  return  x*exp2<Number>(y);
};


// Returns the value of the
// natural log plus one
template<typename Number>
FORCE_INLINE Number log1p(const Number x)
{
  Number one(1.0);
  return  one + log<Number>(x);
};


// Returns the Scale significand using
// floating-point base exponent
template<typename Number>
FORCE_INLINE Number scalbn(const Number x, int n)
{
  double nd(n);
  Number N(nd);
  return ldexp(x, N);
};


// Returns the Scale significand using
// floating-point base exponent (long)
template<typename Number>
FORCE_INLINE Number scalbln(const Number x, long int n)
{
  double nd(n);
  Number N(nd);
  return ldexp(x, N);
};
