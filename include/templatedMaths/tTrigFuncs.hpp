#pragma once
#include "../globalMacros.hpp"
#include "tPowFuncs.hpp"

// Returns the Sine function
// of a number using a 15th order
// polynomial with an appropriate
// integer multiple of PI offset
template<typename Number>
FORCE_INLINE Number sin1(const Number theta, const Number a){
  const double PI=3.14159;
  constexpr double coeffs[8]={1.0                  //Term 1   1
                            ,-1.0/6.0              //Term 3   2
                            , 1.0/120              //Term 5   3
                            ,-1.0/5040.0           //Term 7   4
                            , 1.0/362880.0         //Term 9   5
                            ,-1.0/39916800.0       //Term 11  6
                            , 1.0/6227020800.0     //Term 13  7
                            ,-1.0/130764378000.0}; //Term 15  8
  Number x = (theta - a);
  Number x2 = x*x;
  Number sinX(0.0);
  #pragma unroll
  for(unsigned I=0; I<8; I++){
    sinX = sinX + coeffs[I]*x;
    x = x*x2;
  }
  return sinX;
};

// Returns the sin function
// of a variable
template<typename Number>
FORCE_INLINE Number sin(const Number theta){
  const double PI=3.14159265358979323846264338328;
  double a, dec;
  Number A(0.0);
  if(theta >=   PI)  dec =  2.0*PI;
  if(theta <= (-PI)) dec = -2.0*PI;
  for(int I=0; ;I++){
    a=double(I)*dec;
    if(( (theta-a) <= PI )and( (-PI)  <= (theta-a) )) break;
  }
  A=Number(a);
  return sin1<Number>(theta,A);
};

// Returns the Cosine function
// of a number
template<typename Number>
FORCE_INLINE Number cos(const Number theta){
  Number s = sin<Number>(theta);
  return sqrt<Number>(1.0 - s*s);
};

// Returns the tangent
// of a number
template<typename Number>
FORCE_INLINE Number tan(const Number theta){
  Number s = sin<Number>(theta);
  Number c = sqrt<Number>(1.0 - s*s);
  return s/c;
};

// The arc-cosine
// function
template<typename Number>
FORCE_INLINE Number acos(const Number x){
  Number one(1.00), half(0.5), Q(0.00), D(0.00);
  Number acosT(0.00), cs(0.00), sn(0.00);

  //Use iterative formula x - cos(theta) = 0
  //and solve with 10 iterations of
  //Halley's Method
  #pragma unroll
  for(unsigned I=0; I<10; I++){
    sn =  sin<Number>(acosT);
    cs =  cos<Number>(acosT);
    Q  = (x - cs)*sn;
    D  = (one - half*cs*cs - half*x*cs);
    acosT = acosT - Q/D;
  }
  return acosT;
};

// The arc-sine
// function
template<typename Number>
FORCE_INLINE Number asin(const Number x){
  Number one(1.00), half(0.5), Q(0.00), D(0.00);
  Number asinT(0.00), cs(0.00), sn(0.00);

  //Use iterative formula x - cos(theta) = 0
  //and solve with 10 iterations of
  //Halley's Method
  #pragma unroll
  for(unsigned I=0; I<10; I++){
    sn =  sin<Number>(asinT);
    cs =  cos<Number>(asinT);
    Q  = (sn - x)*cs;
    D  = (one - half*sn*sn - half*x*cs);
    asinT = asinT - Q/D;
  }
  return asinT;
};

// The arc-tangent
// function
template<typename Number>
FORCE_INLINE Number atan(const Number x){
  Number one(1.00), half(0.5), Q(0.00), D(0.00);
  Number atanT(0.00), cs(0.00), sn(0.00);
  return atanT;
};

// Not really implemented yet
//
template<typename Number>
FORCE_INLINE Number atan2(const Number theta){
  return Number(0.0);
};
