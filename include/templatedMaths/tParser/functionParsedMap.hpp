#pragma once
#include <map>
#include <vector>
#include <string>
#include <functional>
#include "../tCmath.hpp"

/*****************************************\
!
! Generates a map of strings to functions
! using the templated cmath library for
! single input functions
!
\*****************************************/
template<typename REAL, typename UINT>
void GenTCMathFuncMap(std::map<std::string, std::function<REAL(const REAL)> & SingleInputFuncs)
{
  //templated TrigFuncs
  SingleInputFuncs["cos"]   = cos<REAL>;   //Cosine function (Radians)
  SingleInputFuncs["sin"]   = sin<REAL>;   //Sine function (Radians)
  SingleInputFuncs["tan"]   = tan<REAL>;   //Tangent function (Radians)
  SingleInputFuncs["acos"]  = acos<REAL>;  //Arc-cosine function (Radians)
  SingleInputFuncs["asin"]  = asin<REAL>;  //Arc-sine function (Radians)
  SingleInputFuncs["atan"]  = atan<REAL>;  //Arc-tangent function (Radians)
  SingleInputFuncs["atan2"] = atan2<REAL>; //Arc-tangent 2nd quadrant function (Radians)

  //templated HypFuncs
  SingleInputFuncs["cosh"]  = cosh<REAL>;  //Hyperbolic cosine function
  SingleInputFuncs["sinh"]  = sinh<REAL>;  //Hyperbolic sine function
  SingleInputFuncs["tanh"]  = tanh<REAL>;  //Hyperbolic tangent function
  SingleInputFuncs["acosh"] = acosh<REAL>; //Hyperbolic arc-cosine function
  SingleInputFuncs["asinh"] = asinh<REAL>; //Hyperbolic arc-sine function
  SingleInputFuncs["atanh"] = atanh<REAL>; //Hyperbolic arc-tangent function

  //templated ExpLogFuncs
  SingleInputFuncs["exp"]   = exp<REAL>;   //Natural exponential
  SingleInputFuncs["log"]   = log<REAL>;   //Natural log function
  SingleInputFuncs["log10"] = log10<REAL>; //Log base 10 function
  SingleInputFuncs["exp2"]  = exp2<REAL>;  //Binary exponential function
  SingleInputFuncs["expm1"] = expm1<REAL>; //Natural exponential - 1
  SingleInputFuncs["log2"]  = log2<REAL>;  //Log base 2 function
  SingleInputFuncs["logb"]  = logb<REAL>;  //Log base 2 function
  SingleInputFuncs["log1p"] = log1p<REAL>; //Natural log function + 1

  //templated PowFuncs
  SingleInputFuncs["sqrt"]  = sqrt<REAL>;  //The sqaure-root function
  SingleInputFuncs["cbrt"]  = cbrt<REAL>;  //The cube-root function
};
