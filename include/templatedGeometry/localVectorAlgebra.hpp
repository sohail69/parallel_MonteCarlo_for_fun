#pragma once
#include "geometricPrimitives.hpp"

template<typename real, size_t sdim>
real pointEuclideanDistance(const Point<real,sdim> & p0, const Point<real,sdim> & p1)
{
  real dist(0.00);
  for(size_t I=0; I<sdim; I++) dist += (p1[I]-p0[I])*(p1[I]-p0[I]);
  return std::sqrt(dist);
};


/*
//Facades used for making
//the code more clear
template<typename real, size_t sdim>
using Point = std::array<real,sdim>

template<typename real, size_t sdim>
using VecND = std::array<real,sdim>
*/
