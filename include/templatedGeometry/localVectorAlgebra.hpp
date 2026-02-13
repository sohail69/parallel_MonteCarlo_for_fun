#pragma once
#include <cmath>
#include "geometricPrimitives.hpp"

// Euclidean distance between two
// points of arbitrary positive
// spatial dimension
template<typename real, size_t sdim>
real pointEuclideanDistance(const Point<real,sdim> & p0, const Point<real,sdim> & p1)
{
  real dist(0.00);
  for(size_t I=0; I<sdim; I++) dist += (p1[I]-p0[I])*(p1[I]-p0[I]);
  return std::sqrt(dist);
};

// Finds a unit n-spherical unit vector
// from a list of angles (right now only
// works for: 3-D, 2-D and 1-D)
template<typename real, size_t sdim>
VecND<real,sdim> getUnitVectorFromAngles(VecND<real,sdim-1> angles)
{
  if(sdim=1) return VecND<real,sdim>({1.00});
  if(sdim=2) return VecND<real,sdim>({ std::cos(angles[0]), std::sin(angles[0]) });
  if(sdim=3) return VecND<real,sdim>({ std::sin(angles[0])*std::cos(angles[1])
                                     , std::sin(angles[1])*std::cos(angles[0])
                                     , std::cos(angles[0]) });
  return VecND<real,sdim>(); //Default values
};
