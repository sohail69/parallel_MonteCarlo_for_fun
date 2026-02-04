#pragma once
#include <array>
#include <vector>
#include "templatedMaths/tCmath.hpp"

// use array to implement 2D and 3D vectors
template<typename REAL>
using Vec2D = std::array<REAL,2>;

template<typename REAL>
using Vec3D = std::array<REAL,3>;

//Poly lines (used for 2D boundary representations)
template<typename REAL>
using Polyline2D = std::vector<Vec2D<REAL>>;

//Tri Soup (used for 3D triangulated boundary representations)
template<typename REAL>
using TriSoup3D = std::vector<std::array<Vec3D<REAL>,3>>;

//Point Soup (used for 3D point clouds)
template<typename REAL>
using PointSoup3D = std::vector<Vec3D<REAL>>;

// Operators acting on Vectors
template<typename REAL>
Vec2D<REAL> operator+(const Vec2D<REAL> & u, const Vec2D<REAL> & v){
  return Vec2D<REAL>({u[0]+v[0], u[1]+v[1]});
};

template<typename REAL>
Vec3D<REAL> operator+(const Vec3D<REAL> & u, const Vec3D<REAL> & v){
  return Vec3D<REAL>({u[0]+v[0], u[1]+v[1], u[2]+v[2]});
};

template<typename REAL>
Vec2D<REAL> operator-(const Vec2D<REAL> & u, const Vec2D<REAL> & v){
  return Vec2D<REAL>({u[0]-v[0], u[1]-v[1]});
};

template<typename REAL>
Vec3D<REAL> operator-(const Vec3D<REAL> & u, const Vec3D<REAL> & v){
  return Vec3D<REAL>({u[0]-v[0], u[1]-v[1], u[2]-v[2]});
};

//The length function
//for 2D and 3D vectors
template<typename REAL>
REAL length( Vec2D<REAL> u ) { return sqrt<REAL>( u[0]*u[0] + u[1]*u[1]); }

template<typename REAL>
REAL length( Vec3D<REAL> u ) { return sqrt<REAL>( u[0]*u[0] + u[1]*u[1] + u[2]*u[2]); }

//Calculate the angle of
//a 2D vector
template<typename REAL>
REAL angleOf(Vec2D<REAL> u) { return arg(u); }

//Rotate a 2D vector by
//90-degrees on the plane
template<typename REAL>
void rotate90(const Vec2D<REAL> & u , Vec2D<REAL> & v) { v[0]=-u[1]; v[1]=u[0]; };

// Calculate the inner product
// of 2D or 3D vectors
template<typename REAL>
REAL dot(Vec2D<REAL> u, Vec2D<REAL> v) { return u[0]*v[0] + u[1]*v[1]; };

template<typename REAL>
REAL dot(Vec3D<REAL> u, Vec3D<REAL> v) { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; };

//Find the magnitude of the cross
//product on a planar 2-D vector
template<typename REAL>
REAL cross(Vec2D<REAL> u, Vec2D<REAL> v) { return u[0]*v[1] - u[1]*v[0];}


//Find the polar angle of
// a 2D vector
template<typename REAL>
double angleOf2DVec(Vec2D<REAL> u) { return atan2(u[0], u[1]); };
