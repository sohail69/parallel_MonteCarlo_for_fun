#pragma once
/**************************************\
! Queries the closest point on the
! boundary surface
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
! Notes on templated parameters:
! real : real number (potential any type)
! uint : Iteratable number (preferrably unsigned)
! sDIM : Spatial dimension of problem
! eDIM : Entity dimension of the boundary
\**************************************/
#include "../templatedGeometry/localVectorAlgebra.hpp"
#include "../templatedGeometry/geometricPrimitives.hpp"

template<typename real, typename uint, uint sDIM, uint eDIM>
real findClosestDistance(const boundary<real,sDIM,eDIM> & DirchBC, const Point<real,sDIM> & p0)
{
  Point<real,sDIM> sfp;
  findClosestSurfacePoint(DirchBC, p0, sfp);
  return pointEuclideanDistance(p0, sfp);
}


template<typename real, typename uint, uint sDIM, uint eDIM>
void findClosestSurfacePoint(const boundary<real,sDIM> & DirchBC
                           , const Point<real,sDIM>    & p0
                           , Point<real,sDIM>          & sfp)
{



}
