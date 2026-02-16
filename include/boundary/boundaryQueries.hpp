#pragma once
/**************************************\
! Queries the the boundary for various
! points of interest
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
! Notes on templated parameters:
! real : real number (potential any type)
! sDIM : Spatial dimension of problem
! eDIM : Entity dimension of the boundary
\**************************************/
#include "../templatedGeometry/localVectorAlgebra.hpp"
#include "../templatedGeometry/geometricPrimitives.hpp"
#include "boundary/IntersectionDetection2D.hpp"


/**************************************\
!
! Dirchelet boundary intersection
! query
!
\**************************************/
// Find the nearest boundary point
// intersection (direct boundary
// intersection)
template<typename real, size_t sdim, size_t edim>
real findClosestDistance(const boundary<real,sDIM,eDIM> & DirchBC, const Point<real,sDIM> & p0)
{
  Point<real,sDIM> sfp;
  findClosestSurfacePoint(DirchBC, p0, sfp);
  return pointEuclideanDistance(p0, sfp);
}

// Find the boundary which
// is closest to the given
// point
template<typename real, size_t sdim, size_t edim>
void findClosestSurfacePoint(const boundary<real,sDIM,edim> & DirchBC
                           , const Point<real,sDIM>    & p0
                           , Point<real,sDIM>          & sfp)
{


}

/**************************************\
!
! Neumann boundary silhouette point
! intersection query
!
\**************************************/
// Find the nearest boundary 
// silhouette point to solve
// the problem
template<typename real, size_t sdim, size_t edim>
real findSilhouettePointDistance(const boundary<real,sDIM,eDIM> & NeumBC, const Point<real,sDIM> & p0)
{
  Point<real,sDIM> sfp;
  findClosestSilhouettePoint(DirchBC, p0, sfp);
  return pointEuclideanDistance(p0, sfp);
};

// Find the boundary which
// is closest to the given
// point
template<typename real, size_t sdim, size_t edim>
void findClosestSilhouettePoint(const boundary<real,sDIM,edim> & NeumBC
                              , const Point<real,sDIM>    & p0
                              , Point<real,sDIM>          & sfp)
{


};

/**************************************\
!
! Closest ray point 
!
\**************************************/
// Find the nearest boundary 
// silhouette point to solve
// the problem
template<typename real, size_t sdim, size_t edim>
real findClosestRayIntersectionPointDistance(const boundary<real,sDIM,eDIM> & NeumBC, const Point<real,sDIM> & p0)
{
  Point<real,sDIM> sfp;
  findClosestSilhouettePoint(DirchBC, p0, sfp);
  return pointEuclideanDistance(p0, sfp);
};

// Find the boundary which
// is closest to the given
// point
template<typename real, size_t sdim, size_t edim>
void findClosestRayIntersectionPoint(const boundary<real,sDIM,edim> & NeumBC
                                    , const Point<real,sDIM>    & p0
                                    , Point<real,sDIM>          & sfp)
{


};



/**************************************\
!
! Neumann boundary silhouette
! intersection query
!
\**************************************/
// Checks for the intersection
// of the the point on the
// closest line
template<typename REAL, typename UINT>
FORCE_INLINE Vec2D<REAL> intersectPolylines( Vec2D<REAL> x, Vec2D<REAL> v, REAL r,
                         const std::vector<Polyline2D<REAL>>& P,
                         Vec2D<REAL> & n, bool & onBoundary ){
   REAL tMin = r; // smallest hit time so far
   n = Vec2D<REAL>({ 0.0, 0.0 }); // first hit normal
   onBoundary = false; // will be true only if the first hit is on a segment
   for( UINT i = 0; i < P.size(); i++ ) { // iterate over polylines
      for( UINT j = 0; j < P[i].size()-1; j++ ) { // iterate over segments
         const REAL c = 1e-5; // ray offset (to avoid self-intersection)
         REAL t = rayIntersection<REAL>(x + c*v, v, P[i][j], P[i][j+1] );
         if( t < tMin ) { // closest hit so far
            tMin = t;
            rotate90<REAL>(P[i][j+1] - P[i][j], n); // get normal
            n /= length(n); // make normal unit length
            onBoundary = true;
         }
      }
   }
   return x + tMin*v; // first hit location
}

