#pragma once
#include <array>
#include <vector>
#include "../globalMacros.hpp"
#include "localVectorAlgebra.hpp"
#include "../templatedMaths/tCmath.hpp"


// returns the closest point to x on a segment with endpoints a and b
template<typename REAL>
FORCE_INLINE Vec2D<REAL> closestPoint( Vec2D<REAL> x, Vec2D<REAL> a, Vec2D<REAL> b ) {
   Vec2D<REAL> u = b-a;
   REAL t = std::clamp( dot<REAL>(x-a,u)/dot<REAL>(u,u), 0.0, 1.0 );
   return (1.0-t)*a + t*b;
};

// returns true if the point b on the polyline abc is a silhoutte relative to x
template<typename REAL>
FORCE_INLINE bool isSilhouette(Vec2D<REAL> x, Vec2D<REAL> a, Vec2D<REAL> b, Vec2D<REAL> c ) {
   return ( cross(b-a,x-a) * cross(c-b,x-b) ) < 0;
};

// returns the time t at which the ray x+tv intersects segment ab,
// or infinity if there is no intersection
template<typename REAL>
FORCE_INLINE REAL rayIntersection( Vec2D<REAL> x, Vec2D<REAL> v, Vec2D<REAL> a, Vec2D<REAL> b ) {
   Vec2D<REAL> u = b - a;
   Vec2D<REAL> w = x - a;
   REAL d = cross(v,u);
   REAL s = cross(v,w) / d;
   REAL t = cross(u,w) / d;
   if (t > 0. && 0. <= s && s <= 1.) {
      return t;
   }
   return std::numeric_limits<REAL>::infinity();
}

// Returns distance from x to closest
// point on the given polylines P
template<typename REAL, typename UINT>
FORCE_INLINE REAL distancePolylines( Vec2D<REAL> x, const std::vector<Polyline2D<REAL>>& P ) {
   REAL d = std::numeric_limits<REAL>::infinity(); // minimum distance so far
   for( UINT i = 0; i < P.size(); i++ ) { // iterate over polylines
      for( UINT j = 0; j < P[i].size()-1; j++ ) { // iterate over segments
         Vec2D<REAL> y = closestPoint( x, P[i][j], P[i][j+1] ); // distance to segment
         d = std::min( d, length(x-y) ); // update minimum distance
      }
   }
   return d;
};

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

// returns distance from x to closest silhouette point on the given polylines P
template<typename REAL, typename UINT>
FORCE_INLINE REAL silhouetteDistancePolylines(Vec2D<REAL> x, const std::vector<Polyline2D<REAL>> &P ){
   REAL d = std::numeric_limits<REAL>::infinity(); // minimum distance so far
   for( UINT i = 0; i < P.size(); i++ ) { // iterate over polylines
      for( UINT j = 1; j < P[i].size()-1; j++ ) { // iterate over segment pairs
         if( isSilhouette( x, P[i][j-1], P[i][j], P[i][j+1] )) {
            d = std::min( d, length(x-P[i][j]) ); // update minimum distance
         }
      }
   }
   return d;
}

// these routines are not used by WoSt itself, but are rather used to check
// whether a given evaluation point is actually inside the domain
template<typename REAL, typename UINT>
FORCE_INLINE REAL signedAngle(const Vec2D<REAL> & x, const std::vector<Polyline2D<REAL>>& P )
{
  REAL Theta = REAL(0.00);
  for(UINT i = 0; i < P.size(); i++ ){
    for(UINT j = 0; j < P[i].size()-1; j++ ){
      Vec2D<REAL> y1 = P[i][j+1]-x;
      Vec2D<REAL> y2 = P[i][j]-x;
      REAL y2_norm = y2[0]*y2[0] + y2[1]*y2[1];
      Vec2D<REAL> y_normed = {(y1[0]*y2[0] + y1[1]*y2[1])/y2_norm, (y1[1]*y2[0] - y1[0]*y2[1])/y2_norm};
      Theta += angleOf2DVec<REAL>(y_normed);
    }
  }
  return Theta;
}

// Returns true if the point x is contained in the region bounded by the Dirichlet
// and Neumann curves.  We assume these curves form a collection of closed polygons,
// and are given in a consistent counter-clockwise winding order.
template<typename REAL, typename UINT>
FORCE_INLINE bool insideDomain(const Vec2D<REAL> & x
                , const std::vector<Polyline2D<REAL>>& boundaryDirichlet
                , const std::vector<Polyline2D<REAL>>& boundaryNeumann )
{
   REAL Theta = signedAngle<REAL,UINT>(x, boundaryDirichlet) + signedAngle<REAL,UINT>(x, boundaryNeumann);
   const REAL delta = 1e-4; // numerical tolerance
   return (std::abs(Theta-2.*M_PI) < delta); // boundary winds around x exactly once
}
