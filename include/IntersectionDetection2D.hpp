#pragma once
#include <vector>
#include <array>
#include "localVectorAlgebra.hpp"

// returns the closest point to x on a segment with endpoints a and b
template<typename REAL>
Vec2D<REAL> closestPoint( Vec2D<REAL> x, Vec2D<REAL> a, Vec2D<REAL> b ) {
   Vec2D<REAL> u = b-a;
   REAL t = clamp( dot(x-a,u)/dot(u,u), 0.0, 1.0 );
   return (1.0-t)*a + t*b;
};

// returns true if the point b on the polyline abc is a silhoutte relative to x
bool isSilhouette( Vec2D x, Vec2D a, Vec2D b, Vec2D c ) {
   return ( cross(b-a,x-a) * cross(c-b,x-b) ) < 0;
};

// returns the time t at which the ray x+tv intersects segment ab,
// or infinity if there is no intersection
template<typename REAL>
REAL rayIntersection( Vec2D<REAL> x, Vec2D<REAL> v, Vec2D<REAL> a, Vec2D<REAL> b ) {
   Vec2D<REAL> u = b - a;
   Vec2D<REAL> w = x - a;
   REAL d = cross(v,u);
   REAL s = cross(v,w) / d;
   REAL t = cross(u,w) / d;
   if (t > 0. && 0. <= s && s <= 1.) {
      return t;
   }
   return infinity;
}

// Returns distance from x to closest
// point on the given polylines P
template<typename REAL, typename UINT>
REAL distancePolylines( Vec2D<REAL> x, const vector<Polyline>& P ) {
   REAL d = infinity; // minimum distance so far
   for( UINT i = 0; i < P.size(); i++ ) { // iterate over polylines
      for( UINT j = 0; j < P[i].size()-1; j++ ) { // iterate over segments
         Vec2D<REAL> y = closestPoint( x, P[i][j], P[i][j+1] ); // distance to segment
         d = min( d, length(x-y) ); // update minimum distance
      }
   }
   return d;
};


// Checks for the intersection
// of the the point on the
// closest line
template<typename REAL, typename UINT>
Vec2D<REAL> intersectPolylines( Vec2D<REAL> x, Vec2D<REAL> v, REAL r,
                         const vector<Polyline>& P,
                         Vec2D<REAL> & n, bool & onBoundary ){
   REAL tMin = r; // smallest hit time so far
   n = Vec2D<REAL>({ 0.0, 0.0 }); // first hit normal
   onBoundary = false; // will be true only if the first hit is on a segment
   for( UINT i = 0; i < P.size(); i++ ) { // iterate over polylines
      for( UINT j = 0; j < P[i].size()-1; j++ ) { // iterate over segments
         const REAL c = 1e-5; // ray offset (to avoid self-intersection)
         REAL t = rayIntersection( x + c*v, v, P[i][j], P[i][j+1] );
         if( t < tMin ) { // closest hit so far
            tMin = t;
            n = rotate90( P[i][j+1] - P[i][j] ); // get normal
            n /= length(n); // make normal unit length
            onBoundary = true;
         }
      }
   }
   return x + tMin*v; // first hit location
}
