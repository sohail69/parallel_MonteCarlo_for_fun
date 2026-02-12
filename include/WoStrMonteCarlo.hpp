#pragma once
#include <random>
#include <cstdio>
#include <vector>
#include <functional>
#include "globalMacros.hpp"
#include "blockHyperMesh.hpp"
#include "Random/RNG.hpp"
#include "boundary/localVectorAlgebra.hpp"
#include "boundary/IntersectionDetection2D.hpp"



/**************************************\
! solves a Laplace equation Delta u = 0
! at x0, where the Dirichlet and Neumann
! boundaries are each given by a 
! collection of polylines, the Neumann
! boundary conditions are all zero, and
! the Dirichlet boundary conditions
! are given by a function g that can
! be evaluated at any point in space
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
template<typename REAL>
FORCE_INLINE REAL lines( Vec2D<REAL> x ) {
   const REAL s = 8.0;
   return std::fmod( std::floor(s*x[0]), 2.0 );
}


template<typename REAL, typename UINT, typename RNGData>
FORCE_INLINE REAL WoStr(Vec2D<REAL> x0,                                  // evaluation point
                        std::vector<Polyline2D<REAL>> boundaryDirichlet, // absorbing part of the boundary
                        std::vector<Polyline2D<REAL>> boundaryNeumann,   // reflecting part of the boundary
                        std::function<REAL(Vec2D<REAL>)> g,              // Dirichlet boundary values
                        RNGData seedData)
{
   const REAL eps = 0.0001;     // stopping tolerance
   const REAL rMin = 0.0001;    // minimum step size
   const UINT maxSteps = 65536; // maximum walk length

   uint32_t rqd_seed = 0UL + uint32_t(seedVal);
   Vec2D<REAL> x = x0;        // start walk at the evaluation point
   Vec2D<REAL> n={ 0.0, 0.0 };// assume x0 is an interior point, and has no normal
   bool onBoundary = false;   // flag whether x is on the interior or boundary

   REAL r, dDirichlet, dSilhouette; //radii used to define star shaped region
   UINT steps = 0;
   for(steps=0; (dDirichlet > eps) && (steps < maxSteps); steps++){
     // loop until the walk hits the Dirichlet boundary
     // compute the radius of the largest star-shaped region
     dDirichlet = distancePolylines<REAL,UINT>( x,boundaryDirichlet);
     dSilhouette = silhouetteDistancePolylines<REAL,UINT>( x, boundaryNeumann );
     r = std::max( rMin, std::min(dDirichlet,dSilhouette) );

     // intersect a ray with the star-shaped region boundary
     REAL theta = random<double>( -M_PI, M_PI );
//     rqd_seed = randqd_uint32(rqd_seed);
//     REAL theta = my_random<double>( -M_PI, M_PI, rqd_seed);

     // sample from a hemisphere around the normal
     if( onBoundary ) theta = theta/REAL(2.) + angleOf2DVec<REAL>(n);
     Vec2D<REAL> v{ cos(theta), sin(theta) }; // unit ray direction
     x = intersectPolylines<REAL,UINT>( x, v, r, boundaryNeumann, n, onBoundary );
     steps++;
   } //stop if we hit the Dirichlet boundary, or the walk is too long
   if( steps >= maxSteps ) printf("Hit max steps \n");
   return g(x); //Single random walk contribution of the boundary value
}
