#pragma once
#include <cstdio>
#include <vector>
#include <functional>
#include "../include/RNG.hpp"
#include "../include/localVectorAlgebra.hpp"
#include "../include/IntersectionDetection2D.hpp"

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
REAL lines( Vec2D<REAL> x ) {
   const REAL s = 8.0;
   return std::fmod( std::floor(s*x[0]), 2.0 );
}


template<typename REAL, typename UINT>
REAL solve(Vec2D<REAL> x0,                                  // evaluation point
           std::vector<Polyline2D<REAL>> boundaryDirichlet, // absorbing part of the boundary
           std::vector<Polyline2D<REAL>> boundaryNeumann,   // reflecting part of the boundary
           std::function<REAL(Vec2D<REAL>)> g,              // Dirichlet boundary values
           UINT seedVal)
{
   const REAL eps = 0.0001;     // stopping tolerance
   const REAL rMin = 0.0001;    // minimum step size
   const UINT nWalks = 65536;   // number of Monte Carlo samples
   const UINT maxSteps = 65536; // maximum walk length

   printf("Hello I am seed:  %d",seedVal);
   REAL sum = 0.0; // running sum of boundary contributions
   for( UINT i = 0; i < nWalks; i++ ) {
      Vec2D<REAL> x = x0;        // start walk at the evaluation point
      Vec2D<REAL> n={ 0.0, 0.0 };// assume x0 is an interior point, and has no normal
      bool onBoundary = false;   // flag whether x is on the interior or boundary

      REAL r, dDirichlet, dSilhouette; // radii used to define star shaped region
      UINT steps = 0;
      uint32_t rqd_seed = 0UL + uint32_t(seedVal);
      do { // loop until the walk hits the Dirichlet boundary
         // compute the radius of the largest star-shaped region
         dDirichlet = distancePolylines<REAL,UINT>( x, boundaryDirichlet );
         dSilhouette = silhouetteDistancePolylines<REAL,UINT>( x, boundaryNeumann );
         r = std::max( rMin, std::min( dDirichlet, dSilhouette ));

         // intersect a ray with the star-shaped region boundary
//         REAL theta = random<double>( -M_PI, M_PI ); //Doesn't work for me
         REAL theta = my_random<double>( -M_PI, M_PI, rqd_seed);
         rqd_seed = randqd_uint32(rqd_seed);

         if( onBoundary ) { // sample from a hemisphere around the normal
            theta = theta/2. + angleOf2DVec<REAL>(n);
         }
         Vec2D<REAL> v{ cos(theta), sin(theta) }; // unit ray direction
         x = intersectPolylines<REAL,UINT>( x, v, r, boundaryNeumann, n, onBoundary );

         steps++;
      }
      while(dDirichlet > eps && steps < maxSteps);
      //stop if we hit the Dirichlet boundary, or the walk is too long
      if( steps >= maxSteps ) printf("Hit max steps \n");
      sum += g(x); // accumulate contribution of the boundary value
   }
   return sum/nWalks; // Monte Carlo estimate
}
