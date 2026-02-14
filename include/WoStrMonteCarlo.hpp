#pragma once
#include <random>
#include <cstdio>
#include <vector>
#include <functional>
#include "globalMacros.hpp"
#include "blockHyperMesh.hpp"
#include "Random/RNG.hpp"
#include "boundary/localVectorAlgebra.hpp"
#include "boundary/boundaryQueries.hpp"


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
template<typename real, typename uint, typename RNGData, size_t nVars, size_t sdim, size_t edim>
FORCE_INLINE VecND<real,nVars> WoStr(const Point<real,sdim>    & x0,            // evaluation point
                                     const boundary<real,sdim,edim> & Dirichlet,// absorbing part of the boundary
                                     const boundary<real,sdim,edim> & Neumann,  // reflecting part of the boundary
                                     std::function<VecND<real,nVars>(Point<real,sdim>)> g, // Greens function
                                     std::function<void(RNGData&)> rngUpdate,              // RN-update
                                     std::function<real(const RNGData&)> rngNormalised,    // RN-normalise
                                     RNGData & seedData)                                   // RN-data
{
   const real eps = 0.0001;     // stopping tolerance
   const real rMin = 0.0001;    // minimum step size
   const uint maxSteps = 65536; // maximum walk length

   Point<real,sdim> x = x0;      // start walk at the evaluation point
   VecND<real,sdim> n=real(0.00);// assume x0 is an interior point, and has no normal
   bool onBoundary = false;      // flag whether x is on the interior or boundary

   real r, dDirichlet, dSilhouette; //radii used to define star shaped region
   uint steps = 0;
   for(steps=0; (dDirichlet > eps) && (steps < maxSteps); steps++){
     // loop until the walk hits the Dirichlet boundary
     // compute the radius of the largest star-shaped region
     dDirichlet = distancePolylines<real,uint>( x,boundaryDirichlet);
     dSilhouette = silhouetteDistancePolylines<real,uint>( x, boundaryNeumann );
     r = std::max( rMin, std::min(dDirichlet,dSilhouette) );

     // intersect a ray with the star-shaped region boundary
     rngUpdate(seedData);
     real theta = RNG_reNormalise<real,RNGData>(rngNormalised, seedData, -M_PI, M_PI);

     if( onBoundary ) theta = theta/real(2.) + angleOf2DVec<real>(n); // sample from hemisphere around the normal

     VecND<real,sdim-1> angles
//     Vec2D<real> v{ cos(theta), sin(theta) };

     VecND<real,sdim> v = getUnitVectorFromAngles(VecND<real,sdim-1> angles);
     x = intersectPolylines<real,uint>( x, v, r, boundaryNeumann, n, onBoundary );
     steps++;
   } //stop if we hit the Dirichlet boundary, or the walk is too long
   if( steps >= maxSteps ) printf("Hit max steps \n");
   return g(x); //Single random walk contribution of the boundary value
}
