#pragma once
#include <random>
#include <cstdio>
#include <vector>
#include <functional>
#include "globalMacros.hpp"
#include "blockHyperMesh.hpp"
#include "Random/RNG.hpp"
#include "Random/RandomSampling.hpp"
#include "boundary/localVectorAlgebra.hpp"
#include "boundary/boundaryQueries.hpp"


/**************************************\
! The walk on stars algorithm
! for a potentially arbitrary
! Partial differential equation
\**************************************/
template<typename real, typename RNGData, size_t sdim, size_t edim>
FORCE_INLINE Point<real,sdim> WoStr_walk(const Point<real,sdim>                  & x0         // evaluation point
                                       , std::vector<boundary<double,sdim,edim>> & boundaries // The boundaries
//                                     , std::unique_ptr<GreensFnBall<sdim>>     & greensFn   // Greens function
                                       , std::function<void(RNGData&)> rngUpdate              // RN-update
                                       , RNGData & seedData)                                  // RN-data
{
  //Star radius
  float starRadius;
//  ComputeStarRadius<real>(starRadius);

  //update the ball center and radius
//  greensFn->updateBall(state.currentPt, starRadius);

  // sample a direction uniformly
//  VecND<real,sdim> direction = sampleUnitSphereUniform<real,RNGData,sdim>(rngUpdate, seedData);



}


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
! dDirichlet = distancePolylines<real,unsigned>( x,boundaryDirichlet);
! dSilhouette = silhouetteDistancePolylines<real,unsigned>( x, boundaryNeumann );
! std::function<VecND<real,nVars>(Point<real,sdim>)> g, // Greens function
\**************************************/
template<typename real, typename RNGData, size_t sdim, size_t edim>
FORCE_INLINE Point<real,sdim> WoStr_point(const Point<real,sdim>    & x0           // evaluation point
                                        , const boundary<real,sdim,edim> & DirchBC // absorbing boundary
                                        , const boundary<real,sdim,edim> & NeumBC  // reflecting boundary
                                        , std::function<void(RNGData&)> rngUpdate  // RN-update
                                        , RNGData & seedData)                      // RN-data
{
   const real eps = 0.0001;     // stopping tolerance
   const real rMin = 0.0001;    // minimum step size
   const unsigned maxSteps = 65536; // maximum walk length

   Point<real,sdim> x = x0;      // start walk at the evaluation point
   VecND<real,sdim> n=real(0.00);// assume x0 is an interior point, and has no normal
   bool onBoundary = false;      // flag whether x is on the interior or boundary

   real r, dDirichlet, dSilhouette; //radii used to define star shaped region
   unsigned steps = 0;
   do{
     // loop until the walk hits the Dirichlet boundary
     // compute the radius of the largest star-shaped region
     dDirichlet  = findClosestDistance<real,sdim,edim>(DirchBC,x);
     dSilhouette = findSilhouettePointDistance<real,sdim,edim>(NeumBC,x);
/*
     r = std::max( rMin, std::min(dDirichlet,dSilhouette) );

     // intersect a ray with the star-shaped region boundary
     VecND<real,sdim-1> angles;
     rngUpdate(seedData); 
     angles[0] = RNG_reNormalise<real,RNGData>(seedData, -M_PI, M_PI);
     if( onBoundary ) theta = theta/real(2.) + angleOf2DVec<real>(n); // sample from hemisphere around the normal
     VecND<real,sdim> v = getUnitVectorFromAnglesVecND<real,sdim>(angles);
     x = intersectPolylines<real,unsigned>( x, v, r, boundaryNeumann, n, onBoundary );*/
     steps++;
   }while((dDirichlet > eps) && (steps < maxSteps));
   //stop if we hit the Dirichlet boundary, or the walk is too long
   if( steps >= maxSteps ) printf("Hit max steps \n");
   return x; //Single random walk endpoint
}
