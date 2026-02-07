#pragma once
#include <random>
#include <cstdio>
#include <vector>
#include <functional>
#include "../include/globalMacros.hpp"
#include "../include/blockHyperMesh.hpp"
#include "../include/Random/RNG.hpp"
#include "../include/localVectorAlgebra.hpp"
#include "../include/IntersectionDetection2D.hpp"

//Simple block mesh data
//for uniform grids
template<typename REAL, typename UINT, UINT dim>
void simpleBlockMeshBuild(const REAL & dx, const UINT & size, blockHyperMeshData<REAL,UINT,dim> & BHMeshData){
  for(UINT I=0; I<dim; I++){
    BHMeshData.sizes[I]  = size;
    BHMeshData.offset[I] = REAL(0.00);
    for(UINT J=0; J<dim; J++){
      BHMeshData.dx[I*dim + J] = ((I==J) ? dx : REAL(0.00));
    }
  }
}

//Simple imbedded boundary
//description data
template<typename REAL>
void simple2DBoundary(std::vector<Polyline2D<REAL>> & bcDirch
                    , std::vector<Polyline2D<REAL>> & bcNeum)
{
  // for simplicity, in this code we assume that the Dirichlet and Neumann
  // boundary polylines form a collection of closed polygons (possibly with holes),
  // and are given with consistent counter-clockwise orientation
  bcDirch.clear();
  bcNeum.clear();
  bcDirch.push_back({{ Vec2D<REAL>({0.2, 0.2}), Vec2D<REAL>({0.6, 0.0}), Vec2D<REAL>({1.0, 0.2}) }});
  bcDirch.push_back({{ Vec2D<REAL>({1.0, 1.0}), Vec2D<REAL>({0.6, 0.8}), Vec2D<REAL>({0.2, 1.0}) }});
  bcNeum.push_back({{  Vec2D<REAL>({1.0, 0.2}), Vec2D<REAL>({0.8, 0.6}), Vec2D<REAL>({1.0, 1.0}) }});
  bcNeum.push_back({{  Vec2D<REAL>({0.2, 1.0}), Vec2D<REAL>({0.0, 0.6}), Vec2D<REAL>({0.2, 0.2}) }});
};


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


template<typename REAL, typename UINT>
FORCE_INLINE REAL poissonWalk(Vec2D<REAL> x0,                                  // evaluation point
                              std::vector<Polyline2D<REAL>> boundaryDirichlet, // absorbing part of the boundary
                              std::vector<Polyline2D<REAL>> boundaryNeumann,   // reflecting part of the boundary
                              std::function<REAL(Vec2D<REAL>)> g,              // Dirichlet boundary values
                              UINT seedVal)
{
   const REAL eps = 0.0001;     // stopping tolerance
   const REAL rMin = 0.0001;    // minimum step size
   const UINT nWalks = 65536;   // number of Monte Carlo samples
   const UINT maxSteps = 65536; // maximum walk length


   uint32_t rqd_seed = 0UL + uint32_t(seedVal);
   REAL sum = 0.0; // running sum of boundary contributions
//   for( UINT i = 0; i < nWalks; i++ ) {
      Vec2D<REAL> x = x0;        // start walk at the evaluation point
      Vec2D<REAL> n={ 0.0, 0.0 };// assume x0 is an interior point, and has no normal
      bool onBoundary = false;   // flag whether x is on the interior or boundary

      REAL r, dDirichlet, dSilhouette; // radii used to define star shaped region
      UINT steps = 0;
      for(steps=0; (dDirichlet > eps) && (steps < maxSteps); steps++){
         // loop until the walk hits the Dirichlet boundary
         // compute the radius of the largest star-shaped region
         dDirichlet = distancePolylines<REAL,UINT>( x, boundaryDirichlet );
         dSilhouette = silhouetteDistancePolylines<REAL,UINT>( x, boundaryNeumann );
         r = std::max( rMin, std::min(dDirichlet,dSilhouette) );

         // intersect a ray with the star-shaped region boundary
         REAL theta = random<double>( -M_PI, M_PI );
//         rqd_seed = randqd_uint32(rqd_seed);
//         REAL theta = my_random<double>( -M_PI, M_PI, rqd_seed);

         if( onBoundary ) { // sample from a hemisphere around the normal
            theta = theta/REAL(2.) + angleOf2DVec<REAL>(n);
         }
         Vec2D<REAL> v{ cos(theta), sin(theta) }; // unit ray direction
         x = intersectPolylines<REAL,UINT>( x, v, r, boundaryNeumann, n, onBoundary );
         steps++;
      } //stop if we hit the Dirichlet boundary, or the walk is too long
      if( steps >= maxSteps ) printf("Hit max steps \n");
      sum += g(x); // accumulate contribution of the boundary value
//   }
   return sum; // Single random walk
}
