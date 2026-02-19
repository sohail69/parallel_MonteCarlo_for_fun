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

/**************************************\
!
! Queries if a point is inside the
! domain
!
\**************************************/
// Returns true is if a point x is interior
// to a set of boundaries on a closed geometry
// false if it falls outside
template<typename real, size_t sdim, size_t edim>
FORCE_INLINE bool insideDomain(const Point<real,sdim> & x
                             , const std::vector<boundary<real,sdim,edim>> & boundaries)
{
//   REAL Theta = signedAngle<REAL,UINT>(x, boundaryDirichlet) + signedAngle<REAL,UINT>(x, boundaryNeumann);
   real Theta = 0.001;
   const real delta = 1e-4; // numerical tolerance
   return (std::abs(Theta-2.*M_PI) < delta); // boundary winds around x exactly once
}


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
real findClosestDistance(const boundary<real,sdim,edim> & DirchBC, const Point<real,sdim> & p0)
{
  Point<real,sdim> sfp;
  findClosestSurfacePoint(DirchBC, p0, sfp);
  return pointEuclideanDistance(p0, sfp);
}

// Find the boundary which
// is closest to the given
// point
template<typename real, size_t sdim, size_t edim>
void findClosestSurfacePoint(const boundary<real,sdim,edim> & DirchBC
                           , const Point<real,sdim>         & p0
                           , Point<real,sdim>               & sfp)
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
real findSilhouettePointDistance(const boundary<real,sdim,edim> & NeumBC, const Point<real,sdim> & p0)
{
  Point<real,sdim> sfp;
  findClosestSilhouettePoint(NeumBC, p0, sfp);
  return pointEuclideanDistance(p0, sfp);
};

// Find the boundary which
// is closest to the given
// point
template<typename real, size_t sdim, size_t edim>
void findClosestSilhouettePoint(const boundary<real,sdim,edim> & NeumBC
                              , const Point<real,sdim>         & p0
                              , Point<real,sdim>               & sfp)
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
real findClosestRayIntersectionPointDistance(const boundary<real,sdim,edim> & NeumBC, const Point<real,sdim> & p0)
{
  Point<real,sdim> sfp;
  findClosestSilhouettePoint(NeumBC, p0, sfp);
  return pointEuclideanDistance(p0, sfp);
};

// Find the boundary which
// is closest to the given
// point
template<typename real, size_t sdim, size_t edim>
void findClosestRayIntersectionPoint(const boundary<real,sdim,edim> & NeumBC
                                    , const Point<real,sdim>    & p0
                                    , Point<real,sdim>          & sfp)
{


};


/**************************************\
!
! Compute the star-shaped region
! radii
!
\**************************************/
/*
template<typename real>
void ComputeStarRadius(float & starRadius)
{
  if (firstStep && firstSphereRadius > 0.0f) {
    starRadius = firstSphereRadius;
  }else{
    // for problems with double-sided boundary conditions, flip the current
    // normal orientation if the geometry is front-facing
    flipNormalOrientation = false;
    if (walkSettings.solveDoubleSided && state.onNeumannBoundary) {
      if (state.prevDistance > 0.0f && state.prevDirection.dot(state.currentNormal) < 0.0f) {
        state.currentNormal *= -1.0f;
        flipNormalOrientation = true;
      }
    }

    if (walkSettings.stepsBeforeUsingMaximalSpheres <= state.walkLength) {
      starRadius = dirichletDist;
    } else {
      // NOTE: using dirichletDist as the maximum radius for the closest silhouette
      // query can result in a smaller than maximal star-shaped region: should ideally
      // use the distance to the closest visible Dirichlet point
      starRadius = queries.computeStarRadius(state.currentPt, walkSettings.minStarRadius
                                            ,dirichletDist, walkSettings.silhouettePrecision
                                            ,flipNormalOrientation);

      // shrink the radius slightly for numerical robustness---using a conservative
      // distance does not impact correctness
      if (walkSettings.minStarRadius <= dirichletDist) {
        starRadius = std::max(RADIUS_SHRINK_PERCENTAGE*starRadius, walkSettings.minStarRadius);
      }
    }
  }
};*/


/**************************************\
!
! Compute the Neumann boundary
! condition
!
\**************************************/
/*
if (!walkSettings.ignoreNeumannContribution) {
  // compute the non-zero Neumann contribution inside the star-shaped region;
  // define the Neumann value to be zero outside this region
  BoundarySample<DIM> neumannSample;
  for (int i = 0; i < DIM; i++) randNumsForNeumannSampling[i] = sampler.nextFloat();
  if (queries.sampleNeumann(state.currentPt, starRadius, randNumsForNeumannSampling, neumannSample)) {
    Vector<DIM> directionToSample = neumannSample.pt - state.currentPt;
    float distToSample = directionToSample.norm();
    float alpha = state.onNeumannBoundary ? 2.0f : 1.0f;
    bool estimateBoundaryNormalAligned = false;

    if (walkSettings.solveDoubleSided) {
      // normalize the direction to the sample, and flip the sample normal
      // orientation if the geometry is front-facing; NOTE: using a precision
      // parameter since unlike direction sampling, samples can lie on the same
      // halfplane as the current walk location
      directionToSample /= distToSample;
      if (flipNormalOrientation) {
        neumannSample.normal *= -1.0f;
        estimateBoundaryNormalAligned = true;
      } else if (directionToSample.dot(neumannSample.normal) < -walkSettings.silhouettePrecision) {
        bool flipNeumannSampleNormal = true;
        if (alpha > 1.0f) {
          // on concave boundaries, we want to sample back-facing neumann
          // values on front-facing geometry below the hemisphere, so we
          // avoid flipping the normal orientation in this case
          flipNeumannSampleNormal = directionToSample.dot(state.currentNormal) < -walkSettings.silhouettePrecision;
        }
        if (flipNeumannSampleNormal) {
          neumannSample.normal *= -1.0f;
          estimateBoundaryNormalAligned = true;
        }
      }
    }
    if (neumannSample.pdf > 0.0f && distToSample < starRadius && 
        !queries.intersectsWithNeumann(state.currentPt, neumannSample.pt, state.currentNormal, 
                                       neumannSample.normal, state.onNeumannBoundary, true))
    {
      float G = greensFn->evaluate(state.currentPt, neumannSample.pt);
      Value<T, DIM> h = walkSettings.solveDoubleSided ? 
                        pde.neumannDoubleSided(neumannSample.pt, estimateBoundaryNormalAligned) :             
                        pde.neumann(neumannSample.pt);
      state.totalNeumannContribution += state.throughput*alpha*G*h/neumannSample.pdf;
    }
  }
}*/

/**************************************\
!
! Neumann boundary silhouette
! intersection query
!
\**************************************/
// Checks for the intersection
// of the the point on the
// closest line
/*
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
}*/

