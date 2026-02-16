// performs a single reflecting random walk starting at the input point
WalkCompletionCode walk(const PDE<T, DIM>& pde
                      , const WalkSettings<T>& walkSettings,
                      , float dirichletDist
                      , float firstSphereRadius
                      , bool flipNormalOrientation
                      , pcg32& sampler
                      , std::unique_ptr<GreensFnBall<DIM>>& greensFn
                      , WalkState<T, DIM>& state) const
{
  // recursively perform a random walk till it reaches the Dirichlet boundary
  bool firstStep = true;
  float randNumsForNeumannSampling[DIM];

  while (dirichletDist > walkSettings.epsilonShell) {
    // compute the star radius
    float starRadius;
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
        if(walkSettings.stepsBeforeUsingMaximalSpheres <= state.walkLength){
          starRadius = dirichletDist;
        }else{
          // NOTE: using dirichletDist as the maximum radius for the closest silhouette
          // query can result in a smaller than maximal star-shaped region: should ideally
          // use the distance to the closest visible Dirichlet point
          starRadius = queries.computeStarRadius(state.currentPt
                                               , walkSettings.minStarRadius
                                               , dirichletDist
                                               , walkSettings.silhouettePrecision
                                               , flipNormalOrientation);
          // shrink the radius slightly for numerical robustness---using a conservative
          // distance does not impact correctness
          if (walkSettings.minStarRadius <= dirichletDist) {
            starRadius = std::max(RADIUS_SHRINK_PERCENTAGE*starRadius, walkSettings.minStarRadius);
          }
        }
      }
      // update the ball center and radius
      greensFn->updateBall(state.currentPt, starRadius);
      // sample a direction uniformly
      Vector<DIM> direction = sampleUnitSphereUniform<DIM>(sampler);
      // perform hemispherical sampling if on the Neumann boundary, which cancels
      // the alpha term in our integral expression
      if(state.onNeumannBoundary && state.currentNormal.dot(direction) > 0.0f){
        direction *= -1.0f;
      }
      // check if there is an intersection with the Neumann boundary along the ray:
      // currentPt + starRadius * direction
      IntersectionPoint<DIM> intersectionPt;
      bool intersectedNeumann = queries.intersectWithNeumann(state.currentPt
                                                           , state.currentNormal
                                                           , direction
                                                           , starRadius
                                                           , state.onNeumannBoundary
                                                           , intersectionPt);
      // check if there is no intersection with the Neumann boundary
      if(!intersectedNeumann){
        // apply small offset to the current pt for numerical robustness if it on
        // the Neumann boundary---the same offset is applied during ray intersections
        Vector<DIM> currentPt = state.onNeumannBoundary ?
                                queries.offsetPointAlongDirection(state.currentPt, -state.currentNormal):
                                state.currentPt;
        // set intersectionPt to a point on the spherical arc of the ball
        intersectionPt.pt = currentPt + starRadius*direction;
        intersectionPt.dist = starRadius;
      }
      if (!walkSettings.ignoreNeumannContribution) {
        // compute the non-zero Neumann contribution inside the star-shaped region;
        // define the Neumann value to be zero outside this region
        BoundarySample<DIM> neumannSample;
        for (int i = 0; i < DIM; i++) randNumsForNeumannSampling[i] = sampler.nextFloat();
        if(queries.sampleNeumann(state.currentPt, starRadius, randNumsForNeumannSampling, neumannSample)){
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
            if(flipNormalOrientation){
              neumannSample.normal *= -1.0f;
              estimateBoundaryNormalAligned = true;
            } else if (directionToSample.dot(neumannSample.normal) < -walkSettings.silhouettePrecision) {
              bool flipNeumannSampleNormal = true;
            if(alpha > 1.0f){
              // on concave boundaries, we want to sample back-facing neumann
              // values on front-facing geometry below the hemisphere, so we
              // avoid flipping the normal orientation in this case
             flipNeumannSampleNormal = directionToSample.dot(state.currentNormal) < (-walkSettings.silhouettePrecision);
           }
            if(flipNeumannSampleNormal){
              neumannSample.normal *= -1.0f;
              estimateBoundaryNormalAligned = true;
            }
          }
        }
        if(  (neumannSample.pdf > 0.0f)
          && (distToSample < starRadius)
          && !queries.intersectsWithNeumann(state.currentPt
                                          , neumannSample.pt
                                          , state.currentNormal
                                          , neumannSample.normal
                                          , state.onNeumannBoundary
                                          , true)
        )
        {
          float G = greensFn->evaluate(state.currentPt, neumannSample.pt);
          Value<T, DIM> h = walkSettings.solveDoubleSided ?
                             pde.neumannDoubleSided(neumannSample.pt, estimateBoundaryNormalAligned) :
                             pde.neumann(neumannSample.pt);
          state.totalNeumannContribution += state.throughput*alpha*G*h/neumannSample.pdf;
        }
      }
    }
    if(!walkSettings.ignoreSourceContribution) {
      // compute the source contribution inside the star-shaped region;
      // define the source value to be zero outside this region
      float sourcePdf;
      Vector<DIM> sourcePt = greensFn->sampleVolume(direction, sampler, sourcePdf);
      if (greensFn->r <= intersectionPt.dist) {
        // NOTE: hemispherical sampling causes the alpha term to cancel when
        // currentPt is on the Neumann boundary; in this case, the green's function
        // norm remains unchanged even though our domain is a hemisphere;
        // for double-sided problems in watertight domains, both the current pt
        // and source pt lie either inside or outside the domain by construction
        Value<T, DIM> sourceContribution = greensFn->norm()*pde.source(sourcePt);
        state.totalSourceContribution += state.throughput*sourceContribution;
      }
    }

    // update walk position
    state.prevDistance = intersectionPt.dist;
    state.prevDirection = direction;
    state.currentPt = intersectionPt.pt;
    state.currentNormal = intersectionPt.normal; // NOTE: stale unless intersectedNeumann is true
    state.onNeumannBoundary = intersectedNeumann;

    // check if the current pt lies outside the domain; for interior problems,
    // this tests for walks that escape due to numerical error
    if (!state.onNeumannBoundary && queries.outsideBoundingDomain(state.currentPt)) {
      if (walkSettings.printLogs) {
        std::cout << "Walk escaped domain!" << std::endl;
      }
      return WalkCompletionCode::EscapedDomain;
    }

    // update the walk throughput and use russian roulette to decide whether
    // to terminate the walk
    state.throughput *= greensFn->directionSampledPoissonKernel(state.currentPt);
    if (state.throughput < walkSettings.russianRouletteThreshold) {
      float survivalProb = state.throughput/walkSettings.russianRouletteThreshold;
      if (survivalProb < sampler.nextFloat()) {
        state.throughput = 0.0f;
        return WalkCompletionCode::TerminatedWithRussianRoulette;
      }
      state.throughput = walkSettings.russianRouletteThreshold;
    }

    // update the walk length and break if the max walk length is exceeded
    state.walkLength++;
    if (state.walkLength > walkSettings.maxWalkLength) {
      if (walkSettings.printLogs && !getTerminalContribution) {
        std::cout << "Maximum walk length exceeded!" << std::endl;
      }
      return WalkCompletionCode::ExceededMaxWalkLength;
    }

    // check whether to start applying Tikhonov regularization
    if (pde.absorption > 0.0f && walkSettings.stepsBeforeApplyingTikhonov == state.walkLength) {
      greensFn = std::make_unique<YukawaGreensFnBall<DIM>>(pde.absorption);
    }

    // compute the distance to the dirichlet boundary
    dirichletDist = queries.computeDistToDirichlet(state.currentPt, false);
    firstStep = false;
  }
  return WalkCompletionCode::ReachedDirichletBoundary;
}
