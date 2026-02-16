
// estimates only the solution of the given PDE at the input point
void estimateSolution(const PDE<T, DIM>& pde,
                      const WalkSettings<T>& walkSettings,
                      int nWalks,
                      int nRecursiveWalks,
                      SamplePoint<T, DIM>& samplePt) const
{
  // initialize statistics if there are no previous estimates
  bool hasPrevEstimates = samplePt.statistics != nullptr;
  if (!hasPrevEstimates) samplePt.statistics = std::make_shared<SampleStatistics<T, DIM>>(walkSettings.initVal);

  // check if the sample pt is on the Dirichlet boundary (and does not require recursive walks)
  if (samplePt.type == SampleType::OnDirichletBoundary && nRecursiveWalks <= 0) {
    if (!hasPrevEstimates) {
      // record the known boundary value
      Value<T, DIM> totalContribution = walkSettings.initVal;
      if (!walkSettings.ignoreDirichletContribution){
        totalContribution = walkSettings.solveDoubleSided ?
                            pde.dirichletDoubleSided(samplePt.pt, samplePt.estimateBoundaryNormalAligned) :
                            pde.dirichlet(samplePt.pt);
      }
				
      // update statistics and set the first sphere radius to 0
      samplePt.statistics->addSolutionEstimate(totalContribution);
      samplePt.firstSphereRadius = 0.0f;
    }
    // no need to run any random walks
    return;

  }else if(samplePt.dirichletDist <= walkSettings.epsilonShell){
    // run just a single walk since the sample pt is inside the epsilon shell
    nWalks = 1;
  }

  // for problems with double-sided boundary conditions, initialize the direction
  // of approach for walks, and flip the current normal orientation if the geometry
  // is front-facing
  Vector<DIM> currentNormal = samplePt.normal;
  Vector<DIM> prevDirection = samplePt.normal;
  float prevDistance = std::numeric_limits<float>::max();
  bool flipNormalOrientation = false;

  if (walkSettings.solveDoubleSided && samplePt.type == SampleType::OnNeumannBoundary) {
    if (samplePt.estimateBoundaryNormalAligned) {
      currentNormal *= -1.0f;
      prevDirection *= -1.0f;
      flipNormalOrientation = true;
    }
  }

  // precompute the first sphere radius for all walks
  if (!hasPrevEstimates) {
    if(samplePt.dirichletDist > walkSettings.epsilonShell && walkSettings.stepsBeforeUsingMaximalSpheres != 0){
      // compute the star radius; NOTE: using dirichletDist as the maximum radius for
      // the closest silhouette query can result in a smaller than maximal star-shaped
      // region: should ideally use the distance to the closest visible Dirichlet point
      float starRadius = queries.computeStarRadius(samplePt.pt
                                                 , walkSettings.minStarRadius
                                                 , samplePt.dirichletDist
                                                 , walkSettings.silhouettePrecision
                                                 , flipNormalOrientation);

      // shrink the radius slightly for numerical robustness---using a conservative
      // distance does not impact correctness
      if (walkSettings.minStarRadius <= samplePt.dirichletDist) {
        starRadius = std::max(RADIUS_SHRINK_PERCENTAGE*starRadius, walkSettings.minStarRadius);
      }
      samplePt.firstSphereRadius = starRadius;
    }else{
      samplePt.firstSphereRadius = samplePt.dirichletDist;
    }
  }

  std::shared_ptr<ProductEstimate<T>> productEstimate;
  if(walkSettings.solutionWeightedDifferentialBatchSize > 0){
    productEstimate = std::make_shared<ProductEstimate<T>>(
    walkSettings.solutionWeightedDifferentialBatchSize, walkSettings.initVal);
  }

  // perform random walks
  for (int w = 0; w < nWalks; w++) {
    // initialize the greens function
    std::unique_ptr<GreensFnBall<DIM>> greensFn = nullptr;
    if(pde.absorption > 0.0f && walkSettings.stepsBeforeApplyingTikhonov == 0){
      greensFn = std::make_unique<YukawaGreensFnBall<DIM>>(pde.absorption);
    }else{
      greensFn = std::make_unique<HarmonicGreensFnBall<DIM>>();
    }

    // initialize the walk state
    WalkState<T, DIM> state(samplePt.pt
                          , currentNormal
                          , prevDirection
                          , prevDistance,
                          , 1.0f
                          , samplePt.type == SampleType::OnNeumannBoundary
                          , 0
                          , walkSettings.initVal);

    // perform walk
    WalkCompletionCode code = walk(pde
                                 , walkSettings
                                 , samplePt.dirichletDist
                                 , samplePt.firstSphereRadius
                                 , flipNormalOrientation
                                 , samplePt.sampler
                                 , greensFn
                                 , state);

    if ((code == WalkCompletionCode::ReachedDirichletBoundary ||
        (code == WalkCompletionCode::TerminatedWithRussianRoulette) ||
        (code == WalkCompletionCode::ExceededMaxWalkLength && getTerminalContribution))
    {
      RecursiveBoundaryData<T> recursiveBoundaryData(walkSettings.initVal);
      if (code == WalkCompletionCode::ReachedDirichletBoundary){
        estimateRecursiveBoundaryData(pde
                                    , walkSettings
                                    , nRecursiveWalks
                                    , state.currentPt
                                    , recursiveBoundaryData);
      }

      // compute the walk contribution
      setTerminalContribution(code, pde, walkSettings, state, recursiveBoundaryData);
      Value<T, DIM> totalContribution = state.throughput*state.terminalContribution
                                      + state.totalNeumannContribution
                                      + state.totalSourceContribution;
				
      // update statistics
      samplePt.statistics->addSolutionEstimate(totalContribution);
      samplePt.statistics->addWalkLength(state.walkLength);

      // compute unbiased estimate of the product u * u'
      if (productEstimate) {
        productEstimate->recordContributions(w, totalContribution.data, totalContribution.differential);
        if (productEstimate->isCacheFull(w)) {
          samplePt.statistics->addSolutionWeightedDifferential(productEstimate->compute(walkSettings.initVal));
          productEstimate->reset(walkSettings.initVal);
        }
      }
    }
  }
}
