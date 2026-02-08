// (Slow) implementation of the "Walk on Stars" algorithm for Laplace equations.
// Corresponds to the estimator given in Equation 18 of Sawhney et al,
// "Walk on Stars: A Grid-Free Monte Carlo Method for PDEs with Neumann Boundary
// Conditions" (2023), assuming no source term and zero-Neumann conditions.
// NOTE: this code makes a few shortcuts for the sake of code brevity; may
// be more suitable for tutorials than for production code/evaluation.
#include <algorithm>
#include <array>
#include <complex>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <cstdio>
#include <random>

//Parallel libs
#include <omp.h>
#include "include/MPIcomm.hpp"
#include "include/partitioner.hpp"

//Geometry libs
#include "include/blockHyperMesh.hpp"
#include "include/boundary/localVectorAlgebra.hpp"
#include "include/boundary/IntersectionDetection2D.hpp"
#include "problemOperators/poissonGreenFunctionOperator2D.hpp"

//
// This version assumes that the number of device cores
// exceeds the number of 
//

int main(){
  //The MPI communicator
  bool IS_MPI_ON=false;
  MPIComm mpiComm(IS_MPI_ON);

  //Problem base size
  const int nWalks = 65536;              //Total number of Monte Carlo samples
  const int nWalksPerThread = 65536/542; //Number of Monte Carlo samples per thread
  const int s = 128;                     //Image length-width
  const int nSize = s*s;                 //Total image size
  const double dx= 1.0/double(s);        //Image increment

  //Generate the blockMesh from
  //the base size (2D-Quad) and
  //a 2D polyline boundary
  blockHyperMeshData<double,int,2> BHMeshData;
  std::vector<Polyline2D<double>> bcDirch, bcNeum;
  simpleBlockMeshBuild<double,int,2>(dx, s, BHMeshData);
  simple2DBoundary<double>(bcDirch,  bcNeum);

  //Partitioning the data on multiple
  //levels
  int MPI_lSize, MPI_ITstart, MPI_ITend;
  int ndev_cores=1000;

  //MPI global partitioning (Geometric paritioning of the domain)
  MPI_ITstart = firstIterator<int>(mpiComm.getProcID(), mpiComm.getNProcs(), nSize);
  MPI_ITend   = lastIterator<int>( mpiComm.getProcID(), mpiComm.getNProcs(), nSize);
  MPI_lSize   = MPI_ITend - MPI_ITstart;

  //Generate a vector for storing
  //the local solution data
  std::vector<int>    InDomainFlag(MPI_lSize);
  std::vector<double> u_sol(MPI_lSize);
  std::vector<double> accumulator(ndev_cores);
  double zero(0.00);

////////  #pragma omp shared() private(i)
  #pragma omp target
  {
    int threadID = omp_get_thread_num();
    int ParentAccumulator = FindAccumParentID<int>(threadID, ndev_cores, MPI_lSize);

    //First find all points that fall
    //inside the domain
    if(threadID < MPI_lSize){
      int I = threadID + MPI_ITstart;
      Vec2D<double> x0;
      BHMeshPoint<double,int,2>(x0.data(), I, BHMeshData);
      InDomainFlag[threadID] = (insideDomain<double,int>(x0, bcDirch, bcNeum) ? 1:0);
    }
    #pragma omp barrier

    //Sample the Monte-Carlo problem
    //on the accumulators
    accumulator[threadID] = double(0.00);
    for(int IWalk=0; IWalk<nWalksPerThread; IWalk++){
      int I = ParentAccumulator + MPI_ITstart;
      int rnd_seed = 0;
      Vec2D<double> x0;
      BHMeshPoint<double,int,2>(x0.data(), I, BHMeshData);
      double val = poissonWalk<double,int>(x0, bcDirch, bcNeum, lines<double>, rnd_seed);
      accumulator[threadID] += (InDomainFlag[ParentAccumulator]==1)? val : zero;
    }

    //Aggregate all the accumulators on
    //the solution vector
    #pragma omp barrier
    if(threadID < MPI_lSize){
      int nAccums  = FindLocalNumberOfAccums<int>(threadID, ndev_cores, MPI_lSize);
      int devStart = FindPartIDFirstLocalAccumPos<int>(threadID, ndev_cores, MPI_lSize);
      u_sol[threadID] = double(0.00);
      for(int I=devStart; I<(devStart+nAccums); I++){
        u_sol[threadID] += accumulator[I];
      }
      u_sol[threadID] /= double(nWalks);
    }
  }


  // Printing some extra data to
  // console (to check for if
  // program launched correctly)
  if(mpiComm.getProcID() == 0) std::cout << "Iterators" << std::endl;
  std::cout << std::setw(10) << mpiComm.getProcID()
            << std::setw(10) << mpiComm.getNProcs() 
            << std::setw(10) << MPI_lSize
            << std::setw(10) << MPI_ITstart
            << std::setw(10) << MPI_ITend   << std::endl;


  // Output the data into
  // a file (the original
  // used CSV's, IO tbd)


  //clean-up
  InDomainFlag.clear();
  accumulator.clear();
  u_sol.clear();
  MPI_Finalize();
  return 0;
}
