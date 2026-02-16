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
//#include "include/boundary/localVectorAlgebra.hpp"
//#include "include/boundary/IntersectionDetection2D.hpp"
#include "problemOperators/poissonGreenFunctionOperator2D.hpp"

//
// This version assumes that the number of device cores
// exceeds the number of 
//

int main(){
/******************************************************\
!  Sample usage region
\******************************************************/
  static const size_t nVars=1; //Number of variables being solved for
  static const size_t sdim=2;  //Number of spatial dimensions
  static const size_t edim=1;  //Dimension of boundary entities
/******************************************************\
\******************************************************/


  //The MPI communicator
  bool IS_MPI_ON=false;
  MPIComm mpiComm(IS_MPI_ON);

  //Problem base size
  const int ndev_cores=1000;
  const int nWalks = 65536;              //Total number of Monte Carlo samples
  const int s = 128;                     //Image length-width
  const int nSize = s*s;                 //Total image size
  const double dx= 1.0/double(s);        //Image increment

  //Generate the blockMesh from
  //the base size (2D-Quad) and
  //a 2D polyline boundary
  blockHyperMeshData<double,int,sdim> BHMeshData;
  simpleBlockMeshBuild<double,int,2>(dx, s, BHMeshData);
//  std::vector<Polyline2D<double>>  bcDirch, bcNeum;
//  simple2DBoundary<double>(bcDirch,  bcNeum);


  //MPI global partitioning (Geometric paritioning of the domain)
  WoStr_Partition<int> wostr_part;
  partitionProblem<double,int,2>(mpiComm, BHMeshData, nWalks, wostr_part);


  //Iterators and other thread local
  //variables
  int rnd_seed=0, threadID=0, parentMPI_ID=0, nAccums=0, AccumStart=0;
  int nTotAccumulators=0;

  //Generate a vector for storing
  //the local solution data

/*
  std::vector<int>    InDomainFlag(MPI_lSize);
  std::vector<double> u_sol(MPI_lSize);



  nTotAccumulators=MPI_lSize;
  std::vector<double> accumulator(nTotAccumulators);
  double zero(0.00);

  #pragma omp target
  #pragma omp private(threadID, parentMPI_ID, rnd_seed, nAccums, AccumStart)
  {
    for(int I=Istart; I<IEnd; I++){

    }

    threadID = omp_get_thread_num();
    parentMPI_ID = FindAccumParentID<int>(threadID, nTotAccumulators, MPI_lSize);

    //First find all initial points
    //that fall inside the domain
    if(threadID < MPI_lSize){
      I = threadID + MPI_ITstart;
      Vec2D<double> x0;
      BHMeshPoint<double,int,2>(x0.data(), I, BHMeshData);
      InDomainFlag[threadID] = (insideDomain<double,int>(x0, bcDirch, bcNeum) ? 1:0);
    }
    #pragma omp barrier

    //Sample the Monte-Carlo problem
    //on the accumulators
    accumulator[threadID] = double(0.00);
    for(int IWalk=0; IWalk<nWalksPerThread; IWalk++){
      I = parentMPI_ID + MPI_ITstart;
      rnd_seed = 0;
      Vec2D<double> x0;
      BHMeshPoint<double,int,2>(x0.data(), I, BHMeshData);
      double val = poissonWalk<double,int>(x0, bcDirch, bcNeum, lines<double>, rnd_seed);
*/
/******************************************************\
!  Sample usage region
\******************************************************/

//Point<double, sdim> x0;
//Point<real,sdim> x = WoStr_point<double,XORSHIFT256_rngData,sdim,edim>
//                                   (x0,Dirichlet,Neumann,lines<double>,XORSHIFT256_rngUpdate, rngData);
/******************************************************\
\******************************************************/
/*
      accumulator[threadID] += (InDomainFlag[parentMPI_ID]==1)? val : zero;
    }

    //Aggregate all the accumulators on
    //the solution vector and average the
    //solution
    #pragma omp barrier
    if(threadID < MPI_lSize){
      nAccums    = FindLocalNumberOfAccums<int>(threadID, nTotAccumulators, MPI_lSize);
      AccumStart = FindPartIDFirstLocalAccumPos<int>(threadID, nTotAccumulators, MPI_lSize);
      u_sol[threadID] = double(0.00);
      for(I=AccumStart; I<(AccumStart+nAccums); I++){
        u_sol[threadID] += accumulator[I];
      }
      u_sol[threadID] /= double(nWalks);
    }
  }
*/
  // Printing some extra data to
  // console (to check for if
  // program launched correctly)
  if(mpiComm.getProcID() == 0) std::cout << "Iterators" << std::endl;
  std::cout << std::setw(10) << wostr_part.procID
            << std::setw(10) << wostr_part.nProcs
            << std::setw(10) << wostr_part.mpi_Istart
            << std::setw(10) << wostr_part.mpi_Iend
            << std::setw(10) << wostr_part.mpi_lsize  << std::endl;



  // Output the data into
  // a file (the original
  // used CSV's, IO tbd)


  //clean-up
  //InDomainFlag.clear();
  //accumulator.clear();
  //u_sol.clear();
  MPI_Finalize();
  return 0;
}
