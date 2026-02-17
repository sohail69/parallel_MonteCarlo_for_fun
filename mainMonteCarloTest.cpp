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
#include "include/templatedGeometry/geometricPrimitives.hpp"
#include "problemOperators/poissonGreenFunctionOperator2D.hpp"

//boundary libs
#include "include/boundary/boundary.hpp"
#include "include/boundary/boundaryQueries.hpp"
//#include "include/boundary/localVectorAlgebra.hpp"
//#include "include/boundary/IntersectionDetection2D.hpp"

/*****************************************\
! A MPI-OpenMP parallel implementation
! of the walk on stars algorithm used
! to solve PDE
\*****************************************/
int main(){
  // Problem size
  static const size_t nVars=1; //Number of variables being solved for
  static const size_t sdim=2;  //The spatial dimension
  static const size_t edim=1;  //The boundary entities dimension

  // The MPI communicator
  bool IS_MPI_ON=false;
  MPIComm mpiComm(IS_MPI_ON);

  // Problem base size
  const int nWalks = 65536;              //Total number of Monte Carlo samples
  const int s = 128;                     //Image length-width
  const double dx= 1.0/double(s);        //Image increment

  /*****************************************\
  ! Read in any boundary data build the
  ! initial geometric sample space and
  ! partition it among the MPI processes
  \*****************************************/
  // Generate the blockMesh from the base
  // size (2D-Quad) and 2D polyline boundary
  blockHyperMeshData<double,int,sdim> BHMeshData;
  simpleBlockMeshBuild<double,int,sdim>(dx, s, BHMeshData);
  std::vector<boundary<double,sdim,edim>> boundaries;

  // Setup the boundaries
  boundary<double,sdim,edim> bcDirch, bcNeum;
  boundaries.push_back(bcDirch);
  boundaries.push_back(bcNeum);
//  std::vector<Polyline2D<double>>  bcDirch, bcNeum;
//  simple2DBoundary<double>(bcDirch,  bcNeum);

  // MPI global partitioning (Geometric paritioning of the domain)
  WoStr_Partition<int> wostr_part;
  partitionProblem<double,int,2>(mpiComm, BHMeshData, nWalks, wostr_part);

  // Generate a vector for storing
  // the local solution data
  std::vector<int>    InDomainFlag(wostr_part.mpi_lsize);
  std::vector<double> u_sol(wostr_part.mpi_lsize);
  std::vector<double> accumulator(wostr_part.nAccumsPerMPI);

  // First find all initial points
  // that fall inside the domain
  for(unsigned iPos=0; iPos< wostr_part.mpi_lsize; iPos++){
    Point<double,sdim> x0;
    BHMeshPoint<double,int,sdim>(x0.data(), iPos + wostr_part.mpi_Istart, BHMeshData);
    InDomainFlag[iPos] = insideDomain<double,sdim,edim>(x0,boundaries) ? 1:0;
  }

  /*****************************************\
  ! Run the Monte-Carlo random walks for
  ! running the PDE problem
  \*****************************************/
/*double zero(0.00);
  #pragma omp target map(tofrom:C[0:Ndim*Mdim]) map(to:B[0:Pdim*Mdim],A[0:Ndim*Pdim])
  #pragma omp private(threadID, parentMPI_ID, rnd_seed, nAccums, AccumStart)
  {
    threadID = omp_get_thread_num();
    parentMPI_ID = FindAccumParentID<int>(threadID, nTotAccumulators, MPI_lSize);

    //Sample the Monte-Carlo problem
    //and aggregate on the accumulators
    for(unsigned iAccum=0; iAccum<(wostr_part.nAccumsPerThread); iAccum++){
      for(unsigned iWalks=0; iWalks<; iWalks++){
        unsigned jAccum = threadID*(wostr_part.nAccumsPerThread) + iAccum;
        unsigned iPart = 
        Point<double, sdim> x0, x;
        BHMeshPoint<double,int,2>(x0.data(), I, BHMeshData);
        accumulator[jAccum] += g(x);
      }
    }

    for(int iSubAccum=0; iSubAccum<nAccum; iSubAccum++){
      accumulator[threadID] = double(0.00);
      for(int IWalk=0; IWalk<nWalksPerThread; IWalk++){
        I = parentMPI_ID + MPI_ITstart;
        rnd_seed = 0;
        Point<double, sdim> x0;
        BHMeshPoint<double,int,2>(x0.data(), I, BHMeshData);
        Point<real,sdim> x = WoStr_point<double,XORSHIFT256_rngData,sdim,edim>
                                           (x0,Dirichlet,Neumann,lines<double>,XORSHIFT256_rngUpdate, rngData);
        accumulator[threadID] += g(x);
      }
    }
  }*/

  /*****************************************\
  ! Aggregate the solution partial
  ! accumulators on to the solution vector
  \*****************************************/
  // First find all initial points
  // that fall inside the domain

/*
  for(unsigned iPos=0; iPos< wostr_part.mpi_lsize; iPos++){
    unsigned nAccums  = wostr_part.nAccumsPerPart;
    unsigned AccStart = iPos * wostr_part.nAccumsPerPart;
    u_sol[iPos] = double(0.00);
    for(unsigned iAccum=AccStart; iAccum<(AccStart+nAccums); iAccum++) u_sol[iPos] += accumulator[iAccum];
    u_sol[iPos] /= double(nWalks);
  }*/


  /*****************************************\
  ! Output the solution
  \*****************************************/
  // Printing some extra data to
  // console (to check for if
  // program launched correctly)
  if(mpiComm.getProcID() == 0) std::cout << "Iterators" << std::endl;
  std::cout << std::setw(10) << wostr_part.procID
            << std::setw(10) << wostr_part.nProcs
            << std::setw(10) << wostr_part.nAccumsPerMPI
            << std::setw(10) << wostr_part.mpi_Istart
            << std::setw(10) << wostr_part.mpi_Iend
            << std::setw(10) << wostr_part.mpi_lsize  << std::endl;

  // Output the data into
  // a file (the original
  // used CSV's, IO tbd)


  //clean-up
  accumulator.clear();
  u_sol.clear();
  InDomainFlag.clear();
  MPI_Finalize();
  return 0;
}
