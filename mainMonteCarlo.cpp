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

#include "include/MPIcomm.hpp"
#include "include/partitioner.hpp"
#include "include/blockHyperMesh.hpp"
#include "include/localVectorAlgebra.hpp"
#include "include/IntersectionDetection2D.hpp"
#include "problemOperators/poissonGreenFunctionOperator2D.hpp"

int main(){
  //The MPI communicator
  bool IS_MPI_ON=false;
  MPIComm mpiComm(IS_MPI_ON);

  //Problem base size
  const int nWalksPerThread = 65536; // number of Monte Carlo samples per thread
  const int s = 128;                 //Image length-width
  const int nSize = s*s;             //Total image size
  const double dx= 1.0/double(s);    //Image increment

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
  int dev_lSize, dev_ITstart, dev_ITend;

  //MPI global partitioning (Geometric paritioning of the domain)
  MPI_ITstart = firstIterator<int>(mpiComm.getProcID(), mpiComm.getNProcs(), nSize);
  MPI_ITend   = lastIterator<int>( mpiComm.getProcID(), mpiComm.getNProcs(), nSize);
  MPI_lSize   = MPI_ITend - MPI_ITstart;

  //device local partitioning
  int dev_id=0, ndev_cores=1, It1D=0, I=0;

  //Generate a vector for storing
  //the local solution data
  std::vector<int> InDomainFlag;   InDomainFlag.resize(MPI_lSize);
  std::vector<double> u_sol;       u_sol.resize(MPI_lSize);
  std::vector<int> accumMapping;
  std::vector<double> accumulator;

  //Get the partitioned size
  //of the accumulator
  int LocalAccumSize = nWalks*dev_lSize;
  accumulator.resize(LocalAccumSize);
  accumMapping.resize(LocalAccumSize);
  double zero(0.00);

  //Mapping for addition of the accumulators
  for(dev_id=0; dev_id<ndev_cores; dev_id++){
    dev_ITstart = firstIterator<int>(dev_id, ndev_cores, LocalAccumSize);
    dev_ITend   = lastIterator<int>( dev_id, ndev_cores, LocalAccumSize);
    dev_lSize   = dev_ITend - dev_ITstart;
    for(int I=dev_ITstart; I <= dev_ITend; I++){ 
      accumMapping[I] = I
    }
  }

//  #pragma omp default(shared) private(dev_lSize, dev_ITstart, dev_ITend, dev_id, It1D, I)
//  #pragma omp target
  {
    dev_ITstart = firstIterator<int>(dev_id, ndev_cores, LocalAccumSize);
    dev_ITend   = lastIterator<int>( dev_id, ndev_cores, LocalAccumSize);
    dev_lSize   = dev_ITend - dev_ITstart;

    printf("dev_ITstart :  %d   dev_ITend :  %d   dev_lSize :  %d \n", dev_ITstart, dev_ITend, dev_lSize);

    int ValMax=0, ValMin=0;
    for(It1D=dev_ITstart; It1D <= dev_ITend; It1D++){
      I = It1D + MPI_ITstart;
      Vec2D<double> x0;
      BHMeshPoint<double,int,2>(x0.data(), I, BHMeshData);
      InDomainFlag[It1D] = (insideDomain<double,int>(x0, bcDirch, bcNeum) ? 1:0);
    }

    for(It1D=dev_ITstart; It1D<=dev_ITend; It1D++){
      I = It1D + MPI_ITstart;
      int rnd_seed = 0;
      Vec2D<double> x0;
      BHMeshPoint<double,int,2>(x0.data(), I, BHMeshData);
      u_sol[It1D]=(InDomainFlag[It1D]==1)?solve<double,int>(x0, bcDirch, bcNeum, lines<double>, rnd_seed):zero;
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
  std::ofstream out( "out"+std::to_string(mpiComm.getProcID())+".csv" );
  for(int I=0; I<s; I++){
    for(int J=0; J<s; J++){
      int Iters[2] = {I,J};
      int Iter1D =  BHMeshFwdIter<double,int,2>(Iters, BHMeshData);
      out <<  ((InDomainFlag[Iter1D] == 1)? 1234:0);
      if( J < s-1 ) out << ",";
    }
    out << std::endl;
  }
  out.close();

  //clean-up
  InDomainFlag.clear();
  u_sol.clear();
  MPI_Finalize();
  return 0;
}
