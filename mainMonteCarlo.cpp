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

#include "include/MPIcomm.hpp"
#include "include/partitioner.hpp"
#include "include/blockHyperMesh.hpp"
#include "include/localVectorAlgebra.hpp"
#include "include/IntersectionDetection2D.hpp"
#include "problemOperators/poissonGreenFunctionOperator2D.hpp"

//Simple block mesh data
//for uniform quad
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
  bcDirch.push_back(Polyline2D<REAL>({Vec2D<REAL>({0.2, 0.2})
                                         ,Vec2D<REAL>({0.6, 0.0})
                                         ,Vec2D<REAL>({1.0, 0.2})  }));
  bcDirch.push_back(Polyline2D<REAL>({Vec2D<REAL>({1.0, 1.0})
                                         ,Vec2D<REAL>({0.6, 0.8})
                                         ,Vec2D<REAL>({0.2, 1.0})  }));

  bcNeum.push_back(Polyline2D<REAL>({Vec2D<REAL>({1.0, 0.2})
                                      , Vec2D<REAL>({0.8, 0.6})
                                      , Vec2D<REAL>({1.0, 1.0})  }));
  bcNeum.push_back(Polyline2D<REAL>({Vec2D<REAL>({0.2, 1.0})
                                      , Vec2D<REAL>({0.0, 0.6})
                                      , Vec2D<REAL>({0.2, 0.2})} ));
};


int main(){
  //The MPI communicator
  bool IS_MPI_ON=false;
  MPIComm mpiComm(IS_MPI_ON);

  //Problem base size
  const int s = 10;//128;              //Image length-width
  const int nSize = s*s;          //Total image size
  const double dx= 1.0/double(s); //Image increment

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

  //MPI global partitioning
  MPI_ITstart = firstIterator<int>(mpiComm.getProcID(), mpiComm.getNProcs(), nSize);
  MPI_ITend   = lastIterator<int>( mpiComm.getProcID(), mpiComm.getNProcs(), nSize);
  MPI_lSize   = MPI_ITend - MPI_ITstart;

  //Generate a vector for storing
  //the local solution data
  std::vector<int> InDomainFlag; InDomainFlag.resize(MPI_lSize);
  std::vector<double> u_sol;     u_sol.resize(MPI_lSize);

  //device core local partitioning
  //all of the code from here should
  //be run on the device
  int dev_id=0, ndev_cores=1, It1D=0, I=0;
  double zero(0.00);

//  #pragma omp default(shared) private(dev_lSize, dev_ITstart, dev_ITend, dev_id, It1D, I)
 // #pragma omp target
  {
    dev_ITstart = firstIterator<int>(dev_id, ndev_cores, MPI_lSize);
    dev_ITend   = lastIterator<int>( dev_id, ndev_cores, MPI_lSize);
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
      int rnd_seed = I;
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


  //clean-up
  InDomainFlag.clear();
  u_sol.clear();
  MPI_Finalize();
  return 0;
}
