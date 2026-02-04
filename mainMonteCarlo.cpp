// (Slow) implementation of the "Walk on Stars" algorithm for Laplace equations.
// Corresponds to the estimator given in Equation 18 of Sawhney et al,
// "Walk on Stars: A Grid-Free Monte Carlo Method for PDEs with Neumann Boundary
// Conditions" (2023), assuming no source term and zero-Neumann conditions.
// NOTE: this code makes a few shortcuts for the sake of code brevity; may
// be more suitable for tutorials than for production code/evaluation.
// To compile: c++ -std=c++17 -O3 -pedantic -Wall WoStLaplace2D.cpp -o wost
#include <algorithm>
#include <array>
#include <complex>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>

#include "include/MPIcomm.hpp"
#include "include/partitioner.hpp"
#include "include/blockHyperMesh.hpp"
#include "include/localVectorAlgebra.hpp"
using namespace std;

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


int main(){
  //The MPI communicator
  bool IS_MPI_ON=false;
  MPIComm mpiComm(IS_MPI_ON);

  //Problem base size
  const int s = 128;       //Image length-width
  const int nSize = s*s;   //Total image size
  const double dx=0.01;

  //Generate the blockMesh from
  //the base size (2D-Quad) 
  blockHyperMeshData<double,int,2> BHMeshData;
  simpleBlockMeshBuild<double,int,2>(dx, s, BHMeshData);

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
  std::vector<int> InDomainFlag(MPI_lSize);
  std::vector<double> u_sol(MPI_lSize);

  //device core local partitioning
  //all of the code from here should
  //be run on the device
  #pragma omp target
  {
    int dev_id=1, ndev_cores=10, I=0;
    dev_ITstart = firstIterator<int>(dev_id, ndev_cores, nSize);
    dev_ITend   = lastIterator<int>( dev_id, ndev_cores, nSize);
    dev_lSize   = dev_ITend - dev_ITstart;

    for(Iter1D=dev_ITstart; I<dev_ITend; I++){
      Vec2D<double> x0;
      BHMeshPoint<double,int,2>(x0.data(), Iter1D, BHMeshData);
      InDomainFlag[I] = (insideDomain(x0, boundaryDirichlet, boundaryNeumann) ? 1:0);
    }
  }

  //Printing some extra
  //data for fun
  if(mpiComm.getProcID() == 0) std::cout << "Iterators" << std::endl;
  std::cout << std::setw(10) << mpiComm.getProcID()
            << std::setw(10) << mpiComm.getNProcs() 
            << std::setw(10) << MPI_lSize
            << std::setw(10) << MPI_ITstart
            << std::setw(10) << MPI_ITend   << std::endl;
  return 0;
}

/*
int main( int argc, char** argv ) {
   srand( time(NULL) );
   ofstream out( "out.csv" );

   int s = 128; // image size
   for( int j = 0; j < s; j++ )
   {
      cerr << "row " << j << " of " << s << endl;
      for( int i = 0; i < s; i++ )
      {
         Vec2D x0( ((double)i+.5)/((double)s),
                   ((double)j+.5)/((double)s) );
         double u = 0.;
         if( insideDomain(x0, boundaryDirichlet, boundaryNeumann) )
            u = solve( x0, boundaryDirichlet, boundaryNeumann, lines );
         out << u;
         if( i < s-1 ) out << ",";
      }
      out << endl;
   }
   return 0;
}*/

