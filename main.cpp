//Parallel libraries
#include <omp.h>
#include <mpi.h>

//Standard template libraries
#include <cmath>
#include <iomanip>
#include <iostream>
#include <functional>

//My source libraries
#include "include/globalMacros.hpp"
#include "include/partitioner.hpp"
#include "include/hierachicalTensorFE.hpp"

//My problem libraries
#include "problemIntegrators/solid_mixedUP_element.hpp"


using namespace std;

// mpic++ -fopenmp  -o main main.cpp
//


int main(){
  const unsigned nDevs=2, nthreads=2;
  const unsigned N=20;
  unsigned x_vec[N], y_vec[N], b_vec[N];

  omp_set_num_threads(nthreads);
  unsigned num_devices = omp_get_num_devices();
  unsigned I_start, I_end, I, Idev;
  unsigned nteams, tid;
  unsigned devIstart[nDevs], devIend[nDevs];

  const unsigned nsampl1D=3, pOrder=2, DIM=3;
  const unsigned nsamplND=pow(nsampl1D,DIM);
  const unsigned NDofs=pow((pOrder+1),DIM);

  //Set up the finite elements
  GaussLegendreInt<unsigned,double,nsampl1D,1>   GLInteg1D;
  GaussLegendreInt<unsigned,double,nsamplND,DIM> GLIntegND;
  CalcGaussLegendreXW1D<unsigned,double,nsampl1D>(GLInteg1D);


  for(int I=0; I<nsamplND; I++){
    cout << setw(20) << GLIntegND.Xi[0][I]
         << setw(20) << GLIntegND.Xi[1][I]
         << setw(20) << GLIntegND.Xi[2][I]
         << setw(20) << GLIntegND.Wi[I]    << endl;
  }

//////////  #pragma omp target private(I_start, I_end, I)
  #pragma omp parallel default(shared) private(tid, Idev, I, I_start, I_end, devIstart,  devIend)
  {
//    nthreads = omp_get_num_threads();
    nteams = omp_get_num_teams();
    tid = omp_get_thread_num();
//  I_start = firstIterator<unsigned>(tid, nthreads, N);
//  I_end   = lastIterator<unsigned>(tid, nthreads, N);
    NestedPartitioner<unsigned,nDevs,nthreads>(N, tid, I_start, I_end, devIstart,  devIend);

    #pragma omp barrier
  }
  return 0;
};














