#pragma once
#include "globalMacros.hpp"

//Miminimum data needed
//to generate and oriented
//block hyper mesh.
template<typename real, typename uint, uint dim>
struct PACKSTRUCT blockHyperMeshData{
   uint sizes[dim];
   real offset[dim], dx[dim*dim];
};

//The hyper-block mesh
//forward iterator
template<typename real, typename uint, uint dim>
FORCE_INLINE uint BHMeshFwdIter(const uint Iters[dim]
                              , const blockHyperMeshData<real,uint,dim> & BHMdata)
{
  uint Pos(0);
  #pragma unroll
  for(uint I=0; I<dim; I++){
    uint pos_b=Iters[I];
    #pragma unroll
    for(uint J=I+1; J<dim; J++){
      uint K=dim-J-1;
      pos_b *= BHMdata.sizes[K];
    }
    Pos += pos_b;
  }
  return Pos;
}

//The hyper-block mesh
//Inverse iterator
template<typename real, typename uint, uint dim>
FORCE_INLINE void BHMeshInvIter(uint ItersND[dim]
                              , const uint & Iter1D
                              , const blockHyperMeshData<real,uint,dim> & BHMdata)
{
  uint a = Iter1D;
  for(uint I=0; I<dim; I++){
    uint temp = uint(a/BHMdata.sizes[dim-I-1]);
    ItersND[dim-I-1] = a - temp*BHMdata.sizes[dim-I-1];
    a = temp;
  }
};


//Gives a physical point
//from a hyper-block mesh
//iterator
template<typename real, typename uint, uint dim>
FORCE_INLINE void BHMeshPoint(real Point[dim], const uint & Iter1D, const blockHyperMeshData<real,uint,dim> & BHMdata){
  uint ItersND[dim];
  BHMeshInvIter(ItersND, Iter1D, BHMdata);
  for(uint I=0; I<dim; I++){
    Point[I] =  BHMdata.offset[I];
    for(uint J=0; J<dim; J++){
      Point[I] += real(double(ItersND[I]) + 0.50)*BHMdata.dx[I*dim + J];
    }
  }
};
