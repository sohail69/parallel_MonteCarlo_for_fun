#pragma once

//Miminimum data needed
//to generate and oriented
//block hyper mesh.
template<typename REAL, typename UINT, UINT dim>
struct blockHyperMeshData{
   UINT sizes[dim];
   REAL offset[dim], dx[dim*dim];
};

//The hyper-block mesh
//forward iterator
template<typename REAL, typename UINT, UINT dim>
UINT BHMeshFwdIter(const UINT Iters[dim], const blockHyperMeshData<REAL,UINT,dim> & BHMdata){
  UINT Pos(0);
  #pragma unroll
  for(UINT I=0; I<dim; I++){
    UINT pos_b=Iters[I];
    #pragma unroll
    for(UINT J=I+1; J<dim; J++){
      UINT K=dim-J-1;
      pos_b *= BHMdata.sizes[K];
    }
    Pos += pos_b;
  }
  return Pos;
}

//The hyper-block mesh
//Inverse iterator
template<typename REAL, typename UINT, UINT dim>
void BHMeshInvIter(UINT ItersND[dim], const UINT & Iter1D, const blockHyperMeshData<REAL,UINT,dim> & BHMdata){
  UINT a = Iter1D;
  for(UINT I=0; I<dim; I++){
    UINT temp = UINT(a/BHMdata.sizes[dim-I-1]);
    ItersND[I] = a - temp*BHMdata.sizes[dim-I-1];
    a = temp;
  }
};


//Gives a physical point
//from a hyper-block mesh
//iterator
template<typename REAL, typename UINT, UINT dim>
void BHMeshPoint(REAL Point[dim], const UINT & Iter1D, const blockHyperMeshData<REAL,UINT,dim> & BHMdata){
  UINT ItersND[dim];
  BHMeshInvIter(ItersND, Iter1D, BHMdata);
  for(UINT I=0; I<dim; I++){
    Point[I] =  BHMdata.offset[I];
    for(UINT J=0; J<dim; J++){
      Point[I] += REAL(double(ItersND[I]))*BHMdata.dx[I*dim + J];
    }
  }
};
