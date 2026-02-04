#ifndef PARTITIONER_HPP
#define PARTITIONER_HPP

#include "globalMacros.hpp"

/***************************\
! Linear partitioning
! routines
!
! Starts numbering from 0
\***************************/
//Minimum partition size
template<typename UINT>
FORCE_INLINE UINT minPartitionSize(const UINT nthreads, const UINT nsize){
  return nsize/nthreads;
};

//Remainder of uniform
//minimal partitioning
template<typename UINT>
FORCE_INLINE UINT partitionRemainder(const UINT nthreads, const UINT nsize){
  return (nsize - nthreads*minPartitionSize(nthreads, nsize));
};

//Starting iterator of the
//processing unit
template<typename UINT>
FORCE_INLINE UINT firstIterator(const UINT tiD, const UINT nthreads, const UINT nsize){
  UINT rem  = partitionRemainder(nthreads, nsize);
  UINT Nmin = minPartitionSize(nthreads, nsize);
  UINT tstart0 = tiD*(Nmin+1) + 1;
  UINT tstart1 = tiD*Nmin + rem + 1;
  UINT tstart2 = tiD*Nmin + 1;
  return ((rem!=0) ? ((tiD<rem)? tstart0 : tstart1) : tstart2) - 1;
};

//Last iterator of the
//processing unit
template<typename UINT>
FORCE_INLINE UINT lastIterator(const UINT tiD, const UINT nthreads, const UINT nsize){
  UINT rem  = partitionRemainder(nthreads, nsize);
  UINT Nmin = minPartitionSize(nthreads, nsize);
  UINT tnend0 = (tiD+1)*(Nmin+1);
  UINT tnend1 = (tiD+1)*Nmin+ rem;
  UINT tnend2 = (tiD+1)*Nmin;
  return ((rem!=0) ? ((tiD<rem)? tnend0 : tnend1) : tnend2) - 1;
};


/***************************\
! Nested Linear partitioning
! routine
! partitioning across
! compute-MPI-nodes and
! device-level
\***************************/
template<typename UINT, UINT nDevs, UINT nNodes>
FORCE_INLINE void NestedPartitioner(const UINT & nsize, UINT & mpiID, UINT & Istart, UINT & Iend
                                  , UINT devIstart[nDevs], UINT devIend[nDevs])
{
  //MPI Proc partitioning
  Istart = firstIterator(mpiID, nNodes, nsize);
  Iend = lastIterator(mpiID, nNodes, nsize);

  //Device block partitioning
  UINT nodeSize = (Iend - Istart) + 1;
  for(UINT Idev=0; Idev<nDevs; Idev++){
    devIstart[Idev] = Istart + firstIterator(Idev, nDevs, nodeSize);
    devIend[Idev]   = Istart + lastIterator(Idev, nDevs, nodeSize);
  }
}

/***************************\
! Forward and inverse 
! Iterators
!
! used for finding the owners
! of the partitioned entities
! whether they be compute-
! nodes or devices
\***************************/
//Finds the owner process
//or compute node of a
//given iterator in a
//partitioned problem
template<typename UINT>
FORCE_INLINE UINT findIteratorOwner(const UINT GnodeID, const UINT nthreads, const UINT nsize){
  UINT rem    = partitionRemainder(nthreads, nsize);
  UINT Nmin   = minPartitionSize(nthreads, nsize);
  UINT Nmax   = (rem == 0) ? Nmin : (Nmin + 1);
  UINT rNsize = rem*Nmax;
  UINT A      = GnodeID - rNsize;
  UINT pidA   = GnodeID/Nmax;
  UINT B      = A/Nmin;
  UINT pidB   = rem + B;
  return ( ((GnodeID <=  rNsize)and(rem != 0)) ? pidA : pidB);
};

/*

! Calculates the local 
! numberingID from global
! numberingID of a particular
! node/element
*/

/*
template<typename UINT>
FORCE_INLINE UINT globalToLocalIter(const UINT GnodeID, const UINT tiD, const UINT nthreads, const UINT nsize){
  UINT rem    = partitionRemainder(nthreads, nsize);
  UINT Nmin   = minPartitionSize(nthreads, nsize);
  UINT Nmax   = (rem == 0) ? Nmin : (Nmin + 1);
  UINT rNsize = rem*Nmax;

  UINT rem = partitionRemainder(nthreads, nsize);
  UINT A   = GnodeID - rem*(minPartitionSize(nthreads, nsize)+1);

  return ((rem!=0) ? ((tiD<rem)? tnend0 : tnend1) : tnend2) - 1;
};


PURE FUNCTION GLOBAL_TO_LOCAL1(G_nodeID,nn_pp1,nn_pp2,num,nprocs) RESULT(L_nodeID)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: G_nodeID, num, nn_pp1, nn_pp2, nprocs;
  INTEGER            :: L_nodeID, procID, threshhold;

  threshhold = num*nn_pp1;
  procID = FIND_NODE_PROC1(G_nodeID,nn_pp1,nn_pp2,num,nprocs)
  procID = procID - 1;
  L_nodeID = G_nodeID - procID*nn_pp1
  IF(num /= 0)THEN
    IF(G_nodeID > threshhold)THEN
      L_nodeID = G_nodeID - num*nn_pp1 - (procID - num)*nn_pp2
    ENDIF
  ENDIF
ENDFUNCTION GLOBAL_TO_LOCAL1
*/
#endif
