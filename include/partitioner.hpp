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
FORCE_INLINE UINT minPartitionSize(const UINT nPartitions, const UINT nsize){
  return nsize/nPartitions;
};

//Remainder of uniform
//minimal partitioning
template<typename UINT>
FORCE_INLINE UINT partitionRemainder(const UINT nPartitions, const UINT nsize){
  return (nsize - nPartitions*minPartitionSize(nPartitions, nsize));
};

//Starting iterator of the
//processing unit
template<typename UINT>
FORCE_INLINE UINT firstIterator(const UINT tiD, const UINT nPartitions, const UINT nsize){
  UINT rem  = partitionRemainder(nPartitions, nsize);
  UINT Nmin = minPartitionSize(nPartitions, nsize);
  UINT tstart0 = tiD*(Nmin+1) + 1;
  UINT tstart1 = tiD*Nmin + rem + 1;
  UINT tstart2 = tiD*Nmin + 1;
  return ((rem!=0) ? ((tiD<rem)? tstart0 : tstart1) : tstart2) - 1;
};

//Last iterator of the
//processing unit
template<typename UINT>
FORCE_INLINE UINT lastIterator(const UINT tiD, const UINT nPartitions, const UINT nsize){
  UINT rem  = partitionRemainder(nPartitions, nsize);
  UINT Nmin = minPartitionSize(nPartitions, nsize);
  UINT tnend0 = (tiD+1)*(Nmin+1);
  UINT tnend1 = (tiD+1)*Nmin+ rem;
  UINT tnend2 = (tiD+1)*Nmin;
  return ((rem!=0) ? ((tiD<rem)? tnend0 : tnend1) : tnend2) - 1;
};

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
FORCE_INLINE UINT findIteratorOwner(const UINT GnodeID, const UINT nPartitions, const UINT nsize){
  UINT rem    = partitionRemainder(nPartitions, nsize);
  UINT Nmin   = minPartitionSize(nPartitions, nsize);
  UINT Nmax   = (rem == 0) ? Nmin : (Nmin + 1);
  UINT rNsize = rem*Nmax;
  UINT A      = GnodeID - rNsize;
  UINT pidA   = GnodeID/Nmax;
  UINT B      = A/Nmin;
  UINT pidB   = rem + B;
  return ( ((GnodeID <=  rNsize)and(rem != 0)) ? pidA : pidB);
};


/***************************\
! Partitioning a problem
! that uses accumulators 
! for reduce operations
! assuming an existing 
! geometric partitioning
\***************************/
//Find optimal~heuristic
//accumulator size
template<typename UINT>
FORCE_INLINE UINT findAccumPartitionSize(const UINT GnodeID, const UINT nPartitions, const UINT nAccum, const UINT nAccumSize){
//  if(nAccum < nPartitions)
  return UINT(0);
};



#endif







