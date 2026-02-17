#pragma once
#include "globalMacros.hpp"
#include "blockHyperMesh.hpp"

/***************************\
! Partitioner struct that
! has all data relating
! to a partitioned dataset
\***************************/
template<typename uint>
struct PACKSTRUCT WoStr_Partition
{
  // The problem total size
  uint total_samplePoints, nWalks0;

  // MPI Partitioning strictly
  // partitions sample points
  uint nProcs, procID;
  uint mpi_Istart, mpi_Iend, mpi_lsize;

  // Device\thread Partitioning
  uint nThreadsPerMPI=1;

  // Accumulator sizing
  uint ThreadPart_GCD=1;
  uint nAccumsPerMPI=0, nAccumsPerPart=1;
  uint nAccumsPerThread=1;

  // Random walk samples
  // partitioned
  uint nWalksPerAccum=0, nWalks=0;
};

/***************************\
! Linear partitioning
! routines for problems
! where the geometric size
! is bigger than the number
! of threads/MPI-procs 
! etc...
!
! Starts numbering from 0
\***************************/
//Starting iterator of the
//processing unit
template<typename uint>
FORCE_INLINE uint firstIterator(const uint tiD, const uint nPartitions, const uint nsize){
  uint rem  =  nsize%nPartitions;
  uint Nmin = nsize/nPartitions;
  uint tstart0 = tiD*(Nmin+1);
  uint tstart1 = tiD*Nmin + rem;
  uint tstart2 = tiD*Nmin;
  return ((rem!=0) ? ((tiD<rem)? tstart0 : tstart1) : tstart2);
};

//Last iterator of the
//processing unit
template<typename uint>
FORCE_INLINE uint lastIterator(const uint tiD, const uint nPartitions, const uint nsize){
  uint rem  =  nsize%nPartitions;
  uint Nmin = nsize/nPartitions;
  uint tnend0 = (tiD+1)*(Nmin+1);
  uint tnend1 = (tiD+1)*Nmin+ rem;
  uint tnend2 = (tiD+1)*Nmin;
  return ((rem!=0) ? ((tiD<rem)? tnend0 : tnend1) : tnend2) - 1;
};


//Finds the greatest common
//divisor for a pair of numbers
template<typename uint>
FORCE_INLINE uint GCD_pairAlg(const uint &A, const uint &B){
  uint a=A, b=B, t;
  if( (a>0) && (b>0) ){
    for(;B != 0;){
      t = b;
      b = a%b;
      a = t;
    }
  }
  return a;
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
template<typename uint>
FORCE_INLINE uint findIteratorOwner(const uint GnodeID, const uint nPartitions, const uint nsize){
  uint rem    =  nsize%nPartitions;
  uint Nmin   = nsize/nPartitions;
  uint Nmax   = (rem == 0) ? Nmin : (Nmin + 1);
  uint rNsize = rem*Nmax;
  uint A      = GnodeID - rNsize;
  uint pidA   = GnodeID/Nmax;
  uint B      = A/Nmin;
  uint pidB   = rem + B;
  return ( ((GnodeID <=  rNsize)and(rem != 0)) ? pidA : pidB);
};


/***************************\
! Partitioning a problem
! that uses accumulators 
! for reduce operations
! this is for cases where
! the cores exceed the
! problem size geometrically
! but not the number of
! samples at each point
\***************************/
//Find the number of accumulators
//sample-sets belonging to a
//given partition
template<typename uint>
FORCE_INLINE uint FindLocalNumberOfAccums(const uint &PartID, const uint &nDevThreads, const uint &nLocalParts){
  uint remThreads = nDevThreads%nLocalParts;
  uint minThreads = nDevThreads/nLocalParts;
  return ((remThreads==0)or(PartID >= remThreads))? minThreads : minThreads+1;
};

//Find the geometric parentID of a
//given accumlator sample-set
template<typename uint>
FORCE_INLINE uint FindAccumParentID(const uint &threadID, const uint &nDevThreads, const uint &nLocalParts){
  uint remThreads = nDevThreads%nLocalParts;
  uint minThreads = nDevThreads/nLocalParts;
  uint maxThreads = (remThreads == 0) ? minThreads : (minThreads+1);
  uint Dividend   = (threadID - remThreads*maxThreads);

  uint Guess1 = threadID/maxThreads;
  uint Guess2 = (Dividend/minThreads) + remThreads;
  uint parent_partition = (Dividend <= 0) ? Guess1 : Guess2;
  return parent_partition;
};

//Find the first position of an 
//accumulator belonging to a given
//parentID partition
template<typename uint>
FORCE_INLINE uint FindPartIDFirstLocalAccumPos(const uint &PartID
                                             , const uint &nDevThreads
                                             , const uint &nLocalParts)
{
  uint remThreads = nDevThreads%nLocalParts;
  uint minThreads = nDevThreads/nLocalParts;
  uint maxThreads = (remThreads == 0) ? minThreads : (minThreads+1);

  uint Guess1 = PartID*maxThreads;
  uint Guess2 = (PartID-remThreads)*minThreads + remThreads*maxThreads;
  uint AccumStartPos = ((remThreads != 0)or(PartID < remThreads)) ? Guess1 : Guess2;
  return AccumStartPos;
};



/***************************\
! Partition from the block
! mesh 
\***************************/
//For a device/thread level parallelism further partitioning
//is needed, however for this case there are 2 possibilities:
// 1: There is a device/multi-cores and are less partitions than threads (split partitions between threads)
// 2: There is a device/multi-cores and are more partitions than threads (split walks between threads)
// 3: There is no device, so just use original MPI-partitions
template<typename real, typename uint, size_t sdim>
FORCE_INLINE void partitionProblem(MPIComm & comm
                                 , const blockHyperMeshData<real,uint,sdim> & BHMeshData
                                 , const uint & nWalks
                                 , WoStr_Partition<uint> & part)
{
  //Get the total problem sample
  //point size from mesh
  part.total_samplePoints=1;
  for(uint i=0; i<sdim; i++) part.total_samplePoints *= BHMeshData.sizes[i];
  part.nWalks0 = nWalks;

  //Get the MPI process data
  //and from the communicator
  part.nProcs = comm.getNProcs();
  part.procID = comm.getProcID();

  //Partition the problem (Geometrically)
  //using the MPI-partitioner
  part.mpi_Istart = firstIterator<uint>(part.procID, part.nProcs, part.total_samplePoints);
  part.mpi_Iend   = lastIterator<uint>( part.procID, part.nProcs, part.total_samplePoints);
  part.mpi_lsize  = part.mpi_Iend - part.mpi_Istart + 1;

  //Partition the full local problem
  //among the threads\devices this is
  //for both the geometry and the walks
  part.ThreadPart_GCD   = GCD_pairAlg( part.mpi_lsize, part.nThreadsPerMPI);
  part.nAccumsPerMPI    = part.mpi_lsize * (part.nThreadsPerMPI/part.ThreadPart_GCD);
  part.nAccumsPerPart   = (part.nThreadsPerMPI/part.ThreadPart_GCD);
  part.nAccumsPerThread = (part.mpi_lsize/part.ThreadPart_GCD);

  //Partition the walks among the
  //accumulators and aim for approximately
  //the initial number of walks
  part.nWalksPerAccum = part.nWalks0/part.nAccumsPerPart; 
  part.nWalks = part.nWalksPerAccum * part.nAccumsPerPart;
};
