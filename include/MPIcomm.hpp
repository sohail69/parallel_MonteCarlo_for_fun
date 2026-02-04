#pragma once
#include <mpi.h>

class MPIComm{
  private:
  int procID, nProcs;
  int MPIerr;
  MPI_Comm comm=MPI_COMM_WORLD;

  public:
    //The constructor
    MPIComm(bool IS_MPI_ON);

    //The destructor
    ~MPIComm();

    //Get MPI process data
    int & getProcID(){ return procID;};
    int & getNProcs(){ return nProcs;};
};

//The constructor
MPIComm::MPIComm(bool IS_MPI_ON)
{
  //Initialise MPI
  if(IS_MPI_ON != true) MPI_Init(NULL, NULL);
  MPI_Comm_size(comm, &nProcs);
  MPI_Comm_rank(comm, &procID);
};

//The destructor
MPIComm::~MPIComm()
{
  //Finalize MPI
  MPI_Finalize();
};
