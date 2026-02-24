#pragma once

//Streams
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

//containers
#include <vector>
#include <string>

//Operations
#include <cmath>
#include <mpi.h>

//My Libs
#include "../globalMacros.hpp"
#include "../blockHyperMesh.hpp"
#include "../MPIcomm.hpp"
#include "../partitioner.hpp"
#include <mpi.h>

// Calculates a^b
// for integer numbers
template<typename uint1, typename uint2>
constexpr uint1 intPow(const uint1 & a, const uint2 & b){
  uint1 k=1;
  for(uint i=0; i<b; i++) k *= a;
  return k;
}

// Finds the mapping for a
// hyperCube element needed
// for topology
template<int sdim, int nNodes>
void hyperCubeOffsets(int Offsets[nNodes*sdim]){
  for(int i=0; i<nNodes; i++){
    for(int j=0; j<sdim; j++){
      Offsets[i*sdim + j] = (i/intPow<int,int>(2,j))%2;
      std::cout << std::setw(8) << Offsets[i*sdim + j];
    }
    std::cout << std::endl;
  }
}

/*****************************************\
! Definition of a ParaView VTK writer
! class for simple problems on the 
! hyperBlockMesh/Grid
\*****************************************/
template<typename real, size_t sdim>
class ParaViewWriter
{
  private:
    MPI_Request req;
    MPI_Status status;
    int MPIerr;

    //Referenced objects
    const blockHyperMeshData<real,int,sdim> & BHMeshData;
    const WoStr_Partition<int> & wostr_part;
    MPIComm & mpiComm;
  public:
    //Constructor
    ParaViewWriter(const blockHyperMeshData<real,int,sdim> & BHMeshData_
                 , const WoStr_Partition<int> & wostr_part_
                 , MPIComm & mpiComm_);

    //Destructor
    ~ParaViewWriter();

    //Write data
    void VTKwrite(const std::string & fName, const std::vector<real> & data, const size_t nVars);
};

/*****************************************\
! Class constructor
!
\*****************************************/
template<typename real, size_t sdim>
ParaViewWriter<real,sdim>::ParaViewWriter(const blockHyperMeshData<real,int,sdim> & BHMeshData_
                                        , const WoStr_Partition<int> & wostr_part_
                                        , MPIComm & mpiComm_) : BHMeshData(BHMeshData_), mpiComm(mpiComm_)
                                        , wostr_part(wostr_part_){};

/*****************************************\
! Class destructor
!
\*****************************************/
template<typename real, size_t sdim>
ParaViewWriter<real,sdim>::~ParaViewWriter(){};

/*****************************************\
! Output some data
!
\*****************************************/
template<typename real, size_t sdim>
void ParaViewWriter<real,sdim>::VTKwrite(const std::string & fName
                                       , const std::vector<real> & data
                                       , const size_t nVars)
{
  //Data objects that only exist for
  //the duration of the write
  std::ofstream OutFile;  //The output 
  int nnTot=1, nCells=1;
  for(int i=0; i<sdim; i++) nnTot  *= BHMeshData.sizes[i];
  for(int i=0; i<sdim; i++) nCells *= (BHMeshData.sizes[i] -1);

  //If process 0, then open
  //file and set the headers
  //output points and elms
  if(wostr_part.procID == 0){
    OutFile.open(fName + ".vtk");
    OutFile << "# vtk DataFile Version 2.0"      << std::endl;
    OutFile << "Block mesh"                      << std::endl;
    OutFile << "ASCII"                           << std::endl;
    OutFile << "DATASET UNSTRUCTURED_GRID"       << std::endl << std::endl;
    OutFile << "POINTS " + std::to_string(nnTot) + " float"   << std::endl;

    //Output the Grid points
    std::array<real,sdim> Point;
    for(int i=0; i<nnTot; i++){
      BHMeshPoint<real,int,sdim>(Point.data(), i, BHMeshData);
      for(int j=0; j<sdim; j++) OutFile << std::setw(15) << Point[j];
      for(int j=sdim; j<3; j++) OutFile << std::setw(15) << 0.00;
      OutFile << std::endl;
    }

    //Output the cell headers
    std::string HName="CELLS  " + std::to_string(nCells) + "  " + std::to_string(nVars*nCells);
    OutFile << std::endl << HName << std::endl;

    //Output the cell topologies
    //for n-Dimensional hyperCube
    const int nNodes = intPow<int,size_t>(2, sdim);
    int ItersND[sdim], ItersND_base[sdim], Nodes[nNodes], Offsets[nNodes*sdim];
    hyperCubeOffsets<sdim,nNodes>(Offsets);
    for(int iCell=0; iCell<nCells; iCell++){
      BHMeshInvIter<real,int,sdim>(ItersND_base, iCell, BHMeshData);
      for(int i=0; i<nNodes; i++){
        for(int k=0; k<sdim; k++) ItersND[k] = ItersND_base[k] + Offsets[i*sdim + k];
        Nodes[i] = BHMeshFwdIter<real,int,sdim>(ItersND, BHMeshData);
      }
      for(int i=0; i<nNodes; i++) OutFile << std::setw(10) << Nodes[i];
      OutFile << std::endl;
    }
  }


  if(wostr_part.procID == 0) OutFile.close();
};

