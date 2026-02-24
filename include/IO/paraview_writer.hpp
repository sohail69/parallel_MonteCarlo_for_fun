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
    for(int j=0; j<sdim; j++) Offsets[i*sdim + j] = (i/intPow<int,int>(2,j))%2;
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
    //Referenced objects
    const blockHyperMeshData<real,int,sdim> & BHMeshData;
    const WoStr_Partition<int> & wostr_part;
    MPIComm & mpiComm;

    //Write the data of the VTK file
    void DataWrite(std::ofstream& oFile
                 , const std::string & datName
                 , const std::vector<real> & data
                 , const size_t nVars);

  public:
    //Constructor
    ParaViewWriter(const blockHyperMeshData<real,int,sdim> & BHMeshData_
                 , const WoStr_Partition<int> & wostr_part_
                 , MPIComm & mpiComm_);

    //Destructor
    ~ParaViewWriter();

    //Write data file
    void VTKwrite(const std::string & fName
                 , const std::string & datName
                 , const std::vector<real> & data
                 , const size_t nVars);
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
! The simple VTK writer
!
\*****************************************/
template<typename real, size_t sdim>
void ParaViewWriter<real,sdim>::VTKwrite(const std::string & fName
                                       , const std::string & datName
                                       , const std::vector<real> & data
                                       , const size_t nVars)
{
  //Data objects that only exist for
  //the duration of the write
  std::ofstream OutFile;  //The output 
  int nnTot=1, nCells=1;
  for(int i=0; i<sdim; i++) nnTot  *= BHMeshData.sizes[i];
  for(int i=0; i<sdim; i++) nCells *= (BHMeshData.sizes[i] -1);

  //If process 0, then open file and set
  //the headers output points and elms
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
      for(int j=0; j<3; j++) OutFile << std::setw(15) << ((j < sdim) ? Point[j] : 0.00);
      OutFile << std::endl;
    }

    //Calculate the necessary topological
    //information for n-Dimensional hyperCube
    const int nNodes = intPow<int,size_t>(2, sdim);
    int ItersND[sdim], ItersND_base[sdim],  Csizes[sdim], Nodes[nNodes], Offsets[nNodes*sdim];
    for(int i=0; i<sdim; i++) Csizes[i] = BHMeshData.sizes[i] - 1;
    hyperCubeOffsets<sdim,nNodes>(Offsets);

    //Output the cell headers
    std::string HName="CELLS  " + std::to_string(nCells) + "  " + std::to_string((nNodes+1)*nCells);
    OutFile << std::endl << HName << std::endl;

    //Output the cell topologies
    //for n-Dimensional hyperCube
    for(int iCell=0; iCell<nCells; iCell++){
      BHMeshInvIter<int,sdim>(ItersND_base, iCell, Csizes);
      for(int i=0; i<nNodes; i++){
        for(int k=0; k<sdim; k++) ItersND[k] = ItersND_base[k] + Offsets[i*sdim + k];
        Nodes[i] = BHMeshFwdIter<real,int,sdim>(ItersND, BHMeshData);
      }
      OutFile << std::setw(10) << nNodes;
      for(int i=0; i<nNodes; i++) OutFile << std::setw(10) << Nodes[i];
      OutFile << std::endl;
    }

    //Output the Cell types
    OutFile << "CELL_TYPES " + std::to_string(nCells) << std::endl;
    int CTypes[3] = {3, 8, 11};
    int CType = ((sdim <= 3)&&(sdim >= 1)) ? CTypes[sdim-1] : 7;
    for(int iCell=0; iCell<nCells; iCell++) OutFile << std::setw(10) << CType << std::endl;
  }

  //Write the data
  DataWrite(OutFile, datName, data, nVars);

  //Close the file
  if(wostr_part.procID == 0) OutFile.close();
};

/*****************************************\
! Output some data (used by the VTK
! writer)
\*****************************************/
template<typename real, size_t sdim>
void ParaViewWriter<real,sdim>::DataWrite(std::ofstream& oFile
                                        , const std::string & datName
                                        , const std::vector<real> & send_buf
                                        , const size_t nVars)
{
  //If process 0 write
  //out the headers
  int procID = wostr_part.procID;
  int nCellsTot = wostr_part.total_samplePoints;
  if(procID == 0){
    oFile << "POINT_DATA "+std::to_string(nCellsTot)<< std::endl;
    oFile << "SCALARS "+datName+" float "          << std::endl;
    oFile << " LOOKUP_TABLE "+datName+"Tab"        << std::endl;
  }

  //The MPI task requests
  MPI_Request req;
  MPI_Status statu;
  int MPIerr;
  MPI_Datatype MPIDtype=MPI_DOUBLE;

  //Send-recieve and Output the VTK Mesh data
  int part_size = nVars*PartitionSize<int>(procID, wostr_part.nProcs, wostr_part.total_samplePoints);
  if(procID != 0){
    MPIerr=MPI_Isend(send_buf.data(),part_size,MPIDtype,0,procID,mpiComm.getComm(),&req);
    std::cout << "PROC_ID : " << procID << " " << part_size << std::endl;
  }
  MPI_Barrier(mpiComm.getComm());

  std::vector<real> recv_buf(part_size);
  if(procID == 0){
     for(int iProc=0; iProc<wostr_part.nProcs; iProc++){
       int nCells_part = PartitionSize<int>(iProc, wostr_part.nProcs, wostr_part.total_samplePoints);
       int nVarsTot = nCells_part*nVars;
       recv_buf.resize(nVarsTot);
       std::cout << "I_PROC_ID : " << procID << " " << nCells_part << std::endl;

       if(iProc!=0) MPIerr=MPI_Recv(recv_buf.data(),nVarsTot,MPIDtype,iProc,iProc,mpiComm.getComm(),&statu);
       for(int iCell=0; iCell<nCells_part; iCell++){
         for(int iVars=0; iVars<nVars; iVars++){
           oFile << std::setw(12) << ((iProc==0)? send_buf[iCell*nVars+iVars] : recv_buf[iCell*nVars+iVars]);
         }
         oFile << std::endl;
       }
     }
  }
  MPI_Barrier(mpiComm.getComm());
  recv_buf.clear();
};
