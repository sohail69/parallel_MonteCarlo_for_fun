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

/*****************************************\
! Definition of a ParaView VTK writer
! class for simple problems
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
                                        , wostr_part(wostr_part_)
{


};


/*****************************************\
! Class destructor
!
\*****************************************/
template<typename real, size_t sdim>
ParaViewWriter<real,sdim>::~ParaViewWriter()
{

};


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
    OutFile << "# vtk DataFile Version 2.0"       << std::endl;
    OutFile << "Block mesh"                       << std::endl;
    OutFile << "ASCII"                            << std::endl;
    OutFile << "DATASET UNSTRUCTURED_GRID"        << std::endl << std::endl;
    OutFile << "POINTS " + std::to_string(nnTot)+" float"<< std::endl;

    //Output the Grid points
    std::array<real,sdim> Point;
    for(int i=0; i<nnTot; i++){
      BHMeshPoint(Point.data(), i, BHMeshData);
      for(int j=0; j<sdim; j++) OutFile << std::setw(15) << Point[j];
      for(int j=sdim; j<3; j++) OutFile << std::setw(15) << 0.00;
      OutFile << std::endl;
    }

    //Output the cell headers
    std::string HName="CELLS  " + std::to_string(nCells) + "  "+ std::to_string(nVars*nnTot);
    OutFile << std::endl << HName << std::endl;

    //Output the cell topologies
    for(int iCell=0; iCell<nCells; iCell++){
      for(int i=0; i<std::pow(2,sdim); i++){
      //  for(int j=0; j<sdim; j++) OutFile << std::setw(15) << Point[j];
      //  for(int j=sdim; j<3; j++) OutFile << std::setw(15) << 0.00;
      //  OutFile << std::endl;
      }
    }
  }


  if(wostr_part.procID == 0) OutFile.close();
};

