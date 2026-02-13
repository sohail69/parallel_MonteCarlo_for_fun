/**************************************\
! An implementation of the bounding
! volume hierarchy storage for storing 
! geometric entities (specifically the
! boundary but could be other things)
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
#include "../templatedGeometry/geometricPrimitives.hpp"
#include "../templatedGeometry/localVectorAlgebra.hpp"

// A single bounded volume
typename<typename real, typename uint, size_t sdim>
struct BoundedVolume
{
  //Identification
  uint myID, parentBV_id;

  //Owned entities
  uint nOwnedEntities;  //If not a leaf (then give idea as to how many entities owne)
  uint left_childBV_id, rightchildBV_id;

  //If a leaf give details about position of the entity
  uint entityID=-1, entityPos=-1;

  //Bounding box geomtric size
  //and orientation
  Point<real,sdim> startPoint, endpoint;
  VecND<real,sdim> boxAxis;
};

//The bounded volume hierarchy 
template<typename real, typename uint, size_t sdim>
using BVHstore = std::vector<BoundedVolume<real,uint,sdim>>;

//Generate the size of the tree (lol binary trees
//are double the size of the initial data)
template<typename uint>
uint sizeOfBVH(const uint & nEntities){ return 2*nEntities;};


