

/**************************************\
! Tensor element N-D
!
! Takes a 1-D element basis and forms
! higher order basis functions, mappings
! supported for Tensor elements e.g.
! lines, quads, hexs, etc...
! as arguments it takes Legendre
! 1-D polynomial modal shape functions
! and derivative
!
! Author: Sohail Rathore
! Date  : 31/01/2025
!
\**************************************/
namespace TensorElementND{
  //
  //  Form a higher dimension shape
  //  function weight using a 1-D shape
  //  function weight and the DOF position
  //
  template<typename UINT, typename REAL, UINT nSampl, UINT DIM>
  FORCE_INLINE REAL Ni_ND(const std::vector<std::array<REAL,nSampl>> Ni1D
                        , const std::array<UINT,DIM> IJK       //The Node
                        , const std::array<UINT,DIM> IJKGauss) //The Gauss-point
  {
    //Calculate shape function
    REAL Ni = Ni1D[IJK[0]][IJKGauss[0]];
    #pragma unroll
    for(UINT I=1; I<DIM; I++) Ni = Ni*Ni1D[IJK[I]][IJKGauss[I]];
    return Ni;
  };


  //
  //  Form a higher dimension shape function
  //  derivative weight using the 1-D shape
  //  function weight, derivative weight,
  //  DOF position and the Gauss-Point
  //
  template<typename UINT, typename REAL, UINT nSampl, UINT DIM>
  FORCE_INLINE
  std::array<REAL,DIM> dNi_ND(const std::vector<std::array<REAL,nSampl>> Ni1D
                               , const std::vector<std::array<REAL,nSampl>> dNi1D
                               , const std::array<UINT,DIM> IJK       //The Node
                               , const std::array<UINT,DIM> IJKGauss) //The Gauss-point
  {
    REAL Zero = REAL(0.0);
    REAL One  = REAL(1.0);
    std::array<REAL,DIM> dNi, mask;

    //Calculate shape function derivative
    #pragma unroll
    for(UINT I=0; I<DIM; I++){
      #pragma unroll
      for(UINT J=0; J<DIM; J++) mask[J]=Zero;
      mask[I] = One;

      dNi[I] = dNi[I]*( (One - mask[0])*Ni1D[IJK[0]][IJKGauss[0]] 
                              + mask[0]*dNi1D[IJK[0]][IJKGauss[0]] );
      #pragma unroll
      for(UINT J=1; J<DIM; J++){
        dNi[I] = dNi[I]*( (One - mask[J])*Ni1D[IJK[J]][IJKGauss[J]] 
                               + mask[J]*dNi1D[IJK[J]][IJKGauss[J]] );
      }
    }
    return dNi;
  };

  //
  //  For effiency purposes one may want to
  //  store the entire shape-function and shape
  //  function derivative weights array for 
  //  an N-D element so it can be accessed as
  //  oppose to being recalculated on every element
  //  in these cases this function can be used
  //
  template<typename UINT, typename REAL, UINT nSampl, UINT DIM>
  FORCE_INLINE
  void CalcND_SF_SFD(std::vector<std::array<REAL,nSampl>>                 NiND
                   , std::vector<std::array<std::array<REAL,DIM>,nSampl>> dNiND
                   , const std::vector<std::array<REAL,nSampl>>           Ni1D
                   , const std::vector<std::array<REAL,nSampl>>           dNi1D)
  {
  //TODO
  /*for(UINT I=0; I<; I++){
      NiND
      dNiND
    }*/
  };


  //
  // Assume the element shape Jacobian matrix
  // dni/dxj is only a function of the vertex
  // contributions this kind of means that the
  // that the other shape functions are just
  // affine transformations of their modal
  // coord-Shape functions
  //
  template<typename UINT, typename REAL, UINT DIM, UINT nSampl>
  FORCE_INLINE
  void calcElmJacDet_JacInv(Eigen::MatrixXd *JacInv
                          , REAL         *JacDet
                          , const UINT Isample
                          , const UINT nVerts
                          , const std::vector<std::array<std::array<REAL,DIM>,nSampl>>  dNiND
                          , const std::vector<std::array<REAL,DIM>>  coord){
    Eigen::MatrixXd Jac = Eigen::MatrixXd::Zero(DIM, DIM);
    for(int I=0; I<DIM; I++){
      for(int J=0; J<DIM; J++){
        for(int K=0; K<nVerts; K++){
          Jac(I,J) += dNiND[I][Isample][K]*coord[J][K];
        }
      }
    }
    *JacDet = Jac.determinant();
    *JacInv = Jac.inverse();
  }
};//End of TensorElementND namespace


