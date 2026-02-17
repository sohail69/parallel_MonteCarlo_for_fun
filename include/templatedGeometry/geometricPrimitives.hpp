#pragma once
#include <array>
#include "../globalMacros.hpp"

//Facades used for making
//the code more clear


//A n-dimensional vector
//for storing variables
//e.g. Vars (operators for vectors)
template<typename real, size_t dim>
using VecND = std::array<real,dim>;

//A point
template<typename real, size_t sdim>
using Point = VecND<real,sdim>;

//vector initialisation from scalar
/*
template<typename real, size_t dim>
FORCE_INLINE void operator=(VecND<real,dim> & x, const real & a)
{
  for(size_t i=0; i<dim; i++) x[i] = a;
}*/

//vector scalar increment
template<typename real, size_t dim>
FORCE_INLINE void operator+=(VecND<real,dim> & x, const real & a)
{
  for(size_t i=0; i<dim; i++) x[i] += a;
}

//vector scalar decrement
template<typename real, size_t dim>
FORCE_INLINE void operator-=(VecND<real,dim> & x, const real & a)
{
  for(size_t i=0; i<dim; i++) x[i] /= a;
}

//vector scaling
template<typename real, size_t dim>
FORCE_INLINE void operator*=(VecND<real,dim> & x, const real & a)
{
  for(size_t i=0; i<dim; i++) x[i] *= a;
}

//vector inverse scaling
template<typename real, size_t dim>
FORCE_INLINE void operator/=(VecND<real,dim> & x, const real & a)
{
  for(size_t i=0; i<dim; i++) x[i] /= a;
}
