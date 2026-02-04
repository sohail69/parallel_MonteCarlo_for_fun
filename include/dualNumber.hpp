#pragma once
#include "../UtilityObjects/macros.hpp"

/***************************************\
!
!  Dual-Numbers
!
\***************************************/
template<typename value_t, typename gradient_t>
struct PACKSTRUCT dualNumber{
  //Definition of dual number
  value_t     val;
  gradient_t  grad;

  FORCE_INLINE dualNumber(value_t r=0.0, gradient_t eps=0.0): val(r), grad(eps){};
};

/***************************************\
!
!  Dual-Dual number
!  Equivalence/augmentation operations
!
\***************************************/
//Increment operator
template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator+=(dualNumber<v_t,g_t> & a, const dualNumber<v_t,g_t> & b) 
{
  a.val = a.val + b.val;
  a.grad = a.grad + b.grad;
};

template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator+=(dualNumber<v_t,g_t> & a, const double & b) 
{
  a.val = a.val + b;
};

template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator+=(dualNumber<v_t,g_t> & a, const float & b) 
{
  a.val = a.val + b;
};

//Decrement operator
template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator-=(dualNumber<v_t,g_t> & a, const dualNumber<v_t,g_t> & b) 
{
  a.val = a.val - b.val;
  a.grad = a.grad - b.grad;
};

template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator-=(dualNumber<v_t,g_t> & a, const double & b) 
{
  a.val = a.val - b;
};

template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator-=(dualNumber<v_t,g_t> & a, const float & b) 
{
  a.val = a.val - b;
};

//Multiplication equals operator
template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator*=(dualNumber<v_t,g_t> & a, const dualNumber<v_t,g_t> & b) 
{
  a.grad = a.val*b.grad + a.grad*b.val;
  a.val  = a.val*b.val;
};

template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator*=(dualNumber<v_t,g_t> & a, const double & b) 
{
  a.val = b*a.val;
  a.grad = b*a.grad;
};

template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator*=(dualNumber<v_t,g_t> & a, const float & b) 
{
  a.val = b*a.val;
  a.grad = b*a.grad;
};

//Divide equals operator
template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator/=(dualNumber<v_t,g_t> & a, const dualNumber<v_t,g_t> & b) 
{
  a.grad = (a.grad*b.val + a.val*b.grad)/(b.val*b.val);
  a.val  = (a.val/b.val);
};

template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator/=(dualNumber<v_t,g_t> & a, const double & b) 
{
  a.grad = a.grad/b;
  a.val  = a.val/b;
};

template<typename v_t, typename g_t>
FORCE_INLINE constexpr void operator/=(dualNumber<v_t,g_t> & a, const float & b) 
{
  a.grad = (a.grad*b.val + a.val*b.grad)/(b.val*b.val);
  a.val  = (a.val/b.val);
};

/***************************************\
!
!  Dual-Dual number operations
!
\***************************************/
//Multiplication operator
template<typename v_t, typename g_t>
FORCE_INLINE constexpr dualNumber<v_t,g_t> operator*(const dualNumber<v_t,g_t> a, const dualNumber<v_t,g_t> b) 
{
  dualNumber<v_t,g_t> newVal;
  newVal.val  = a.val*b.val;
  newVal.grad = a.val*b.grad + a.grad*b.val;
  return newVal;
};

//Division operator
template<typename v_t, typename g_t>
FORCE_INLINE constexpr dualNumber<v_t,g_t> operator/(const dualNumber<v_t,g_t> a, const dualNumber<v_t,g_t> b) 
{
  dualNumber<v_t, g_t> newVal;
  newVal.val  = (a.val/b.val);
  newVal.grad = (a.grad*b.val + a.val*b.grad)/(b.val*b.val);
  return newVal;
};

//Addition operator
template<typename v_t, typename g_t>
FORCE_INLINE constexpr dualNumber<v_t,g_t> operator+(const dualNumber<v_t,g_t> a, const dualNumber<v_t,g_t> b)
{
  dualNumber<v_t, g_t> newVal;
  newVal.val  = a.val + b.val;
  newVal.grad = a.grad + b.grad;
  return newVal;
};

//Subtraction operator
template<typename v_t, typename g_t>
FORCE_INLINE constexpr dualNumber<v_t,g_t> operator-(const dualNumber<v_t,g_t> a, const dualNumber<v_t,g_t> b)
{
  dualNumber<v_t, g_t> newVal;
  newVal.val  = a.val - b.val;
  newVal.grad = a.grad - b.grad;
  return newVal;
};

/***************************************\
!
!  Dual number comparison operations
!
\***************************************/
//equivalence operator
template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator==(const dualNumber<v_t,g_t> a, const Number b)
{
  return (a.val == b);
};

template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator==(const Number b, const dualNumber<v_t,g_t> a)
{
  return (a.val == b);
};


//inequivalence operator
template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator!=(const Number a, const dualNumber<v_t,g_t> b)
{
  return (a != b.val);
};

template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator!=(const dualNumber<v_t,g_t> a, const Number b)
{
  return (a.val != b);
};


//more than operator
template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator>(const  Number a, const dualNumber<v_t,g_t> b)
{
  return a > b.val;
};

template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator>(const dualNumber<v_t,g_t> a, const Number b)
{
  return a.val > b;
}


//less than operator
template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator<(const Number a, const dualNumber<v_t,g_t> b)
{
  return a < b.val;
};

template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator<(const dualNumber<v_t,g_t> a, const Number b)
{
  return a.val < b;
};


//more than equal operator
template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator>=(const dualNumber<v_t,g_t> a, const Number b)
{
  return a.val >= b;
};

template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator>=(const Number a, const dualNumber<v_t,g_t> b)
{
  return a >= b.val;
};


//less than equal operator
template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator<=(const dualNumber<v_t,g_t> a, const Number b)
{
  return a.val <= b;
};

template<typename v_t, typename g_t, typename Number>
FORCE_INLINE constexpr bool operator<=(const Number a, const dualNumber<v_t,g_t> b)
{
  return a <= b.val;
};


/***************************************\
!
!  Number-Dual number operations
!
\***************************************/
///////////
//Multiplication operator
///////////
template<typename val_t, typename grad_t>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator*(const dualNumber<val_t,grad_t> dNum, const double Num) 
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = Num*dNum.val;
  newVal.grad = Num*dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator*(const dualNumber<val_t,grad_t> dNum, const float Num) 
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = Num*dNum.val;
  newVal.grad = Num*dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t, typename Number>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator*(const Number Num, const dualNumber<val_t,grad_t> dNum) 
{
  return dNum*Num;
};

///////////
//Division operator
///////////
template<typename val_t, typename grad_t>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator/(const dualNumber<val_t,grad_t> dNum, const double Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = newVal.val/Num;
  newVal.grad = newVal.grad/Num;
  return newVal;
};

template<typename val_t, typename grad_t>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator/(const dualNumber<val_t,grad_t> dNum, const float Num)
{
   dualNumber<val_t,grad_t> newVal;
  newVal.val  = newVal.val/Num;
  newVal.grad = newVal.grad/Num;
  return newVal;
};

template<typename val_t, typename grad_t, typename Number>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator/(const Number Num, const dualNumber<val_t,grad_t> dNum)
{
  return Num/dNum;
};

///////////
//Addition operator
///////////
template<typename val_t, typename grad_t>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator+(const dualNumber<val_t,grad_t> dNum, const double Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = dNum.val  + Num;
  newVal.grad = dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator+(const dualNumber<val_t,grad_t> dNum, const float Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = dNum.val  + Num;
  newVal.grad = dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t, typename Number>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator+(const Number Num, const dualNumber<val_t,grad_t> dNum)
{
  return dNum + Num;
};

///////////
//Subtraction operator
///////////
template<typename val_t, typename grad_t>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator-(const dualNumber<val_t,grad_t> dNum, const double Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = dNum.val  - Num;
  newVal.grad = dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator-(const dualNumber<val_t,grad_t> dNum, const float Num)
{
  dualNumber<val_t,grad_t> newVal;
  newVal.val  = dNum.val  - Num;
  newVal.grad = dNum.grad;
  return newVal;
};

template<typename val_t, typename grad_t, typename Number>
FORCE_INLINE constexpr dualNumber<val_t,grad_t> operator-(const Number Num, const dualNumber<val_t,grad_t> dNum)
{
  return -1.0*(dNum - Num);
};
