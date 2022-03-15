/*! \file
  CUDA operators for complex numbers.
  For some reason, these aren't defined in the CUDA headers form NVIDIA

  @author R. G. Edgar
*/

#pragma once

#include "woden_precision_defs.h"

#ifdef DOUBLE_PRECISION
/*! If -DDOUBLE_PRECISION flag is added at compilation,
then cuUserComplex is set to cuDoubleComplex */
typedef cuDoubleComplex cuUserComplex;
#else
/*! If -DDOUBLE_PRECISION flag is NOT added at compilation,
then cuUserComplex is set to cuFloatComplex */
typedef cuFloatComplex cuUserComplex;
#endif

// --------------------------------
//! Unary negative
inline __device__ cuFloatComplex operator-( cuFloatComplex a ) {
  return( make_cuFloatComplex( -a.x, -a.y ) );
}

// --------------------------------
//! Self addition
inline __device__ void operator+=( cuFloatComplex &a, const cuFloatComplex b ) {
  a = cuCaddf( a, b );
}

//! Addition of two complex numbers
inline __device__ cuFloatComplex operator+( const cuFloatComplex a, const cuFloatComplex b ) {
  return( cuCaddf( a, b ) );
}

//! Addition of float to complex
inline __device__ cuFloatComplex operator+( const cuFloatComplex a, const float r ) {
  return( make_cuFloatComplex( a.x+r, a.y ) );
}

//! Addition of complex to float
inline __device__ cuFloatComplex operator+( const float r, const cuFloatComplex z ) {
  return( z + r );
}

// --------------------------------
//! Self subtraction
inline __device__ void operator-=( cuFloatComplex &a, const cuFloatComplex b ) {
  a = cuCsubf( a, b );
}

//! Subtraction
inline __device__ cuFloatComplex operator-( const cuFloatComplex a, const cuFloatComplex b ) {
  return( cuCsubf( a, b ) );
}


// --------------------------------
//! Self multiplication
inline __device__ void operator*=( cuFloatComplex &a, const cuFloatComplex b ) {
  a = cuCmulf( a, b );
}

//! Multiplication of two complex numbers
inline __device__ cuFloatComplex operator*( const cuFloatComplex a, const cuFloatComplex b ) {
  return( cuCmulf( a, b ) );
}

// --------------------------------
//! Self multiplication by real number
inline __device__ void operator*=( cuFloatComplex &z, const float r ) {
  z.x *= r;
  z.y *= r;
}

//! Multiplication by a real number
inline __device__ cuFloatComplex operator*( const cuFloatComplex z, const float r ) {
  cuFloatComplex temp;

  temp=z;
  temp *= r;
  return( temp );
}

//! Multiplication of a real number
inline __device__ cuFloatComplex operator*( const float r, const cuFloatComplex z ) {
  return( z*r );
}


// --------------------------------
//! Division of two complex numbers
inline __device__ cuFloatComplex operator/( const cuFloatComplex a, const cuFloatComplex b ) {
  return( cuCdivf( a, b ) );
}

//! Division of a real number
inline __device__ cuFloatComplex operator/( const float r, const cuFloatComplex b ) {
  return( cuCdivf( make_cuFloatComplex( r, 0 ), b ) );
}

//! Division by a real number
inline __device__ cuFloatComplex operator/( const cuFloatComplex a, const float r ) {
  return( make_cuFloatComplex( a.x / r, a.y / r ) );
}


// ----------------------------------------------------------
//! Exponentiation
//! Computes e^z == e^x ( cos y + i sin y )
inline __device__ cuFloatComplex cuComplexExp( const cuFloatComplex z ) {

  float x = cuCrealf( z );
  float y = cuCimagf( z );

  cuFloatComplex temp = make_cuFloatComplex( cosf(y), sinf(y) );

  return( (expf(x)) * temp );
}

// ----------------------------------------------------------
//! Calculate real and imaginary parts of U(1) variable e^(i*theta)
inline __device__ cuFloatComplex U1polar( const float theta ) {

  cuFloatComplex z;
  sincosf(theta, &(z.y), &(z.x));

  return z;
}





//Do everything again, but for double precision

// --------------------------------
//! Unary negative
inline __device__ cuDoubleComplex operator-( cuDoubleComplex a ) {
  return( make_cuDoubleComplex( -a.x, -a.y ) );
}

// --------------------------------
//! Self addition
inline __device__ void operator+=( cuDoubleComplex &a, const cuDoubleComplex b ) {
  a = cuCadd( a, b );
}

//! Addition of two complex numbers
inline __device__ cuDoubleComplex operator+( const cuDoubleComplex a, const cuDoubleComplex b ) {
  return( cuCadd( a, b ) );
}

//! Addition of double to complex
inline __device__ cuDoubleComplex operator+( const cuDoubleComplex a, const double r ) {
  return( make_cuDoubleComplex( a.x+r, a.y ) );
}

//! Addition of complex to double
inline __device__ cuDoubleComplex operator+( const double r, const cuDoubleComplex z ) {
  return( z + r );
}

// --------------------------------
//! Self subtraction
inline __device__ void operator-=( cuDoubleComplex &a, const cuDoubleComplex b ) {
  a = cuCsub( a, b );
}

//! Subtraction
inline __device__ cuDoubleComplex operator-( const cuDoubleComplex a, const cuDoubleComplex b ) {
  return( cuCsub( a, b ) );
}


// --------------------------------
//! Self multiplication
inline __device__ void operator*=( cuDoubleComplex &a, const cuDoubleComplex b ) {
  a = cuCmul( a, b );
}

//! Multiplication of two complex numbers
inline __device__ cuDoubleComplex operator*( const cuDoubleComplex a, const cuDoubleComplex b ) {
  return( cuCmul( a, b ) );
}

// --------------------------------
//! Self multiplication by real number
inline __device__ void operator*=( cuDoubleComplex &z, const double r ) {
  z.x *= r;
  z.y *= r;
}

//! Multiplication by a real number
inline __device__ cuDoubleComplex operator*( const cuDoubleComplex z, const double r ) {
  cuDoubleComplex temp;

  temp=z;
  temp *= r;
  return( temp );
}

//! Multiplication of a real number
inline __device__ cuDoubleComplex operator*( const double r, const cuDoubleComplex z ) {
  return( z*r );
}

// --------------------------------
//! Division of two complex numbers
inline __device__ cuDoubleComplex operator/( const cuDoubleComplex a, const cuDoubleComplex b ) {
  return( cuCdiv( a, b ) );
}

//! Division of a real number
inline __device__ cuDoubleComplex operator/( const double r, const cuDoubleComplex b ) {
  return( cuCdiv( make_cuDoubleComplex( r, 0 ), b ) );
}

//! Division by a real number
inline __device__ cuDoubleComplex operator/( const cuDoubleComplex a, const double r ) {
  return( make_cuDoubleComplex( a.x / r, a.y / r ) );
}


// ----------------------------------------------------------
//! Exponentiation
//! Computes e^z == e^x ( cos y + i sin y )
inline __device__ cuDoubleComplex cuComplexExp( const cuDoubleComplex z ) {

  double x = cuCreal( z );
  double y = cuCimag( z );

  cuDoubleComplex temp = make_cuDoubleComplex( cos(y), sin(y) );

  return( (exp(x)) * temp );
}

// ----------------------------------------------------------
//! Calculate real and imaginary parts of U(1) variable e^(i*theta)
inline __device__ cuDoubleComplex U1polar( const double theta ) {

  cuDoubleComplex z;
  sincos(theta, &(z.y), &(z.x));

  return z;
}

//Make a CUDA complex to the precision set during compilation
inline __device__ cuUserComplex make_cuUserComplex( user_precision_t real,
                                                    user_precision_t imag) {

  cuUserComplex z;
  z.x = real;
  z.y = imag;

  return z;
}
