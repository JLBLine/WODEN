/*! \file
  CUDA operators for complex numbers.
  For some reason, these aren't defined in the CUDA headers form NVIDIA

  @author R. G. Edgar
*/

#ifndef CUDACOMPLEX_H
#define CUDACOMPLEX_H

// #include "rts_common.h"

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

//CUBE_DEVICE(cuFloatComplex, cuComplexExp, const cuFloatComplex z ) {
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

#endif
