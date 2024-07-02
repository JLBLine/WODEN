/*! \file
  CUDA operators for complex numbers.
  For some reason, these aren't defined in the CUDA headers form NVIDIA

  @author R. G. Edgar; modified with macros by M. Sokolowski
*/

#pragma once

#include "woden_precision_defs.h"
#include "gpu_macros.h"

#ifdef DOUBLE_PRECISION
/*! If -DDOUBLE_PRECISION flag is added at compilation,
then cuUserComplex is set to cuDoubleComplex */
typedef gpuDoubleComplex gpuUserComplex;
#else
/*! If -DDOUBLE_PRECISION flag is NOT added at compilation,
then cuUserComplex is set to cuFloatComplex */
typedef gpuFloatComplex gpuUserComplex;
#endif

// --------------------------------
//! Unary negative
inline __device__ gpuFloatComplex operator-( gpuFloatComplex a ) {
  return( make_gpuFloatComplex( -a.x, -a.y ) );
}

// --------------------------------
//! Self addition
inline __device__ void operator+=( gpuFloatComplex &a, const gpuFloatComplex b ) {
  a = gpuCaddf( a, b );
}

//! Addition of two complex numbers
inline __device__ gpuFloatComplex operator+( const gpuFloatComplex a, const gpuFloatComplex b ) {
  return( gpuCaddf( a, b ) );
}

//! Addition of float to complex
inline __device__ gpuFloatComplex operator+( const gpuFloatComplex a, const float r ) {
  return( make_gpuFloatComplex( a.x+r, a.y ) );
}

//! Addition of complex to float
inline __device__ gpuFloatComplex operator+( const float r, const gpuFloatComplex z ) {
  return( z + r );
}

// --------------------------------
//! Self subtraction
inline __device__ void operator-=( gpuFloatComplex &a, const gpuFloatComplex b ) {
  a = gpuCsubf( a, b );
}

//! Subtraction
inline __device__ gpuFloatComplex operator-( const gpuFloatComplex a, const gpuFloatComplex b ) {
  return( gpuCsubf( a, b ) );
}


// --------------------------------
//! Self multiplication
inline __device__ void operator*=( gpuFloatComplex &a, const gpuFloatComplex b ) {
  a = gpuCmulf( a, b );
}

//! Multiplication of two complex numbers
inline __device__ gpuFloatComplex operator*( const gpuFloatComplex a, const gpuFloatComplex b ) {
  return( gpuCmulf( a, b ) );
}

// --------------------------------
//! Self multiplication by real number
inline __device__ void operator*=( gpuFloatComplex &z, const float r ) {
  z.x *= r;
  z.y *= r;
}

//! Multiplication by a real number
inline __device__ gpuFloatComplex operator*( const gpuFloatComplex z, const float r ) {
  gpuFloatComplex temp;

  temp=z;
  temp *= r;
  return( temp );
}

//! Multiplication of a real number
inline __device__ gpuFloatComplex operator*( const float r, const gpuFloatComplex z ) {
  return( z*r );
}


// --------------------------------
//! Division of two complex numbers
inline __device__ gpuFloatComplex operator/( const gpuFloatComplex a, const gpuFloatComplex b ) {
  return( gpuCdivf( a, b ) );
}

//! Division of a real number
inline __device__ gpuFloatComplex operator/( const float r, const gpuFloatComplex b ) {
  return( gpuCdivf( make_gpuFloatComplex( r, 0 ), b ) );
}

//! Division by a real number
inline __device__ gpuFloatComplex operator/( const gpuFloatComplex a, const float r ) {
  return( make_gpuFloatComplex( a.x / r, a.y / r ) );
}


// ----------------------------------------------------------
//! Exponentiation
//! Computes e^z == e^x ( cos y + i sin y )
inline __device__ gpuFloatComplex gpuComplexExp( const gpuFloatComplex z ) {

  float x = gpuCrealf( z );
  float y = gpuCimagf( z );

  gpuFloatComplex temp = make_gpuFloatComplex( cosf(y), sinf(y) );

  return( (expf(x)) * temp );
}

// ----------------------------------------------------------
//! Calculate real and imaginary parts of U(1) variable e^(i*theta)
inline __device__ gpuFloatComplex U1polar( const float theta ) {

  gpuFloatComplex z;
  sincosf(theta, &(z.y), &(z.x));

  return z;
}





//Do everything again, but for double precision

// --------------------------------
//! Unary negative
inline __device__ gpuDoubleComplex operator-( gpuDoubleComplex a ) {
  return( make_gpuDoubleComplex( -a.x, -a.y ) );
}

// --------------------------------
//! Self addition
inline __device__ void operator+=( gpuDoubleComplex &a, const gpuDoubleComplex b ) {
  a = gpuCadd( a, b );
}

//! Addition of two complex numbers
inline __device__ gpuDoubleComplex operator+( const gpuDoubleComplex a, const gpuDoubleComplex b ) {
  return( gpuCadd( a, b ) );
}

//! Addition of double to complex
inline __device__ gpuDoubleComplex operator+( const gpuDoubleComplex a, const double r ) {
  return( make_gpuDoubleComplex( a.x+r, a.y ) );
}

//! Addition of complex to double
inline __device__ gpuDoubleComplex operator+( const double r, const gpuDoubleComplex z ) {
  return( z + r );
}

// --------------------------------
//! Self subtraction
inline __device__ void operator-=( gpuDoubleComplex &a, const gpuDoubleComplex b ) {
  a = gpuCsub( a, b );
}

//! Subtraction
inline __device__ gpuDoubleComplex operator-( const gpuDoubleComplex a, const gpuDoubleComplex b ) {
  return( gpuCsub( a, b ) );
}


// --------------------------------
//! Self multiplication
inline __device__ void operator*=( gpuDoubleComplex &a, const gpuDoubleComplex b ) {
  a = gpuCmul( a, b );
}

//! Multiplication of two complex numbers
inline __device__ gpuDoubleComplex operator*( const gpuDoubleComplex a, const gpuDoubleComplex b ) {
  return( gpuCmul( a, b ) );
}

// --------------------------------
//! Self multiplication by real number
inline __device__ void operator*=( gpuDoubleComplex &z, const double r ) {
  z.x *= r;
  z.y *= r;
}

//! Multiplication by a real number
inline __device__ gpuDoubleComplex operator*( const gpuDoubleComplex z, const double r ) {
  gpuDoubleComplex temp;

  temp=z;
  temp *= r;
  return( temp );
}

//! Multiplication of a real number
inline __device__ gpuDoubleComplex operator*( const double r, const gpuDoubleComplex z ) {
  return( z*r );
}

// --------------------------------
//! Division of two complex numbers
inline __device__ gpuDoubleComplex operator/( const gpuDoubleComplex a, const gpuDoubleComplex b ) {
  return( gpuCdiv( a, b ) );
}

//! Division of a real number
inline __device__ gpuDoubleComplex operator/( const double r, const gpuDoubleComplex b ) {
  return( gpuCdiv( make_gpuDoubleComplex( r, 0 ), b ) );
}

//! Division by a real number
inline __device__ gpuDoubleComplex operator/( const gpuDoubleComplex a, const double r ) {
  return( make_gpuDoubleComplex( a.x / r, a.y / r ) );
}


// ----------------------------------------------------------
//! Exponentiation
//! Computes e^z == e^x ( cos y + i sin y )
inline __device__ gpuDoubleComplex gpuComplexExp( const gpuDoubleComplex z ) {

  double x = gpuCreal( z );
  double y = gpuCimag( z );

  gpuDoubleComplex temp = make_gpuDoubleComplex( cos(y), sin(y) );

  return( (exp(x)) * temp );
}

// ----------------------------------------------------------
//! Calculate real and imaginary parts of U(1) variable e^(i*theta)
inline __device__ gpuDoubleComplex U1polar( const double theta ) {

  gpuDoubleComplex z;
  sincos(theta, &(z.y), &(z.x));

  return z;
}

//! Make a CUDA complex to the precision set during compilation
inline __device__ gpuUserComplex make_gpuUserComplex( user_precision_t real,
                                                    user_precision_t imag) {

  gpuUserComplex z;
  z.x = real;
  z.y = imag;

  return z;
}
