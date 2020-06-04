# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "legendre_polynomial.h"

/******************************************************************************/

char digit_to_ch ( int digit )

/******************************************************************************/
/*
  Purpose:

    DIGIT_TO_CH returns the character representation of a decimal digit.

  Example:

    DIGIT   C
    -----  ---
      0    '0'
      1    '1'
    ...    ...
      9    '9'
     17    '*'

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int DIGIT, the digit value between 0 and 9.

    Output, char DIGIT_TO_CH, the corresponding character, or '*' if DIGIT
    was illegal.
*/
{
  if ( 0 <= digit && digit <= 9 )
  {
    return ( digit + 48 );
  }
  else
  {
    return '*';
  }
}
/******************************************************************************/

int i4_log_10 ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.

  Example:

        I  I4_LOG_10
    -----  --------
        0    0
        1    0
        2    0
        9    0
       10    1
       11    1
       99    1
      100    2
      101    2
      999    2
     1000    3
     1001    3
     9999    3
    10000    4

  Discussion:

    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number whose logarithm base 10 is desired.

    Output, int I4_LOG_10, the integer part of the logarithm base 10 of
    the absolute value of X.
*/
{
  int i_abs;
  int ten_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    ten_pow = 10;

    i_abs = abs ( i );

    while ( ten_pow <= i_abs )
    {
      value = value + 1;
      ten_pow = ten_pow * 10;
    }
  }
  return value;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

char *i4_to_s ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_TO_S converts an I4 to a string.

  Example:

    INTVAL  S

         1  1
        -1  -1
         0  0
      1952  1952
    123456  123456
   1234567  1234567

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int I, an integer to be converted.

    Output, char *I4_TO_S, the representation of the integer.
*/
{
  int digit;
  int j;
  int length;
  int ten_power;
  char *s;
  static double ten = 10.0;

  length = i4_log_10 ( i );

  ten_power = ( int ) ( pow ( ten, length ) );

  if ( i < 0 )
  {
    length = length + 1;
  }
/*
  Add one position for the trailing null.
*/
  length = length + 1;

  s = ( char * ) malloc ( length * sizeof ( char ) );

  if ( i == 0 )
  {
    s[0] = '0';
    s[1] = '\0';
    return s;
  }
/*
  Now take care of the sign.
*/
  j = 0;
  if ( i < 0 )
  {
    s[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
/*
  Find the leading digit of I, strip it off, and stick it into the string.
*/
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
/*
  Tack on the trailing NULL.
*/
  s[j] = '\0';
  j = j + 1;

  return s;
}
/******************************************************************************/

void imtqlx ( int n, double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    IMTQLX diagonalizes a symmetric tridiagonal matrix.

  Discussion:

    This routine is a slightly modified version of the EISPACK routine to
    perform the implicit QL algorithm on a symmetric tridiagonal matrix.

    The authors thank the authors of EISPACK for permission to use this
    routine.

    It has been modified to produce the product Q' * Z, where Z is an input
    vector and Q is the orthogonal matrix diagonalizing the input matrix.
    The changes consist (essentialy) of applying the orthogonal transformations
    directly to Z as they are generated.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

    Roger Martin, James Wilkinson,
    The Implicit QL Algorithm,
    Numerische Mathematik,
    Volume 12, Number 5, December 1968, pages 377-383.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D(N), the diagonal entries of the matrix.
    On output, the information in D has been overwritten.

    Input/output, double E(N), the subdiagonal entries of the
    matrix, in entries E(1) through E(N-1).  On output, the information in
    E has been overwritten.

    Input/output, double Z(N).  On input, a vector.  On output,
    the value of Q' * Z, where Q is the matrix that diagonalizes the
    input symmetric tridiagonal matrix.
*/
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( fabs ( e[m-1] ) <= prec * ( fabs ( d[m-1] ) + fabs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "IMTQLX - Fatal error!\n" );
        fprintf ( stderr, "  Iteration limit exceeded\n" );
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r =  sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + fabs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( fabs ( g ) <= fabs ( f ) )
        {
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
/*
  Sorting.
*/
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
/******************************************************************************/

double *p_exponential_product ( int p, double b )

/******************************************************************************/
/*
  Purpose:

    P_EXPONENTIAL_PRODUCT: exponential products for P(n,x).

  Discussion:

    Let P(n,x) represent the Legendre polynomial of degree n.

    For polynomial chaos applications, it is of interest to know the
    value of the integrals of products of exp(B*X) with every possible pair
    of basis functions.  That is, we'd like to form

      Tij = Integral ( -1.0 <= X <= +1.0 ) exp(B*X) * P(I,X) * P(J,X) dx

    We will estimate these integrals using Gauss-Legendre quadrature.
    Because of the exponential factor exp(B*X), the quadrature will not
    be exact.

    However, when B = 0, the quadrature is exact, and moreoever, the
    table will be the identity matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int P, the maximum degree of the polyonomial
    factors.  0 <= P.

    Input, double B, the coefficient of X in the exponential factor.

    Output, double P_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of integrals.
*/
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = ( double * ) malloc ( (p+1)*(p+1) * sizeof ( double ) );;

  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = ( 3 * p + 4 ) / 2;

  x_table = ( double * ) malloc ( order * sizeof ( double ) );;
  w_table = ( double * ) malloc ( order * sizeof ( double ) );;

  p_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = p_polynomial_value ( 1, p, x_table + k );
/*
  The following formula is an outer product in H_TABLE.
*/
    for ( j = 0; j <= p; j++ )
    {
      for ( i = 0; i <= p; i++ )
      {
        table[i+j*(p+1)] = table[i+j*(p+1)]
          + w_table[k] * exp ( b * x ) * h_table[i] * h_table[j];
      }
    }
    free ( h_table );
  }

  free ( w_table );
  free ( x_table );

  return table;
}
/******************************************************************************/

double p_integral ( int n )

/******************************************************************************/
/*
  Purpose:

    P_INTEGRAL evaluates a monomial integral associated with P(n,x).

  Discussion:

    The integral:

      integral ( -1 <= x < +1 ) x^n dx

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the exponent.
    0 <= N.

    Output, double P_INTEGRAL, the value of the integral.
*/
{
  double value;

  if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = 2.0 / ( double ) ( n + 1 );
  }
  return value;
}
/******************************************************************************/

double *p_polynomial_coefficients ( int n )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomials P(n,x).

  First terms:

     1
     0     1
    -1/2   0      3/2
     0    -3/2    0     5/2
     3/8   0    -30/8   0     35/8
     0    15/8    0   -70/8    0     63/8
    -5/16  0    105/16  0   -315/16   0    231/16
     0   -35/16   0   315/16   0   -693/16   0    429/16

     1.00000
     0.00000  1.00000
    -0.50000  0.00000  1.50000
     0.00000 -1.50000  0.00000  2.5000
     0.37500  0.00000 -3.75000  0.00000  4.37500
     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.

  Parameters:

    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.

    Output, double P_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients of
    the Legendre polynomials of degree 0 through N.
*/
{
  double *c;
  int i;
  int j;

  if ( n < 0 )
  {
    return NULL;
  }

  c = ( double * ) malloc ( (n+1)*(n+1) * sizeof ( double ) );;

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n <= 0 )
  {
    return c;
  }

  c[1+1*(n+1)] = 1.0;

  for ( i = 2; i <= n; i++ )
  {
    for ( j = 0; j <= i-2; j++ )
    {
      c[i+j*(n+1)] =
          ( double ) (   - i + 1 ) * c[i-2+j*(n+1)] / ( double ) i;
    }
    for ( j = 1; j <= i; j++ )
    {
      c[i+j*(n+1)] = c[i+j*(n+1)]
        + ( double ) ( i + i - 1 ) * c[i-1+(j-1)*(n+1)] / ( double ) i;
    }
  }

  return c;
}
/******************************************************************************/

double *p_polynomial_prime ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_PRIME evaluates the derivative of Legendre polynomials P(n,x).

  Discussion:

    P(0,X) = 1
    P(1,X) = X
    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N

    P'(0,X) = 0
    P'(1,X) = 1
    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N

    Thanks to Dimitriy Morozov for pointing out a memory leak caused by
    not deleting the work array V before return, 19 March 2013.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.

  Parameters:

    Input, int M, the number of evaluation points.

    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.

    Input, double X[M], the evaluation points.

    Output, double P_POLYNOMIAL_PRIME[M*(N+1)], the values of the derivatives
    of the Legendre polynomials of order 0 through N at the points.
*/
{
  int i;
  int j;
  double *v;
  double *vp;

  if ( n < 0 )
  {
    return NULL;
  }

  vp = ( double * ) malloc ( m*(n+1) * sizeof ( double ) );;

  for ( i = 0; i < m; i++ )
  {
    vp[i+0*m] = 0.0;
  }

  if ( n < 1 )
  {
    return vp;
  }

  v = ( double * ) malloc ( m*(n+1) * sizeof ( double ) );;

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
    vp[i+1*m] = 1.0;
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = ( ( double ) ( 2 * j - 1 ) * x[i] * v[i+(j-1)*m]
                 - ( double ) (     j - 1 ) *        v[i+(j-2)*m] )
                 / ( double ) (     j     );

      vp[i+j*m] = ( ( double ) ( 2 * j - 1 ) * ( v[i+(j-1)*m] + x[i] * vp[i+(j-1)*m] )
                  - ( double ) (     j - 1 ) *   vp[i+(j-2)*m]               )
                  / ( double ) (     j     );
    }
  }

  free ( v );

  return vp;
}
/******************************************************************************/

double *p_polynomial_prime2 ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_PRIME2: second derivative of Legendre polynomials P(n,x).

  Discussion:

    P(0,X) = 1
    P(1,X) = X
    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N

    P'(0,X) = 0
    P'(1,X) = 1
    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N

    P"(0,X) = 0
    P"(1,X) = 0
    P"(N,X) = ( (2*N-1)*(2*P'(N-1,X)+X*P"(N-1,X)-(N-1)*P"(N-2,X) ) / N

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.

  Parameters:

    Input, int M, the number of evaluation points.

    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.

    Input, double X[M], the evaluation points.

    Output, double P_POLYNOMIAL_PRIME2[M*(N+1)], the second derivative
    of the Legendre polynomials of order 0 through N at the points.
*/
{
  int i;
  int j;
  double *v;
  double *vp;
  double *vpp;

  if ( n < 0 )
  {
    return NULL;
  }

  vpp = ( double * ) malloc ( m*(n+1) * sizeof ( double ) );;

  for ( i = 0; i < m; i++ )
  {
    vpp[i+0*m] = 0.0;
  }

  if ( n < 1 )
  {
    return vpp;
  }

  v = ( double * ) malloc ( m*(n+1) * sizeof ( double ) );;
  vp = ( double * ) malloc ( m*(n+1) * sizeof ( double ) );;

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
    vp[i+0*m] = 0.0;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
    vp[i+1*m] = 1.0;
    vpp[i+1*m] = 0.0;
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] =
        ( ( double ) ( 2 * j - 1 ) * x[i] * v[i+(j-1)*m]
        - ( double ) (     j - 1 ) *        v[i+(j-2)*m] )
        / ( double ) (     j     );

      vp[i+j*m] =
        ( ( double ) ( 2 * j - 1 ) * ( v[i+(j-1)*m] + x[i] * vp[i+(j-1)*m] )
        - ( double ) (     j - 1 ) *   vp[i+(j-2)*m]               )
        / ( double ) (     j     );

      vpp[i+j*m] =
        ( ( double ) ( 2 * j - 1 ) * ( 2.0 * vp[i+(j-1)*m] + x[i] * vpp[i+(j-1)*m] )
        - ( double ) (     j - 1 ) *   vpp[i+(j-2)*m]               )
        / ( double ) (     j     );
    }
  }

  free ( v );
  free ( vp );

  return vpp;
}
/******************************************************************************/

double *p_polynomial_value ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).

  Discussion:

    P(n,1) = 1.
    P(n,-1) = (-1)^N.
    | P(n,x) | <= 1 in [-1,1].

    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
    quadrature of the integral of a function F(X) with weight function 1
    over the interval [-1,1].

    The Legendre polynomials are orthogonal under the inner product defined
    as integration from -1 to 1:

      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
        = 0 if I =/= J
        = 2 / ( 2*I+1 ) if I = J.

    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.

    A function F(X) defined on [-1,1] may be approximated by the series
      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
    where
      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.

    The formula is:

      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)

  Differential equation:

    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0

  First terms:

    P( 0,x) =      1
    P( 1,x) =      1 X
    P( 2,x) = (    3 X^2 -       1)/2
    P( 3,x) = (    5 X^3 -     3 X)/2
    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256

  Recursion:

    P(0,x) = 1
    P(1,x) = x
    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n

    P'(0,x) = 0
    P'(1,x) = 1
    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.

  Parameters:

    Input, int M, the number of evaluation points.

    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.

    Input, double X[M], the evaluation points.

    Output, double P_POLYNOMIAL_VALUE[M*(N+1)], the values of the Legendre
    polynomials of order 0 through N.
*/
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = ( double * ) malloc ( m*(n+1) * sizeof ( double ) );;

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = ( ( double ) ( 2 * j - 1 ) * x[i] * v[i+(j-1)*m]
                 - ( double ) (     j - 1 ) *        v[i+(j-2)*m] )
                 / ( double ) (     j     );
    }
  }

  return v;
}
/******************************************************************************/

void p_polynomial_values ( int *n_data, int *n, double *x, double *fx )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_VALUES returns values of the Legendre polynomials P(n,x).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, int *N, the order of the function.

    Output, double *X, the point where the function is evaluated.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.2500000000000000E+00,
     -0.4062500000000000E+00,
     -0.3359375000000000E+00,
      0.1577148437500000E+00,
      0.3397216796875000E+00,
      0.2427673339843750E-01,
     -0.2799186706542969E+00,
     -0.1524540185928345E+00,
      0.1768244206905365E+00,
      0.2212002165615559E+00,
      0.0000000000000000E+00,
     -0.1475000000000000E+00,
     -0.2800000000000000E+00,
     -0.3825000000000000E+00,
     -0.4400000000000000E+00,
     -0.4375000000000000E+00,
     -0.3600000000000000E+00,
     -0.1925000000000000E+00,
      0.8000000000000000E-01,
      0.4725000000000000E+00,
      0.1000000000000000E+01 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10,  3,
     3,  3,  3,
     3,  3,  3,
     3,  3,  3,
     3 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.00E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.90E+00,
     1.00E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

double *p_polynomial_zeros ( int nt )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_ZEROS: zeros of Legendre function P(n,x).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NT, the order of the rule.

    Output, double P_POLYNOMIAL_ZEROS[NT], the zeros.
*/
{
  double *bj;
  int i;
  double *t;
  double *wts;

  t = ( double * ) malloc ( nt * sizeof ( double ) );;

  for ( i = 0; i < nt; i++ )
  {
    t[i] = 0.0;
  }

  bj = ( double * ) malloc ( nt * sizeof ( double ) );;
  for ( i = 0; i < nt; i++ )
  {
    bj[i] = ( double ) ( ( i + 1 ) * ( i + 1 ) )
          / ( double ) ( 4 *  ( i + 1 ) * ( i + 1 ) - 1 );
    bj[i] = sqrt ( bj[i] );
  }

  wts = ( double * ) malloc ( nt * sizeof ( double ) );;
  wts[0] = sqrt ( 2.0 );
  for ( i = 1; i < nt; i++ )
  {
    wts[i] = 0.0;
  }

  imtqlx ( nt, t, bj, wts );

  free ( bj );
  free ( wts );

  return t;
}
/******************************************************************************/

double *p_power_product ( int p, int e )

/******************************************************************************/
/*
  Purpose:

    P_POWER_PRODUCT: power products for Legendre polynomial P(n,x).

  Discussion:

    Let P(n,x) represent the Legendre polynomial of degree n.

    For polynomial chaos applications, it is of interest to know the
    value of the integrals of products of X with every possible pair
    of basis functions.  That is, we'd like to form

      Tij = Integral ( -1.0 <= X <= +1.0 ) X^E * P(I,x) * P(J,x) dx

    We will estimate these integrals using Gauss-Legendre quadrature.

    When E is 0, the computed table should be the identity matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int P, the maximum degree of the polyonomial
    factors.  0 <= P.

    Input, int E, the exponent of X in the integrand.
    0 <= E.

    Output, double P_POWER_PRODUCT[(P+1)*(P+1)], the table of integrals.
*/
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = ( double * ) malloc ( (p+1)*(p+1) * sizeof ( double ) );;

  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = p + 1 + ( ( e + 1 ) / 2 );

  x_table = ( double * ) malloc ( order * sizeof ( double ) );
  w_table = ( double * ) malloc ( order * sizeof ( double ) );

  p_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = p_polynomial_value ( 1, p, x_table + k );
/*
  The following formula is an outer product in H_TABLE.
*/
    if ( e == 0 )
    {
      for ( i = 0; i <= p; i++ )
      {
        for ( j = 0; j <= p; j++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] + w_table[k] * h_table[i] * h_table[j];
        }
      }
    }
    else
    {
      for ( i = 0; i <= p; i++ )
      {
        for ( j = 0; j <= p; j++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)]
            + w_table[k] * pow ( x, e ) * h_table[i] * h_table[j];
        }
      }
    }
    free ( h_table );
  }

  free ( w_table );
  free ( x_table );

  return table;
}
/******************************************************************************/

void p_quadrature_rule ( int nt, double t[], double wts[] )

/******************************************************************************/
/*
  Purpose:

    P_QUADRATURE_RULE: quadrature for Legendre function P(n,x).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NT, the order of the rule.

    Output, double T[NT], WTS[NT], the points and weights
    of the rule.
*/
{
  double *bj;
  int i;

  for ( i = 0; i < nt; i++ )
  {
    t[i] = 0.0;
  }

  bj = ( double * ) malloc ( nt * sizeof ( double ) );;
  for ( i = 0; i < nt; i++ )
  {
    bj[i] = ( double ) ( ( i + 1 ) * ( i + 1 ) )
          / ( double ) ( 4 *  ( i + 1 ) * ( i + 1 ) - 1 );
    bj[i] = sqrt ( bj[i] );
  }

  wts[0] = sqrt ( 2.0 );
  for ( i = 1; i < nt; i++ )
  {
    wts[i] = 0.0;
  }

  imtqlx ( nt, t, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = pow ( wts[i], 2 );
  }
  free ( bj );

  return;
}
/******************************************************************************/

double *pm_polynomial_value ( int mm, int n, int m, double x[] )

/******************************************************************************/
/*
  Purpose:

    PM_POLYNOMIAL_VALUE evaluates the Legendre polynomials Pm(n,m,x).

  Differential equation:

    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0

  First terms:

    M = 0  ( = Legendre polynomials of first kind P(N,X) )

    Pm(0,0,x) =    1
    Pm(1,0,x) =    1 X
    Pm(2,0,x) = (  3 X^2 -   1)/2
    Pm(3,0,x) = (  5 X^3 -   3 X)/2
    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16

    M = 1

    Pm(0,1,x) =   0
    Pm(1,1,x) =   1 * SQRT(1-X^2)
    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)

    M = 2

    Pm(0,2,x) =   0
    Pm(1,2,x) =   0
    Pm(2,2,x) =   3 * (1-X^2)
    Pm(3,2,x) =  15 * (1-X^2) * X
    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)

    M = 3

    Pm(0,3,x) =   0
    Pm(1,3,x) =   0
    Pm(2,3,x) =   0
    Pm(3,3,x) =  15 * (1-X^2)^1.5
    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X

    M = 4

    Pm(0,4,x) =   0
    Pm(1,4,x) =   0
    Pm(2,4,x) =   0
    Pm(3,4,x) =   0
    Pm(4,4,x) = 105 * (1-X^2)^2

  Recursion:

    if N < M:
      Pm(N,M,x) = 0
    if N = M:
      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
      all the odd integers less than or equal to N.
    if N = M+1:
      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
    if M+1 < N:
      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

  Parameters:

    Input, int MM, the number of evaluation points.

    Input, int N, the maximum first index of the Legendre
    function, which must be at least 0.

    Input, int M, the second index of the Legendre function,
    which must be at least 0, and no greater than N.

    Input, double X[MM], the point at which the function is to be
    evaluated.

    Output, double PM_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
*/
{
  double fact;
  int i;
  int j;
  int k;
  double *v;

  v = ( double * ) malloc ( mm*(n+1) * sizeof ( double ) );;

  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = 0.0;
    }
  }
/*
  J = M is the first nonzero function.
*/
  if ( m <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+m*mm] = 1.0;
    }

    fact = 1.0;
    for ( k = 0; k < m; k++ )
    {
      for ( i = 0; i < mm; i++ )
      {
        v[i+m*mm] = - v[i+m*mm] * fact * sqrt ( 1.0 - x[i] * x[i] );
      }
      fact = fact + 2.0;
    }
  }
/*
  J = M + 1 is the second nonzero function.
*/
  if ( m + 1 <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+(m+1)*mm] = x[i] * ( double ) ( 2 * m + 1 ) * v[i+m*mm];
    }
  }
/*
  Now we use a three term recurrence.
*/
  for ( j = m + 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = ( ( double ) ( 2 * j     - 1 ) * x[i] * v[i+(j-1)*mm]
                  + ( double ) (   - j - m + 1 ) *        v[i+(j-2)*mm] )
                  / ( double ) (     j - m     );
    }
  }

  return v;
}
/******************************************************************************/

void pm_polynomial_values ( int *n_data, int *n, int *m, double *x, double *fx )

/******************************************************************************/
/*
  Purpose:

    PM_POLYNOMIAL_VALUES returns values of Legendre polynomials Pm(n,m,x).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, int *N, int *M, double *X,
    the arguments of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.0000000000000000E+00,
     -0.5000000000000000E+00,
      0.0000000000000000E+00,
      0.3750000000000000E+00,
      0.0000000000000000E+00,
     -0.8660254037844386E+00,
     -0.1299038105676658E+01,
     -0.3247595264191645E+00,
      0.1353164693413185E+01,
     -0.2800000000000000E+00,
      0.1175755076535925E+01,
      0.2880000000000000E+01,
     -0.1410906091843111E+02,
     -0.3955078125000000E+01,
     -0.9997558593750000E+01,
      0.8265311444100484E+02,
      0.2024442836815152E+02,
     -0.4237997531890869E+03,
      0.1638320624828339E+04,
     -0.2025687389227225E+05  };

  static int m_vec[N_MAX] = {
    0, 0, 0, 0,
    0, 1, 1, 1,
    1, 0, 1, 2,
    3, 2, 2, 3,
    3, 4, 4, 5 };

  static int n_vec[N_MAX] = {
    1,  2,  3,  4,
    5,  1,  2,  3,
    4,  3,  3,  3,
    3,  4,  5,  6,
    7,  8,  9, 10 };

  static double x_vec[N_MAX] = {
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *m = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *m = m_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

double *pmn_polynomial_value ( int mm, int n, int m, double x[] )

/******************************************************************************/
/*
  Purpose:

    PMN_POLYNOMIAL_VALUE: normalized Legendre polynomial Pmn(n,m,x).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

  Parameters:

    Input, int MM, the number of evaluation points.

    Input, int N, the maximum first index of the Legendre
    function, which must be at least 0.

    Input, int M, the second index of the Legendre function,
    which must be at least 0, and no greater than N.

    Input, double X[MM], the evaluation points.

    Output, double PMN_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
*/
{
  double factor;
  int i;
  int j;
  double *v;

  v = pm_polynomial_value ( mm, n, m, x );
/*
  Normalization.
*/
  for ( j = m; j <= n; j++ )
  {
    factor = sqrt ( ( ( double ) ( 2 * j + 1 ) * r8_factorial ( j - m ) )
      / ( 2.0 * r8_factorial ( j + m ) ) );
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = v[i+j*mm] * factor;
    }
  }

  return v;
}
/******************************************************************************/

void pmn_polynomial_values ( int *n_data, int *n, int *m, double *x,
  double *fx )

/******************************************************************************/
/*
  Purpose:

    PMN_POLYNOMIAL_VALUES: normalized Legendre polynomial Pmn(n,m,x).

  Discussion:

    In Mathematica, the unnormalized function can be evaluated by:

      LegendreP [ n, m, x ]

    The function is normalized by dividing by

      sqrt ( 2 * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, int *N, int *M, double *X,
    the arguments of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
    0.7071067811865475E+00, 
    0.6123724356957945E+00, 
   -0.7500000000000000E+00, 
   -0.1976423537605237E+00, 
   -0.8385254915624211E+00, 
    0.7261843774138907E+00, 
   -0.8184875533567997E+00, 
   -0.1753901900050285E+00, 
    0.9606516343087123E+00, 
   -0.6792832849776299E+00, 
   -0.6131941618102092E+00, 
    0.6418623720763665E+00, 
    0.4716705890038619E+00, 
   -0.1018924927466445E+01, 
    0.6239615396237876E+00, 
    0.2107022704608181E+00, 
    0.8256314721961969E+00, 
   -0.3982651281554632E+00, 
   -0.7040399320721435E+00, 
    0.1034723155272289E+01, 
   -0.5667412129155530E+00 };

  static int m_vec[N_MAX] = {
    0, 0, 1, 0,
    1, 2, 0, 1,
    2, 3, 0, 1,
    2, 3, 4, 0,
    1, 2, 3, 4,
    5 };

  static int n_vec[N_MAX] = {
    0,  1,  1,  2,
    2,  2,  3,  3,
    3,  3,  4,  4,
    4,  4,  4,  5,
    5,  5,  5,  5,
    5 };

  static double x_vec[N_MAX] = {
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *m = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *m = m_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

double *pmns_polynomial_value ( int mm, int n, int m, double x[] )

/******************************************************************************/
/*
  Purpose:

    PMNS_POLYNOMIAL_VALUE: sphere-normalized Legendre polynomial Pmn(n,m,x).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

  Parameters:

    Input, int MM, the number of evaluation points.

    Input, int N, the maximum first index of the Legendre
    function, which must be at least 0.

    Input, int M, the second index of the Legendre function,
    which must be at least 0, and no greater than N.

    Input, double X[MM], the evaluation points.

    Output, double PMNS_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
*/
{
  double factor;
  int i;
  int j;
  const double pi = 3.141592653589793;
  double *v;

  v = pm_polynomial_value ( mm, n, m, x );
/*
  Normalization.
*/
  for ( j = m; j <= n; j++ )
  {
    factor = sqrt ( ( ( double ) ( 2 * j + 1 ) * r8_factorial ( j - m ) )
      / ( 4.0 * pi * r8_factorial ( j + m ) ) );
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = v[i+j*mm] * factor;
    }
  }

  return v;
}
/******************************************************************************/

void pmns_polynomial_values ( int *n_data, int *n, int *m, double *x,
  double *fx )

/******************************************************************************/
/*
  Purpose:

    PMNS_POLYNOMIAL_VALUES: sphere-normalized Legendre polynomial Pmns(n,m,x).

  Discussion:

    In Mathematica, the unnormalized function can be evaluated by:

      LegendreP [ n, m, x ]

    The function is normalized for the sphere by dividing by

      sqrt ( 4 * pi * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, int *N, int *M, double *X,
    the arguments of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     0.2820947917738781,
     0.2443012559514600,
    -0.2992067103010745,
    -0.07884789131313000,
    -0.3345232717786446,
     0.2897056515173922,
    -0.3265292910163510,
    -0.06997056236064664,
     0.3832445536624809,
    -0.2709948227475519,
    -0.2446290772414100,
     0.2560660384200185,
     0.1881693403754876,
    -0.4064922341213279,
     0.2489246395003027,
     0.08405804426339821,
     0.3293793022891428,
    -0.1588847984307093,
    -0.2808712959945307,
     0.4127948151484925,
    -0.2260970318780046 };

  static int m_vec[N_MAX] = {
    0, 0, 1, 0,
    1, 2, 0, 1,
    2, 3, 0, 1,
    2, 3, 4, 0,
    1, 2, 3, 4,
    5 };

  static int n_vec[N_MAX] = {
    0,  1,  1,  2,
    2,  2,  3,  3,
    3,  3,  4,  4,
    4,  4,  4,  5,
    5,  5,  5,  5,
    5 };

  static double x_vec[N_MAX] = {
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *m = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *m = m_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

double *pn_pair_product ( int p )

/******************************************************************************/
/*
  Purpose:

    PN_PAIR_PRODUCT: pair products for normalized Legendre polynomial Pn(n,x).

  Discussion:

    Let P(n,x) represent the Legendre polynomial of degree n.

    To check orthonormality, we compute

      Tij = Integral ( -1.0 <= X <= +1.0 ) Pn(i,x) * Pn(j,x) dx

    We will estimate these integrals using Gauss-Legendre quadrature.

    The computed table should be the identity matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int P, the maximum degree of the polyonomial
    factors.  0 <= P.

    Input, int E, the exponent of X in the integrand.
    0 <= E.

    Output, double PN_PAIR_PRODUCT[(P+1)*(P+1)], the table of integrals.
*/
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = ( double * ) malloc ( (p+1)*(p+1) * sizeof ( double ) );;

  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = p + 1;

  x_table = ( double * ) malloc ( order * sizeof ( double ) );;
  w_table = ( double * ) malloc ( order * sizeof ( double ) );;

  p_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = pn_polynomial_value ( 1, p, x_table + k );

    for ( i = 0; i <= p; i++ )
    {
      for ( j = 0; j <= p; j++ )
      {
        table[i+j*(p+1)] = table[i+j*(p+1)] + w_table[k] * h_table[i] * h_table[j];
      }
    }
    free ( h_table );
  }

  free ( w_table );
  free ( x_table );

  return table;
}
/******************************************************************************/

double *pn_polynomial_coefficients ( int n )

/******************************************************************************/
/*
  Purpose:

    PN_POLYNOMIAL_COEFFICIENTS: coefficients of normalized Legendre Pn(n,x).

  Discussion:

    Pn(n,x) = P(n,x) * sqrt ( (2n+1)/2 )

          1       x       x^2     x^3     x^4      x^5    x^6     x^7

    0   0.707
    1   0.000   1.224
    2  -0.790   0.000   2.371
    3   0.000  -2.806   0.000   4.677
    4   0.795   0.000  -7.954   0.000   9.280
    5   0.000   4.397   0.000 -20.520   0.000   18.468
    6  -0.796   0.000  16.731   0.000 -50.193    0.000  36.808
    7   0.000  -5.990   0.000  53.916   0.000 -118.616   0.000  73.429 

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2014

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.

  Parameters:

    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.

    Output, double PN_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients of
    the normalized Legendre polynomials of degree 0 through N.
*/
{
  double *c;
  int i;
  int j;
  double t;

  if ( n < 0 )
  {
    return NULL;
  }
/*
  Compute P(i,x) coefficients.
*/
  c = ( double * ) malloc ( (n+1)*(n+1) * sizeof ( double ) );;

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( 0 < n )
  {
    c[1+1*(n+1)] = 1.0;
  }

  for ( i = 2; i <= n; i++ )
  {
    for ( j = 0; j <= i-2; j++ )
    {
      c[i+j*(n+1)] =
          ( double ) (   - i + 1 ) * c[i-2+j*(n+1)] / ( double ) i;
    }
    for ( j = 1; j <= i; j++ )
    {
      c[i+j*(n+1)] = c[i+j*(n+1)]
        + ( double ) ( i + i - 1 ) * c[i-1+(j-1)*(n+1)] / ( double ) i;
    }
  }
/*
  Normalize them.
*/
  for ( i = 0; i <= n; i++ )
  {
    t = sqrt ( ( double ) ( 2 * i + 1 ) / 2.0 );
    for ( j = 0; j <= i; j++ )
    {
      c[i+j*(n+1)] = c[i+j*(n+1)] * t;
    }
  }

  return c;
}
/******************************************************************************/

double *pn_polynomial_value ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    PN_POLYNOMIAL_POLYNOMIAL evaluates the normalized Legendre polynomials Pn(n,x).

  Discussion:

    The normalized Legendre polynomials are orthonormal under the inner product
    defined as integration from -1 to 1:

      Integral ( -1 <= x <= +1 ) Pn(i,x) * Pn(j,x) dx = delta(i,j)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.

  Parameters:

    Input, int M, the number of evaluation points.

    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.

    Input, double X[M], the evaluation points.

    Output, double PN_POLYNOMIAL_VALUE[M*(N+1)], the values of the Legendre
    polynomials of order 0 through N.
*/
{
  int i;
  int j;
  double norm;
  double *v;

  v = p_polynomial_value ( m, n, x );

  for ( j = 0; j <= n; j++ )
  {
    norm = sqrt ( 2 / ( double ) ( 2 * j + 1 ) );
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = v[i+j*m] / norm;
    }
  }
  return v;
}
/******************************************************************************/

double r8_epsilon ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  const double value = 2.220446049250313E-016;

  return value;
}
/******************************************************************************/

double r8_factorial ( int n )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL computes the factorial of N.

  Discussion:

    factorial ( N ) = product ( 1 <= I <= N ) I

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the argument of the factorial function.
    If N is less than 1, the function value is returned as 1.

    Output, double R8_FACTORIAL, the factorial of N.
*/
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = - 1.0;
  }
  else
  {
    value = + 1.0;
  }
  return value;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

double r8vec_dot_product ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], the two vectors to be considered.

    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
/******************************************************************************/

double *r8vec_linspace_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
 
    In other words, the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the first and last entries.

    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - 1 - i ) * a 
             + ( double ) (         i ) * b ) 
             / ( double ) ( n - 1     );
    }
  }
  return x;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void r8vec2_print ( int n, double a1[], double a2[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_PRINT prints an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A1[N], double A2[N], the vectors to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %4d: %14f  %14f\n", i, a1[i], a2[i] );
  }

  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

