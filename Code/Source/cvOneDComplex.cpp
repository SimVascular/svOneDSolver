/***************************************************************************
* Complex.C                                                                *
*                                                                          *
* Complex numbers, complex arithmetic, and functions of a complex          *
* variable.                                                                *
*                                                                          *
***************************************************************************/
#include "cvOneDComplex.h"

const cvOneDComplex cvOneDComplex::i(0.0,1.0);

cvOneDComplex cos(const cvOneDComplex &z){
  return cvOneDComplex(cos(z.Real()) * cosh(z.Imag()), -sin(z.Real()) * sinh(z.Imag()));
}

cvOneDComplex sin(const cvOneDComplex &z){
  return cvOneDComplex(sin( z.Real() ) * cosh( z.Imag() ), cos( z.Real() ) * sinh( z.Imag()));
}

cvOneDComplex cosh(const cvOneDComplex &z){
  return cvOneDComplex(cosh( z.Real() ) * cos( z.Imag() ),sinh( z.Real() ) * sin( z.Imag()));
}

cvOneDComplex sinh(const cvOneDComplex &z){
  return cvOneDComplex(sinh( z.Real() ) * cos( z.Imag() ),cosh( z.Real() ) * sin( z.Imag()));
}

cvOneDComplex log(const cvOneDComplex &z){
  double r = sqrt( z.Real() * z.Real() + z.Imag() * z.Imag() );
  double t = acos( z.Real() / r );
  if( z.Imag() < 0.0 ) t = 2.0 * 3.1415926 - t;
    return cvOneDComplex( log(r), t );
}

