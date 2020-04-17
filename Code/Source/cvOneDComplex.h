/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CVONEDCOMPLEX_H
#define CVONEDCOMPLEX_H

/***************************************************************************
* Complex.h                                                                *
*                                                                          *
* Complex numbers, complex arithmetic, and functions of a complex          *
* variable.                                                                *
***************************************************************************/

#include <cmath>

class cvOneDComplex{
  public:
    cvOneDComplex(){x = 0;y = 0;}
    cvOneDComplex(double a){x = a;y = 0;}
    cvOneDComplex(double a,double b){x = a;y = b;}
    cvOneDComplex(const cvOneDComplex &z){*this = z;}
    double &Real(){return x;}
    double &Imag(){return y;}
    double Real()const{return x;}
    double Imag()const{return y;}
    inline cvOneDComplex &operator=(const cvOneDComplex &z);
    static const cvOneDComplex i;
  private:
    double x;
    double y;
};

inline cvOneDComplex &cvOneDComplex::operator=(const cvOneDComplex &z){ 
  x = z.Real(); 
  y = z.Imag(); 
  return *this;
}

inline double Real(const cvOneDComplex &z){
  return z.Real();
}

inline double Imag(const cvOneDComplex &z){
  return z.Imag();
}

inline cvOneDComplex conj(const cvOneDComplex &z){
  return cvOneDComplex( z.Real(), -z.Imag() );
}

inline double modsqr(const cvOneDComplex &z){
  return z.Real() * z.Real() + z.Imag() * z.Imag();
}

inline double cmodulus(const cvOneDComplex &z){
  return sqrt( z.Real() * z.Real() + z.Imag() * z.Imag() );
}

inline double arg(const cvOneDComplex &z){
  double t = acos( z.Real() / cmodulus(z) );
  if( z.Imag() < 0.0 ) t = 2.0 * 3.1415926 - t;
  return t;
}

inline cvOneDComplex operator*(const cvOneDComplex &z, double a){
  return cvOneDComplex( a * z.Real(), a * z.Imag());
}

inline cvOneDComplex operator*(double a, const cvOneDComplex &z){
  return cvOneDComplex( a * z.Real(), a * z.Imag() );
}

inline cvOneDComplex operator*(const cvOneDComplex &z, const cvOneDComplex &w){
  return cvOneDComplex(z.Real() * w.Real() - z.Imag() * w.Imag(),z.Real() * w.Imag() + z.Imag() * w.Real());
}

inline cvOneDComplex operator+(const cvOneDComplex &z, const cvOneDComplex &w){
  return cvOneDComplex(z.Real() + w.Real(), z.Imag() + w.Imag());
}

inline cvOneDComplex operator-(const cvOneDComplex &z, const cvOneDComplex &w){
  return cvOneDComplex(z.Real() - w.Real(), z.Imag() - w.Imag());
}

inline cvOneDComplex operator-(const cvOneDComplex &z){
  return cvOneDComplex(-z.Real(), -z.Imag());
}

inline cvOneDComplex operator/(const cvOneDComplex &z, double w){
  return cvOneDComplex(z.Real() / w, z.Imag() / w);
}

inline cvOneDComplex operator/(const cvOneDComplex &z, const cvOneDComplex &w){
  return ( z * conj(w) ) / modsqr(w);
}

inline cvOneDComplex operator/(double a, const cvOneDComplex &w){
  return conj(w) * ( a / modsqr(w) );
}

inline cvOneDComplex &operator+=(cvOneDComplex &z, const cvOneDComplex &w){
  z.Real() += w.Real();
  z.Imag() += w.Imag();
  return z;
}

inline cvOneDComplex &operator*=(cvOneDComplex &z, const cvOneDComplex &w){
  return z = ( z * w );
}

inline cvOneDComplex &operator-=(cvOneDComplex &z, const cvOneDComplex &w){
  z.Real() -= w.Real();
  z.Imag() -= w.Imag();
  return z;
}

inline cvOneDComplex exp(const cvOneDComplex &z){
  double r = exp( z.Real() );
  return cvOneDComplex(r * cos( z.Imag() ), r * sin( z.Imag()));
}

inline cvOneDComplex pow(const cvOneDComplex &z, int n){
  double r = pow( cmodulus( z ), (double)n );
  double t = arg( z );
  return cvOneDComplex( r * cos( n * t ), r * sin( n * t ) );
}

inline cvOneDComplex polar(double r, double theta){
  return cvOneDComplex( r * cos( theta ), r * sin( theta ) );
}

// Returns the square root of this complex number
inline cvOneDComplex sqrt(const cvOneDComplex &z){
  return polar(std::pow(cmodulus(z),0.5),arg(z)*0.5);
}

cvOneDComplex cos ( const cvOneDComplex &z );
cvOneDComplex sin ( const cvOneDComplex &z );
cvOneDComplex cosh( const cvOneDComplex &z );
cvOneDComplex sinh( const cvOneDComplex &z );
cvOneDComplex log ( const cvOneDComplex &z );

#endif // CVONEDCOMPLEX_H

