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

