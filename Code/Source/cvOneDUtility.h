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

#ifndef CVONEDUTILITY_H
#define CVONEDUTILITY_H

//
//  cvOneDUtility.h - Header for Some Utility functions and rules of random point quadratures.
//  ~~~~~~~~~~~   
//  

# include <cmath>
# include <iostream>
# include <cstdarg>
# include <cassert>
# include <string>

# include "cvOneDTypes.h"
# include "cvOneDException.h"

const int MaxChar = 128;

std::vector<std::string> split_string(std::string& s, const std::string& delims);
std::string trim_string(const std::string& s);
std::string upper_string(const std::string& s);

long min(long a, long b);
long max( long a, long b);	
long min( long size, long* values);
long sum( long size, long* values);	
void clear( long size, long* vec);	
// Calculates the modulus of the 2x2 matrix A and put the results
// in modulusA using Cayley-Hamilton theory
void GetModulus(double* A, double* modulusA);

class cvOneDQuadrature{
  public:
    cvOneDQuadrature( int n);
    ~cvOneDQuadrature();
    void Get( double* w, double* x) const;

  private:	
  	// number of integration points in the quadrature
    int npts;  
    double* xi;
    double* weights;
};

int getListIDWithStringKey(string key,cvStringVec list);

#endif // CVONEDUTILITY_H
