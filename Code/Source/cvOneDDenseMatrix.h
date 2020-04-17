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

#ifndef CVONEDDENSEMATRIXX_H
#define CVONEDDENSEMATRIXX_H

//
//  cvOneDDenseMatrix.h - Header for a class to handle dense matrices
//  ~~~~~~~~~~~~~~~
//  
//  SYNOPSIS...This class creates an n x n dense matrix
//             for use in finite element calculations.
//             Typically, this will be used to store local
//             (element) matrices before they are assembled
//             into the global Jacobian.

# include <iostream>
# include <cstdio>

# include "cvOneDSolverDefinitions.h"

using namespace std;

class cvOneDDenseMatrix{
  
  public:

    cvOneDDenseMatrix( long dim, const char* tit = "matrix");
    cvOneDDenseMatrix( long dim, long* eqNumbers, const char* tit = "matrix");
    virtual ~cvOneDDenseMatrix();
    void SetEquationNumbers( long* eqNumbers);
    void Clear();
    long GetDimension() const {return dimension;}
    long* GetEquationNumbers(){return equationNumbers;}
    double* GetPointerToEntries() {return entries;}
    void Set( long row, long column, double value);
    void Add( long row, long column, double value);
    
  private:

    void CreateMatrix( long dim, const char* tit = "matrix");

    long dimension;
    long* equationNumbers;
    double* entries;
    char title[MAX_STRING_SIZE + 1];

};

#endif // CVONEDDENSEMATRIXX_H
