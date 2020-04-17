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

#ifndef CVONEDFEAMATRIX_H
#define CVONEDFEAMATRIX_H

//
//  cvOneDSkylineMatrix.h - Solve Skyline Matrices
//  ~~~~~~~~~~~~~~~~~
//  matrix storage for direct solver, skyline solver
//  

#include <iostream>
#include "cvOneDSolverDefinitions.h"
#include "cvOneDDenseMatrix.h"

using namespace std;

// Matrix is stored in descending sky-line format
class cvOneDFEAMatrix{
  protected:

    char title[MAX_STRING_SIZE + 1];
  
  public:

    cvOneDFEAMatrix(const char* tit);
    virtual ~cvOneDFEAMatrix();

    // VIRTUAL FUNCTIONS
    virtual void Add(cvOneDDenseMatrix& matrix) = 0;
    virtual void AddValue( long row, long column, double value) = 0;
    virtual void SetValue(long row, long column, double value) = 0;
    virtual void Clear() = 0;
    virtual void ClearRow(long row) = 0;
    virtual void ClearColumn(long column) = 0;
    virtual double GetValue(long row, long column) = 0;
    virtual long GetDimension() const = 0;
    // print matrix
    virtual void print(std::ostream &os) = 0;

};

#endif // CVONEDFEAMATRIX_H
