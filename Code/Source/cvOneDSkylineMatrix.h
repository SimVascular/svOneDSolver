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

#ifndef CVONEDSKYLINEMATRIX_H
#define CVONEDSKYLINEMATRIX_H

//
//  cvOneDSkylineMatrix.h - Solve Skyline Matrices
//  ~~~~~~~~~~~~~~~~~
//  matrix storage for direct solver, skyline solver
//  

#include <iostream>
#include "cvOneDSolverDefinitions.h"
#include "cvOneDDenseMatrix.h"
#include "cvOneDFEAMatrix.h"

using namespace std;

// Matrix is stored in descending sky-line format
class cvOneDSkylineMatrix: public cvOneDFEAMatrix{
  protected:

    bool wasSet;
    long dimension;
    long* position;
    
    double* KU; // upper diagonal part
    double* KD; // diagonal components
    double* KL; // lower diagonal part
  
  public:

    cvOneDSkylineMatrix( const char* tit = "matrix");
    cvOneDSkylineMatrix( long dim, long* pos, const char* tit = "matrix");
    virtual ~cvOneDSkylineMatrix();
    

    long GetPosition( long row, long column) const;
    long* GetPosition();
    double* GetDiagonalEntries();
    double* GetUpperDiagonalEntries();
    double* GetLowerDiagonalEntries();    
    void Set( long dim, long* pos);	
    // sets the position array and allocates space
    //  friend ostream & operator<<( ostream & oFile, SkylineMatrix& matrix);		// output using Matlab sparse format

    // the following are used to set up Dirichlet-type boundary conditions
    // the names are a bit misleading since the number of entries correspond
    // the number of entries in a row or column excluding the diagonal entry
    // the row and column entries I refer to do not include the diagonal entry
    long GetNumberOfEntriesIn( long equation) const;
    void GetRowEntries( long row, long* columns, double* values) const;
    long GetRowEntries( long row, long* columns) const;
    void GetColumnEntries( long column, long* rows, double* values) const;
    long GetColumnEntries( long column, long* rows) const;
    // this is to perform a quick test
    void SetArrays( double* upperEntries, double* lowerEntries, double* diagonalEntries);

    // VIRTUAL FUNCTIONS
    virtual void Add(cvOneDDenseMatrix& matrix);   
    virtual void AddValue( long row, long column, double value);
    virtual void Clear();
    virtual void SetValue(long row, long column, double value);
    // helper functions for applying dirichlet b.c.
    virtual void ClearRow(long row);
    virtual void ClearColumn(long column);    
    virtual double GetValue(long row, long column);
    virtual long GetDimension() const;
    // print matrix
    virtual void print(std::ostream &os);//to replace friend ostream temp fix from Jing IV 081403 

};

#endif // CVONEDSKYLINEMATRIX_H
