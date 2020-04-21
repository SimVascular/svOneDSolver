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

#ifndef CVONEDSPARSEMATRIX_H
#define CVONEDSPARSEMATRIX_H

//
//  cvOneDSkylineMatrix.h - Solve Skyline Matrices
//  ~~~~~~~~~~~~~~~~~
//  matrix storage for direct solver, skyline solver
//

# include <iostream>

# include "../cvOneDTypes.h"
# include "../cvOneDSolverDefinitions.h"
# include "../cvOneDDenseMatrix.h"
# include "../cvOneDFEAMatrix.h"

#define ALLOCATE_KENTRIES_INCR 50000
#define SPARSE_OFFSET 1

using namespace std;

// Matrix is stored in descending sky-line format
class cvOneDSparseMatrix: public cvOneDFEAMatrix{
  protected:

    cvOneDKentry* Kentries_;
    int numEntries_;
    int numNonzeros_;
    int allocatedSizeEntries_;
    int dim_;

  private:
    void allocateKentries();
    void siftDownKentries(cvOneDKentry Kentries[], int root, int bottom);

    void MaxHeapify(cvOneDKentry Kentries[], int i, int n);
    void Build_MaxHeap(cvOneDKentry Kentries[], int n);
    void HeapSort(cvOneDKentry Kentries[], int n);    

  public:
    cvOneDSparseMatrix( const char* tit = "matrix");
    cvOneDSparseMatrix( long dim, long* pos, const char* tit = "matrix");
    virtual ~cvOneDSparseMatrix();

    // query function to apply dirichlet b.c.
    int GetColumnEntries( long column, cvOneDKentry** Kentries);

    // functions to set values
    int GetNumberOfEntries() {return numEntries_;}
    int GetNumberOfNonzeros() {return numNonzeros_;}
    cvOneDKentry* GetKentries() {return Kentries_;}

    void CondenseMatrix();

    // VIRTUAL FUNCTIONS
    virtual void Add(cvOneDDenseMatrix& matrix);
    virtual void AddValue(long row, long column, double value);
    virtual void Clear();
    virtual void SetValue(long row, long column, double value);
    // helper functions for applying dirichlet b.c.
    virtual void ClearRow(long row);
    virtual void ClearColumn(long column);
    // query methods
    virtual double GetValue(long row, long column);
    virtual long GetDimension() const;
    // print matrix
    virtual void print(std::ostream &os);
};

#endif // CVONEDSPARSEMATRIX_H
