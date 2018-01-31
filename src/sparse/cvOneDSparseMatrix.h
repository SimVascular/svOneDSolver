#ifndef CVONEDSPARSEMATRIX_H
#define CVONEDSPARSEMATRIX_H

//
//  cvOneDSkylineMatrix.h - Solve Skyline Matrices
//  ~~~~~~~~~~~~~~~~~
//  matrix storage for direct solver, skyline solver
//
//  History:
//  May 1999, J.Wan
//      modification for 1D FEM project
//  Mar. 1999, G.R. Feijoo
//      creation of file
//

# include <iostream>

# include "cvOneDTypes.h"
# include "cvOneDSolverDefinitions.h"
# include "cvOneDDenseMatrix.h"
# include "cvOneDFEAMatrix.h"

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
