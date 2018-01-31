#ifndef CVONEDFEAMATRIX_H
#define CVONEDFEAMATRIX_H

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
