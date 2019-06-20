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
