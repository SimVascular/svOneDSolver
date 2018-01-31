//
//  DenseMatrix.cxx - Source for a class to handle dense matrices
//  ~~~~~~~~~~~~~~~
//  
//  This class creates an n x n dense matrix
//  for use in finite element calculations.
//  Typically, this will be used to store local
//  (element) matrices before they are assembled
//  into the global Jacobian.
//  
//  May 1999, J.Wan, G.R.Feijoo, S.A.Spicer and S.Strohband
//  Creation of file, class project of ME234C of T.J.R. Hughes and C.Taylor

#include <cassert>
#include <iomanip>

#include "cvOneDDenseMatrix.h"

void cvOneDDenseMatrix::CreateMatrix(long dim, const char* tit){
  long i;
  
  assert( dim > 0);
  dimension = dim;
  
  equationNumbers = new long[dimension];
  entries = new double[dimension * dimension];
  assert( equationNumbers != 0 && entries != 0);
  
  i = 0;
  while(i < MAX_STRING_SIZE && tit[i] != '\0'){
    title[i] = tit[i];
    i++;
  }
  title[i] = '\0';
}

cvOneDDenseMatrix::cvOneDDenseMatrix(long dim, const char* tit){
  CreateMatrix(dim,tit);
}

cvOneDDenseMatrix::cvOneDDenseMatrix(long dim, long* eqNumbers, const char* tit){
  CreateMatrix(dim,tit);
  SetEquationNumbers(eqNumbers);
}

void cvOneDDenseMatrix::SetEquationNumbers(long* eqNumbers){
  for(long i = 0; i < dimension; i++){
    equationNumbers[i] = eqNumbers[i];
  }
}

cvOneDDenseMatrix::~cvOneDDenseMatrix(){
  delete [] equationNumbers;
  delete [] entries;
}

void cvOneDDenseMatrix::Clear(){
  double* ptr = entries;
  for( long i = 0; i < dimension*dimension; i++){
    *ptr++ = 0.0;  
  }
}

void cvOneDDenseMatrix::Set(long row, long column, double value){
  assert( row >= 0 && column < dimension);
  entries[ row * dimension + column] = value;
}

void cvOneDDenseMatrix::Add(long row, long column, double value){
  assert( row >= 0 && column < dimension);
  entries[ row * dimension + column] += value;
}
