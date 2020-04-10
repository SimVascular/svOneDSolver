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
