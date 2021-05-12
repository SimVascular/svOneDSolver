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
//  SkylineMatrix.cxx - Solve Skyline Matrices
//

# include <cassert>
# include <iostream>
# include "cvOneDSkylineMatrix.h"
# include "cvOneDDenseMatrix.h"

using namespace std;

cvOneDSkylineMatrix::cvOneDSkylineMatrix(const char* tit): cvOneDFEAMatrix(tit){
  wasSet = false;
}

cvOneDSkylineMatrix::cvOneDSkylineMatrix(long dim, long* pos, const char* tit): cvOneDFEAMatrix(tit){
  wasSet = false;  
  Set( dim, pos);
}

cvOneDSkylineMatrix::~cvOneDSkylineMatrix(){
  if( wasSet){
    delete [] KU;
    delete [] KD;
    delete [] KL;
    delete [] position;
  }
}

long cvOneDSkylineMatrix::GetDimension()const{
  assert( wasSet);
  return dimension;
}

// Cleans up the entries, not the position array
void cvOneDSkylineMatrix::Clear(){
  assert( wasSet);

  long i;
  double *ptr1, *ptr2;
  long nEntries = position[dimension];
  for( i = 0, ptr1 = KD; i < dimension; i++, ptr1++)
    *ptr1 = 0.0;

  for(i = 0, ptr1 = KU, ptr2 = KL; i<nEntries; i++, ptr1++, ptr2++){
    *ptr1 = 0.0;
    *ptr2 = 0.0;
  }
}

void cvOneDSkylineMatrix::Set( long dim, long* pos){
  assert( wasSet == false);
  
  dimension = dim;
  position = new long[ dimension + 1];
  assert( position != 0);
  for( long i = 0; i < dimension + 1; i++)
    position[i] = pos[i];
  
  long nEntries = position[dimension];
  KU = new double[nEntries];
  KL = new double[nEntries];
  KD = new double[dimension];
  
  assert( KU != 0 && KL != 0 && KD != 0);
  
  wasSet = true;
}

void cvOneDSkylineMatrix::Add( cvOneDDenseMatrix& matrix){
  
  assert( wasSet);
  
  long i, j, eq;
  
  // acquires matrix's data
  long eDimension  = matrix.GetDimension();
  const long* eqNumbers = matrix.GetEquationNumbers();
  double* eEntries = matrix.GetPointerToEntries();

  // assembling matrix's upper and lower entries
  for(i=0;i<eDimension;i++){
    for(j=i+1;j<eDimension;j++){
      eq = GetPosition( eqNumbers[i], eqNumbers[j]);
      KU[eq] += eEntries[ i * eDimension + j];
      KL[eq] += eEntries[ j * eDimension + i];
    }
  }
  // assembling matrix's diagonal entries
  for(i=0;i<eDimension;i++){
    KD[eqNumbers[ i]] += eEntries[ i * eDimension + i];
  }
}

long cvOneDSkylineMatrix::GetPosition( long row, long column)const{

  assert( wasSet && row >= 0 && column >= 0 && row <= dimension && column <= dimension);

  if( row < column) // upper diagonal entry
    return( position[column + 1] - (column - row));
  else if( row > column) // lower diagonal entry
    return( position[row + 1] - (row - column));
  else // row == column : diagonal entry
    return( row);
}

long* cvOneDSkylineMatrix::GetPosition(){
  assert( wasSet);
  return position;
}

double* cvOneDSkylineMatrix::GetDiagonalEntries(){
  assert( wasSet);
  return KD;
}

double* cvOneDSkylineMatrix::GetUpperDiagonalEntries(){
  assert( wasSet);
  return KU;
}

double* cvOneDSkylineMatrix::GetLowerDiagonalEntries(){
  assert( wasSet);
  return KL;
}


void cvOneDSkylineMatrix::SetArrays( double* upperEntries, double* lowerEntries, double* diagonalEntries){
  assert( wasSet);
  
  long i;
  long nEntries = position[dimension];
  
  for( i = 0; i < nEntries; i++){
    KU[i] = upperEntries[i];
    KL[i] = lowerEntries[i];
  }

  for( i = 0; i < dimension; i++){
    KD[i] = diagonalEntries[i];
  }
}

long cvOneDSkylineMatrix::GetNumberOfEntriesIn(long equation)const{
  
  // Get the column height
  long numEntries = position[equation + 1] - position[equation];
  
  for( long i = equation + 1; i < dimension; i++){
    // column height is long enough
    if( i - equation <= position[i+1] - position[i]){
      numEntries++; // account for this one...
    }  
  }  
  return numEntries;
}

long cvOneDSkylineMatrix::GetRowEntries(long row, long* columns)const{
  long i;
  
  long numEntries = 0;
  long* ptr_c = columns;
  
  long height = position[row+1] - position[row];
  long col = row - height;
  while( height--){
    *ptr_c++ = col++;
    numEntries++;
  }

  // Traverse the remaining columns
  for( i = row + 1; i < dimension; i++){
    // column height is long enough
    if( i - row <= position[i+1] - position[i]){
      *ptr_c++ = i;
      numEntries++;
    }
  } 
  return numEntries;
}

long cvOneDSkylineMatrix::GetColumnEntries(long column, long* rows)const{
  
  long i;
  
  long numEntries = 0;
  long* ptr_r = rows;

  long height = position[column+1] - position[column];
  long row = column - height;
  while( height--){
    *ptr_r++ = row++;
    numEntries++;
  }
  
  // Traverse the remaining rows
  for(i=column+1;i<dimension;i++){
    // column height is long enough
    if(i - column <= position[i+1] - position[i]){
      *ptr_r++ = i;
      numEntries++;
    }
  }      
  return numEntries;
}

void cvOneDSkylineMatrix::GetRowEntries(long row, long* columns, double* values)const{

  long numEntries = GetRowEntries( row, columns);
  
  // now go for the values
  long p;
  for(long i=0;i<numEntries;i++){
    p = GetPosition( row, columns[i]);
    // lower diagonal part of matrix
    if( row > columns[i]){
      values[i] = KL[p];    
    }else{          
      // otherwise it has to belong to the upper diagonal part of the matrix
      values[i] = KU[p];
    }
  }
}

void cvOneDSkylineMatrix::GetColumnEntries(long column, long* rows, double* values) const{

  long numEntries = GetColumnEntries( column, rows);
  
  // now get the values
  long p;
  for(long i=0;i<numEntries;i++){
    p = GetPosition( rows[i], column);
    // lower diagonal part of matrix
    if( rows[i] > column){
      values[i] = KL[p];    
    }else{
      // obviously, it has to belong to the upper diagonal part of the matrix
      values[i] = KU[p];
    }
  }
}

void cvOneDSkylineMatrix::ClearRow(long row){
  // this method does not modify the column height and therefore the position array
  long numEntries = GetNumberOfEntriesIn(row);
  long* columns = new long[numEntries];
  assert( columns != 0);
  
  GetRowEntries( row, columns);

  for( long i = 0; i < numEntries; i++){
    SetValue( row, columns[i], 0.0);
  }
  delete [] columns;
}

void cvOneDSkylineMatrix::ClearColumn(long column){
  // this method does not modify the column height and therefore the position array
  long numEntries = GetNumberOfEntriesIn( column);
  long* rows = new long[numEntries];
  assert( rows != 0);
  
  GetColumnEntries( column, rows);
  
  for( long i = 0; i < numEntries; i++){
    SetValue( rows[i], column, 0.0);
  }
  delete [] rows;
}

void cvOneDSkylineMatrix::SetValue(long row, long column, double value){
  long p = GetPosition( row, column);
  if( row > column){
    KL[p] = value;
  }else if( row < column){
    KU[p] = value;
  }else{
    KD[row] = value;
  }
}

double cvOneDSkylineMatrix::GetValue(long row, long column){
  long p = GetPosition( row, column);
  double res = 0.0;
  if( row > column){
    res = KL[p];
  }else if( row < column){
    res = KU[p];
  }else{
    res = KD[row];
  }
  return res;
}


void cvOneDSkylineMatrix::AddValue(long row, long column, double value){
  long p = GetPosition(row, column);
  if(row > column){
    KL[p] += value;
  }else if( row < column){
    KU[p] += value;
  }else{
    KD[row] += value;
  }
}

void cvOneDSkylineMatrix::print(std::ostream &os){ 
  
  assert(this->wasSet);

  long dimension = this->dimension;
  char* title = this->title;
  long* position = this->position;
  double* KU = this->KU;
  double* KL = this->KL;
  double* KD = this->KD;
  long i, j; 
  double** fullmatrix = new double*[dimension];
  for(i = 0; i < dimension; i++){
    fullmatrix[i] = new double[dimension];
    for(j = 0; j < dimension; j++){
      fullmatrix[i][j] = 0;
    }
  }
  int * rowid = new int[dimension*dimension];
  int * colid = new int[dimension*dimension];
  double * value = new double[dimension*dimension];
  long height, row, column, pos;
  
  
  int temp = 0;
  // first index corresponds to the row positions
  for( i = 0; i < dimension; i++){ // diagonal entries
    rowid[temp++] = i+1;
  }
  for( i = 0; i < dimension; i++){ // upper diagonal entries
    height = position[i + 1] - position[i]; // number of elements in the skyline
    row = i + 1 - height;
//    printf("Height: %ld, row: %ld\n",height,row);
    while( height--){
      rowid[temp++] = row++;
    }
  }
  for( i = 0; i < dimension; i++){ // lower diagonal entries
    height = position[i + 1] - position[i];
    row = i + 1;
    while( height--){
      rowid[temp++] = row;
    }
  }
  temp = 0;
  // second index corresponds to the column positi
  os << "j" << title << " = [\n";
  for( i = 0; i < dimension; i++){ // diagonal entries
    colid[temp++] = i+1;
  }
  for( i = 0; i < dimension; i++){ // upper diagonal entries
    height = position[i + 1] - position[i];
    column = i + 1;
    while( height--){
      colid[temp++] = column;
    }
  }
  for( i = 0; i < dimension; i++){ // lower diagonal entries
    height = position[i + 1] - position[i];
    column = i + 1 - height;
    while( height--){
      colid[temp++] = column++;
    }
  }
  os << "]; \n";
  temp = 0;
  // now we output the entries of the matrix
  os << "s" << title << " = [ \n" ;
  for( i = 0; i < dimension; i++){ // diagonal entries
    value[temp++] = KD[i];
  }
  for( i = 0; i < dimension; i++){ // upper diagonal entries
    height = position[i + 1] - position[i];
    pos = position[i];
    while( height--){
      value[temp++] = KU[pos++];
    }
  }
  for( i = 0; i < dimension; i++){ // lower diagonal entries   
    height = position[i + 1] - position[i];
    pos = position[i];
    while( height--){
      value[temp++] = KL[pos++];
    }
  }
  // Create Full Matrix
  for(i = 0; i < temp; i++){
//    printf("i: %ld, temp: %d\n",i,temp);
//    printf("rowid-1: %d, colid-1: %d\n",rowid[i]-1,colid[i]-1);
    fullmatrix[rowid[i]-1][colid[i]-1] = value[i];
  }
  // Print Full Matrix
  for(i = 0; i < dimension; i++){
    for(int j = 0; j < dimension; j++){
      os << fullmatrix[i][j] << "  ";
    }
    os << "\n";
  }
  // create sparse matrix using the above arrays
  os << title << " = sparse( i" << title << ", j" << title << ", s" << title << "); \n";
  for(i = 0; i < dimension; i++){
    delete [] fullmatrix[i];
  }
  delete [] fullmatrix;
  delete [] rowid;
  delete [] colid;
  delete [] value;
}
