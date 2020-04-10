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
//  FEAVector.cxx - Source for a class to handle finite element Vectors
//  ~~~~~~~~~~~~~
//  SYNOPSIS...This is a utility class used to handle the finite element
//             Vectors (currentSolution, increment, previousSolution, etc.
//

# include <cassert>
# include <iomanip>
# include <cstring>
# include <cmath>
# include "cvOneDFEAVector.h"


using namespace std;

void cvOneDFEAVector::CreateVector(long dim, const char* tit){
  assert( dim > 0);

  dimension = dim;
  equationNumbers = new long[dimension];
  entries = new double[dimension];
  assert( equationNumbers != 0 && entries != 0);

  for(long i=0;i<dimension;i++){
    equationNumbers[i] = i;
  }

  Rename(tit);
  Clear();
}

cvOneDFEAVector::cvOneDFEAVector(long dim, const char* tit){
  CreateVector(dim,tit);
}

cvOneDFEAVector::cvOneDFEAVector(long dim, long* eqNumbers, const char* tit){
  CreateVector(dim,tit);
  SetEquationNumbers(eqNumbers);
}

cvOneDFEAVector::~cvOneDFEAVector(){
  delete [] equationNumbers;
  delete [] entries;
}

void cvOneDFEAVector::Clear(){
  for(long i=0;i<dimension;i++){
    entries[i] = 0.0;
  }
}

void cvOneDFEAVector::SetEquationNumbers(long* eqNumbers){
  for(long i=0;i<dimension;i++){
    equationNumbers[i] = eqNumbers[i];
  }
}

void cvOneDFEAVector::Add(cvOneDFEAVector& vector){
  long eDimension = vector.GetDimension();
  const long* eEquationNumbers = vector.GetEquationNumbers();

  long pos;
  for(long i=0;i<eDimension;i++){
    pos = eEquationNumbers[i];
    assert( pos >= 0 && pos < dimension);
    entries[pos] += vector[i];
  }
}

void cvOneDFEAVector::Add(long component, double value){
  assert(component >= 0 && component < dimension);
  entries[component] += value;
}

double& cvOneDFEAVector::operator[](long i){
  assert(i>=0 && i<dimension);
  return(entries[i]);
}

void cvOneDFEAVector::Set(long i, double value){
  assert(i>=0 && i<dimension);
  entries[i] = value;
}

cvOneDFEAVector& cvOneDFEAVector::operator =(const cvOneDFEAVector& rhs){
  assert(dimension == rhs.dimension);

  strcpy(title, rhs.title);
  for(long i=0;i<dimension;i++){
    entries[i] = rhs.entries[i];
    equationNumbers[i] = rhs.equationNumbers[i];
  }
  return *this;
}

cvOneDFEAVector& cvOneDFEAVector::operator +=(const cvOneDFEAVector& rhs){
  assert( dimension == rhs.dimension);

  for(long i = 0; i < dimension; i++){
    entries[i] += rhs.entries[i];
  }
  return *this;
}

double cvOneDFEAVector::Get( long i) const{
  assert(i>=0 && i<dimension);
  return(entries[i]);
}

void cvOneDFEAVector::Rename(const char* tit){
  long i=0;
  while(i < MAX_STRING_SIZE && tit[i] != '\0'){
    title[i] = tit[i];
    i++;
  }
  title[i] = '\0';
}

double cvOneDFEAVector::Norm( normType type, int start, int step, int stop_index)const{
  assert( type == L2_norm); //only L2_norm
  double norm = 0.0;
  double* ptr = entries + start;
  int dim = (stop_index == -1) ? dimension : stop_index;
  for( long i = 0; i < dim; i+=step, ptr+=step)
    norm += (*ptr) * (*ptr);
  return sqrt( norm);
}

void cvOneDFEAVector::CheckPositive(int start, int step, int stop){
  assert(start >= 0 && stop <= dimension);
  for(long i= start; i< stop; i+=step){
    assert(cvOneDFEAVector::Get(i) >= 0.0);
    }
}

void cvOneDFEAVector::print(std::ostream &os){
  os.setf(ios::scientific);
  os.precision( OUTPUT_PRECISION);

  os << this->title << " = [" << endl;
  for( long i = 0; i < this->GetDimension(); i++){
    os << " " << this->Get(i) << ";" << endl;
  }
  os << "];" << endl;
}


