//
//  SkylineMatrix.cxx - Solve Skyline Matrices
//  ~~~~~~~~~~~~~~~~~
//

#include <cassert>
#include <iostream>
#include <math.h>
#include "cvOneDSparseMatrix.h"
#include "cvOneDDenseMatrix.h"
#include <time.h>
using namespace std;

cvOneDSparseMatrix::cvOneDSparseMatrix(const char* tit): cvOneDFEAMatrix(tit){
  numEntries_ = 0;
  numNonzeros_= 0;
  allocatedSizeEntries_ = 0;
  Kentries_ = NULL;
  dim_ = 0;
  allocateKentries();
}

cvOneDSparseMatrix::cvOneDSparseMatrix(long dim, long* pos, const char* tit): cvOneDFEAMatrix(tit){
  numEntries_ = 0;
  allocatedSizeEntries_ = dim*40;
  Kentries_ = NULL;
  dim_ = dim;
  allocateKentries();
}

cvOneDSparseMatrix::~cvOneDSparseMatrix(){
  if (Kentries_ != NULL){
    delete [] Kentries_;
  }
}

void cvOneDSparseMatrix::allocateKentries(){
  int i;
  allocatedSizeEntries_ += ALLOCATE_KENTRIES_INCR;
  cvOneDKentry* ptr = new cvOneDKentry[allocatedSizeEntries_];
  if (numEntries_ > 0) {
    for (i = 0; i < numEntries_; i++) {
      ptr[i].row   = Kentries_[i].row;
      ptr[i].col   = Kentries_[i].col;
      ptr[i].value = Kentries_[i].value;
    }
    delete [] Kentries_;
  }
  Kentries_ = ptr;
}

void cvOneDSparseMatrix::AddValue( long row, long column, double value){
  if (fabs(value)>0.0000001) {
    Kentries_[numEntries_].row = row;
    Kentries_[numEntries_].col = column;
    Kentries_[numEntries_].value = value;
    numEntries_++;
    
    if (numEntries_ == allocatedSizeEntries_) {
      printf("Allocated size reached. Increasing size by 50000.\n");
      allocateKentries();
    }
  }
}

void cvOneDSparseMatrix::ClearRow(long row){
  for (int i = 0; i < numEntries_; i++){
    if (Kentries_[i].row == row) {
      Kentries_[i].value = 0.0;
    }
  }
}

void cvOneDSparseMatrix::ClearColumn(long column){
  for (int i = 0; i < numEntries_; i++){
    if (Kentries_[i].col == column) {
      Kentries_[i].value = 0.0;
    }
  }
}

void cvOneDSparseMatrix::Clear(){
  numEntries_ = 0;
}

void cvOneDSparseMatrix::SetValue(long row, long column, double value){
  for (int i = 0; i < numEntries_; i++){
    if ((Kentries_[i].row == row) && (Kentries_[i].col == column)) {
      Kentries_[i].value = 0.0;
    }
  }
  AddValue(row,column,value);
}

double cvOneDSparseMatrix::GetValue(long row, long column){
  double rtnval = 0.0;
  for (int i = 0; i < numEntries_; i++){
   if ((Kentries_[i].row == row) && (Kentries_[i].col == column)) {
     rtnval += Kentries_[i].value;
    }
  }
  return rtnval;
}

int cvOneDSparseMatrix::GetColumnEntries(long column, cvOneDKentry **Krtr){

  cvOneDKentry *Kcol = NULL;

  int i;
  int numInCol = 0;
  for (i = 0; i < numEntries_; i++){
    if (Kentries_[i].col == column) {
      numInCol++;
    }
  }

  if (numInCol == 0) {
    *Krtr = NULL;
    return 0;
  }

  Kcol = new cvOneDKentry[numInCol];
  int rtnnum = 0;
  for (i = 0; i < numEntries_; i++){
    if (Kentries_[i].col == column) {
      Kcol[rtnnum].row = Kentries_[i].row;
      Kcol[rtnnum].col = Kentries_[i].col;
      Kcol[rtnnum].value = Kentries_[i].value;
      rtnnum++;
    }
  }

  *Krtr = Kcol;

  return rtnnum;
}

void cvOneDSparseMatrix::Add(cvOneDDenseMatrix& matrix){
  int i,j;
  // acquires matrix's data
  long eDimension	= matrix.GetDimension();
  const long* eqNumbers = matrix.GetEquationNumbers();
  double* eEntries = matrix.GetPointerToEntries();

  // assembling matrix's upper and lower entries
  for( i = 0; i < eDimension; i++){
    for( j = 0; j < eDimension; j++){
        AddValue(eqNumbers[i],eqNumbers[j],eEntries[i*eDimension+j]);
    }
  }
}

void cvOneDSparseMatrix::siftDownKentries(cvOneDKentry Kentries[], int root, int bottom){
  int done, maxChild;
  int row, col;
  double value;

  done = 0;
  while((root*2 <= bottom) && (!done)){
    if(root*2 == bottom){
      maxChild = root * 2;
    }else if (Kentries[root * 2].row > Kentries[root * 2 + 1].row){
      maxChild = root * 2;
    }else{
      maxChild = root * 2 + 1;
    }

    if (Kentries[root].row < Kentries[maxChild].row){
      row   = Kentries[root].row;
      col   = Kentries[root].col;
      value = Kentries[root].value;
      Kentries[root].row   = Kentries[maxChild].row;
      Kentries[root].col   = Kentries[maxChild].col;
      Kentries[root].value = Kentries[maxChild].value;
      Kentries[maxChild].row   = row;
      Kentries[maxChild].col   = col;
      Kentries[maxChild].value = value;
      root = maxChild;
    }else{
      done = 1;
    }
  }
}

void cvOneDSparseMatrix::CondenseMatrix(){

  int i,j;

  // need to sort entries here, consolidate them into flat structure
  int row, col;
  double value;

  int array_size = numEntries_;

  clock_t tstart_condense;
  clock_t tstart_sum;
  clock_t tend_sum;
  clock_t tend_condense;


  // heap sort by row
  Build_MaxHeap(Kentries_, array_size-1);
  HeapSort(Kentries_,array_size-1);


  // sum duplicate entries, tag duplicate entries and remove entries with col=-1 in stanfordSparse.cpp
  array_size = numEntries_;

  for (i=0;i<numEntries_;i++){
    if (Kentries_[i].col!=-1) {
       for (j=i+1;j<numEntries_;j++){
         if(Kentries_[i].row==Kentries_[j].row ){
             if(Kentries_[i].col==Kentries_[j].col){
         Kentries_[i].value=Kentries_[i].value+Kentries_[j].value;
         //tage duplicate entry
         Kentries_[j].col=-1;
         Kentries_[j].value=0.0;
         array_size=array_size-1;
             }
          } else {
          break;
          }
       }
      }

   }

  numNonzeros_ = array_size;
  
}



void cvOneDSparseMatrix::MaxHeapify(cvOneDKentry Kentries[], int i, int n)
{
	int j, temp_row,temp_col;
	double temp_val;
	temp_row = Kentries[i].row;
	temp_col = Kentries[i].col;
	temp_val = Kentries[i].value;
	j = 2*i;

 	while (j <= n)
	{
		if (j < n && Kentries[j+1].row > Kentries[j].row)
		j = j+1;
		// Break if parent value is already greater than child value.
		if (temp_row > Kentries[j].row)
			break;
		// Switching value with the parent node if temp < a[j].
		else if (temp_row <= Kentries[j].row)
		{
			Kentries[j/2].row = Kentries[j].row;
			Kentries[j/2].col = Kentries[j].col;
			Kentries[j/2].value = Kentries[j].value;

			j = 2*j;
		}
	}
	Kentries[j/2].row = temp_row;
	Kentries[j/2].col = temp_col;
	Kentries[j/2].value = temp_val;
	return;
}
void cvOneDSparseMatrix::Build_MaxHeap(cvOneDKentry Kentries[], int n)
{
	int i;
	for(i = n/2; i >= 1; i--)
		MaxHeapify(Kentries, i, n);
}
void cvOneDSparseMatrix::HeapSort(cvOneDKentry Kentries[], int n)
{

	int i, temp_row,temp_col;
	double temp_val;
	for (i = n; i >= 2; i--){
		// Storing maximum value at the end.
		temp_row = Kentries[i].row;
		temp_col = Kentries[i].col;
		temp_val = Kentries[i].value;



		Kentries[i].row = Kentries[1].row;
		Kentries[i].col = Kentries[1].col;
		Kentries[i].value = Kentries[1].value;

		Kentries[1].row = temp_row;
		Kentries[1].col =temp_col;
		Kentries[1].value=temp_val;
		// Building max heap of remaining element.
		MaxHeapify(Kentries, 1, i - 1);
	}
}

long cvOneDSparseMatrix::GetDimension() const{
  return dim_;
}

void cvOneDSparseMatrix::print(std::ostream &os){ 

  // Init Full Matrix
  cvDoubleMat fullmatrix;
  fullmatrix.resize(dim_);
  for(long loopA=0;loopA<dim_;loopA++){
    fullmatrix[loopA].resize(dim_);
    for(long loopB=0;loopB<dim_;loopA++){
      fullmatrix[loopA][loopB] = 0.0;
    }
  }

  // Fill Matrix
  for(long loopA=0;loopA<numEntries_;loopA++){
    if(Kentries_[loopA].col > -1){
      fullmatrix[Kentries_[loopA].row][Kentries_[loopA].col] = Kentries_[loopA].value;
    }
  }

  // Print Full Matrix
  for(int i = 0; i < dim_; i++){
    for(int j = 0; j < dim_; j++){
      os << fullmatrix[i][j] << "  ";
    }
    os << "\n";
  }

}


