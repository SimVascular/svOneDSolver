#include "csparseSolve.h"
#include "cvOneDException.h"

extern "C" {
  #include "csparse.h"
}

int csparseSolve(cvOneDKentry* Kentries,double *b,int numEntries, int NNZ,int nunknown,double u[]){

  // Set Tolerance
  double tol = 1.0e-12;
  // Set Sparse Ordering AMD
  int sparseOrdering = 1;
  // Init Matrix
  // Aux Triplet Matrix
  cs* T = new cs();
  // Sparse Matrix
  cs* A; 

  // Assign Matrix Entries in Triplet Format
  T->nzmax = numEntries;
  T->nz = NNZ;
  T->m = nunknown;
  T->n = nunknown;
  T->p = new int[NNZ];    // Column Indices
  T->i = new int[NNZ];    // Row Indices
  T->x = new double[NNZ]; // Numerical Values
  int j=0;
  for(int i=0;i<numEntries;i++){
    if(Kentries[i].col!=-1){
      T->i[j] = Kentries[i].row;
      T->p[j] = Kentries[i].col;
      T->x[j] = Kentries[i].value;
      j++;
    }
  }
  
  // Convert Matrix in Compressed Solumn Format
  A = cs_triplet(T);

  // Solve system
  int ok = cs_lusol(A,b,sparseOrdering,tol);
  if (ok == 0){
    std::string errorMsg("Error: Cannot Solve Linear System\n");
    throw cvException(errorMsg.c_str());
  }

  // Copy Solution Back
  for(int loopA=0;loopA<nunknown;loopA++){
    u[loopA] = b[loopA];
  }

  // Return Result
  return ok;
}

