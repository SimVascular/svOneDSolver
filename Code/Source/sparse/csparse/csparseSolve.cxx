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

#include "csparseSolve.h"
#include "../../cvOneDException.h"

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

