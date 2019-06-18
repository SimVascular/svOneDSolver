//
//  LinearSolver.cxx - Source for a Linear Sparse Matrix Solver
//  ~~~~~~~~~~~~~~~~
//
//  This class provides functionality for solving matrix systems
//  presented in the skyline format, and some special manipulations.
//

# include "cvOneDSparseLinearSolver.h"

# ifdef USE_SUPERLU    
  # include "superLUSolve.h"
# endif

# ifdef USE_CSPARSE
  # include "csparseSolve.h"
# endif

#define dmax(a,b) ((a<b) ? (b) : (a))
#include <time.h>

cvOneDSparseLinearSolver::cvOneDSparseLinearSolver(){
}
    
cvOneDSparseLinearSolver::~cvOneDSparseLinearSolver(){
}

void cvOneDSparseLinearSolver::SetLHS(cvOneDFEAMatrix* matrix){
  lhsMatrix = matrix;
}

void cvOneDSparseLinearSolver::SetRHS(cvOneDFEAVector *vector){
  rhsVector = vector;
}

void cvOneDSparseLinearSolver::Solve(cvOneDFEAVector& sol){

  int i;

  // translate from FEAvector to shifted vector for sparse
  int dim = rhsVector->GetDimension();

  clock_t tstart_solve;
  clock_t tstart_LU;
  clock_t tend_LU;
  clock_t tend_solve;

  double *Fglobal = new double[dim];
  double* entries = rhsVector->GetEntries();

  Fglobal[0] = 0.0;
  for (i = 0; i < dim;i++) {
    Fglobal[i] = entries[i];
  }

  double* soln = new double[dim];

  tstart_solve = clock();
  ((cvOneDSparseMatrix*)lhsMatrix)->CondenseMatrix();
  tstart_LU = clock();

# ifdef USE_SUPERLU    
  superLUSolve(((cvOneDSparseMatrix*)lhsMatrix)->GetKentries(), 
               Fglobal, 
               ((cvOneDSparseMatrix*)lhsMatrix)->GetNumberOfEntries(),
               ((cvOneDSparseMatrix*)lhsMatrix)->GetNumberOfNonzeros(), 
               dim, soln);
# endif

# ifdef USE_CSPARSE
  csparseSolve(((cvOneDSparseMatrix*)lhsMatrix)->GetKentries(), 
               Fglobal, 
               ((cvOneDSparseMatrix*)lhsMatrix)->GetNumberOfEntries(),
               ((cvOneDSparseMatrix*)lhsMatrix)->GetNumberOfNonzeros(), 
               dim, soln);
# endif

  tend_LU = clock();

  double* solution = sol.GetEntries();
  for (i = 0; i < dim; i++) {
      solution[i] = soln[i];
  }

  delete [] soln;
}

cvOneDFEAMatrix* cvOneDSparseLinearSolver::GetLHS(){
  return lhsMatrix;
}

cvOneDFEAVector* cvOneDSparseLinearSolver::GetRHS(){
  return rhsVector;
}

void cvOneDSparseLinearSolver::SetSolution(long equation, double value){
  int i;
  cvOneDKentry* columnValues = NULL;

  int numEntries = 0;
  numEntries = ((cvOneDSparseMatrix*)lhsMatrix)->GetColumnEntries( equation, &columnValues);

  lhsMatrix->ClearRow(equation);
  lhsMatrix->ClearColumn(equation);

  lhsMatrix->SetValue(equation, equation, 1.0);
  (*rhsVector)[equation] = value;

  // Subtract values from the right hand side
  for( i = 0; i < numEntries; i++){
    (*rhsVector)[columnValues[i].row] -= value * columnValues[i].value;
  }

  delete [] columnValues;
}

// The values of one dense matrix are to be changed
//                                k_m-1,m^11 | k_m-1,m^12(kr0)    row 1
//                                k_m-1,m^21 | k_m-1,m^22(kr1)    row 2
//  k_m,m-1^11     k_m,m-1^12     k_m,m^11   | k_m,m^12(kr2)      row 3
//  ------------------------------------------
//  k_m,m-1^11(k0) k_m,m-1^12(k1) k_m,m^21(k2) k_m,m^22(k3)      row 4
//    col 1       col 2       col 3        col 4
//
//  At first, multiply row 4 by k_m, then add it to row 3;
//  secondly multiply column 4 by k_m, then add it to column 3;
//  finally the upper 3x3 matrix is kept.

void cvOneDSparseLinearSolver::Minus1dof(long rbEqnNo, double k_m){
  double k[3]; //the array of k1...k3 shown above
  double kr[3]; //the array of kr1...kr3 shown above
  int i;
  for(i = 3; i > 0; i--){
    k[3-i] = lhsMatrix->GetValue(rbEqnNo-i,rbEqnNo);
    kr[3-i] = lhsMatrix->GetValue(rbEqnNo, rbEqnNo-i);
    lhsMatrix->SetValue(rbEqnNo-i,rbEqnNo,0);
    lhsMatrix->SetValue(rbEqnNo, rbEqnNo-i,0);
  }
  
  kr[2] += lhsMatrix->GetValue(rbEqnNo, rbEqnNo-i)*k_m;
  lhsMatrix->SetValue(rbEqnNo,rbEqnNo,1);
  for(i = 3; i > 0; i--){
    lhsMatrix->AddValue(rbEqnNo-1, rbEqnNo-i, k[3-i]*k_m);
    lhsMatrix->AddValue(rbEqnNo-i, rbEqnNo-1, kr[3-i]*k_m);
  }
  (*rhsVector)[rbEqnNo-1] += (*rhsVector)[rbEqnNo]*k_m;
  (*rhsVector)[rbEqnNo] = 0;
}

void cvOneDSparseLinearSolver::DirectAppResistanceBC(long rbEqnNo, double resistance, double dpds, double rhs){
  int i;
  for(i = 3; i > 1; i--){
    lhsMatrix->SetValue(rbEqnNo, rbEqnNo-i,0);
  }
  lhsMatrix->SetValue(rbEqnNo, rbEqnNo-1,dpds);
  lhsMatrix->SetValue(rbEqnNo,rbEqnNo,-resistance);
  (*rhsVector)[rbEqnNo] = rhs;
}

// Assumes 2 nodes/element and 2degrees of freedom/node
void cvOneDSparseLinearSolver::AddFlux(long rbEqnNo, double* OutletLHS11, double* OutletRHS1){

  lhsMatrix->AddValue(rbEqnNo-1, rbEqnNo-1, *OutletLHS11);
  lhsMatrix->AddValue(rbEqnNo-1, rbEqnNo, *(OutletLHS11+1));
  lhsMatrix->AddValue(rbEqnNo, rbEqnNo-1, *(OutletLHS11+2));
  lhsMatrix->AddValue(rbEqnNo, rbEqnNo, *(OutletLHS11+3));
  (*rhsVector)[rbEqnNo-1] += *OutletRHS1;
  (*rhsVector)[rbEqnNo] += *(OutletRHS1+1);
}
