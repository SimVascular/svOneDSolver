#ifndef CVONEDLINEARSOLVER_H
#define CVONEDLINEARSOLVER_H

//
//  cvOneDLinearSolver.h - Header for a Linear Skyline Matrix Solver
//  
//  This class provides functionality for solving matrix systems
//  presented in the skyline format, and some special manipulations. 
//  

# include <cmath>

# include "cvOneDFEAMatrix.h"
# include "cvOneDFEAVector.h"

class cvOneDLinearSolver{

  public:

    static cvOneDFEAMatrix* lhsMatrix;
    static cvOneDFEAVector* rhsVector;

    cvOneDLinearSolver();
    virtual ~cvOneDLinearSolver();

    virtual void SetLHS( cvOneDFEAMatrix* matrix) = 0;
    virtual void SetRHS( cvOneDFEAVector* vector) = 0;
  
    // matrix is overwritten with its LU decomposition
    // solution gets overwritten with the solution of 
    // the linear system of equations
    virtual void Solve(cvOneDFEAVector& solution) = 0;
  
    virtual cvOneDFEAMatrix* GetLHS() = 0;
    virtual cvOneDFEAVector* GetRHS() = 0;
    virtual void SetSolution(long equation, double value) = 0;
    // when one more constraint (dQ = k_m*dS, resistance boundary
    // condition) is added, the basic dense matrix (4x4) is decreased
    // to (3x3)
    virtual void Minus1dof(long rightBottomEquationNumber, double k_m) = 0;
    // direct application of the resistance constraint without reduction, which means 
    // Newton Raphson scheme is applied on equation Q = PR. Therefore the dense matrix 
    // is still 4x4. 
    virtual void DirectAppResistanceBC(long rbEqnNo, double resistance, double dpds, double rhs) = 0;
	virtual void AddFlux(long rbEqnNo, double* OutletLHS11, double* OutletRHS1) = 0;
	
};

#endif // CVONEDLINEARSOLVER_H
