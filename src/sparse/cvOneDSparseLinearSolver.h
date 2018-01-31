#ifndef CVONEDSPARSELINEARSOLVER_H
#define CVONEDSPARSELINEARSOLVER_H

//
//  cvOneDSparseLinearSolver.h - Header for a Linear Sparse Matrix Solver
//  ~~~~~~~~~~~~~~~~
//
//  This class provides functionality for solving matrix systems
//  presented in the skyline format, and some special manipulations.
//
//  History:
//  May 1999, J.Wan
//      modifications for 1D FEM project and special treatments for boundary conditions and joint handling
//  Mar. 1999, G.R. Feijoo
//      creation of file

# include <cmath>
# include <cassert>

# include "cvOneDLinearSolver.h"
# include "cvOneDSparseMatrix.h"
# include "cvOneDFEAMatrix.h"
# include "cvOneDFEAVector.h"

class cvOneDSparseLinearSolver: public cvOneDLinearSolver{

  public:

    cvOneDSparseLinearSolver();
    virtual ~cvOneDSparseLinearSolver();

    virtual void SetLHS( cvOneDFEAMatrix* matrix);
    virtual void SetRHS( cvOneDFEAVector* vector);

    // matrix is overwritten with its LU decomposition
    // solution gets overwritten with the solution of
    // the linear system of equations
    virtual void Solve( cvOneDFEAVector& solution);

    virtual cvOneDFEAMatrix* GetLHS();
    virtual cvOneDFEAVector* GetRHS();
    virtual void SetSolution(long equation, double value);
    // when one more constraint (dQ = k_m*dS, resistance boundary
    // condition) is added, the basic dense matrix (4x4) is decreased
    // to (3x3)
    virtual void Minus1dof( long rightBottomEquationNumber, double k_m);
    // direct application of the resistance constraint without reduction, which means
    // Newton Raphson scheme is applied on equation Q = PR. Therefore the dense matrix
    // is still 4x4.
    virtual void DirectAppResistanceBC(long rbEqnNo, double resistance, double dpds, double rhs);
	virtual void AddFlux(long rbEqnNo, double* OutletLHS11, double* OutletRHS1);
	//AddFlux added by IV 03-26-03, assumes 2 nodes/element and 2degrees of freedom/node
};

#endif // CVONEDSPARSELINEARSOLVER_H
