#ifndef CVONEDSKYLINELINEARSOLVER_H
#define CVONEDSKYLINELINEARSOLVER_H

//
//  cvOneDLinearSolver.h - Header for a Linear Skyline Matrix Solver
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

# include "cvOneDLinearSolver.h"
# include "cvOneDFEAMatrix.h"
# include "cvOneDFEAVector.h"

class cvOneDSkylineLinearSolver: public cvOneDLinearSolver{

  public:
 
    cvOneDSkylineLinearSolver();
    virtual ~cvOneDSkylineLinearSolver();

    virtual void SetLHS( cvOneDFEAMatrix* matrix);
    virtual void SetRHS( cvOneDFEAVector* vector);
  
    // matrix is overwritten with its LU decomposition
    // solution gets overwritten with the solution of 
    // the linear system of equations
    virtual void Solve( cvOneDFEAVector& solution);
  
    virtual cvOneDFEAMatrix* GetLHS();
    virtual cvOneDFEAVector* GetRHS();
    virtual void SetSolution( long equation, double value);
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
  
  private:

    static int SolNonSymSysSkyLine(double*, double*, double*, 
				                   double*, long*, double*, 
				                   long, int, double);

    static void solvLT(double*, double*, long*, long);
    static void solvUT(double*, double*, double*, double*, 
			           long*, long);

    static double scalv(double*, double*, long);
};

#endif // CVONEDSKYLINELINEARSOLVER_H
