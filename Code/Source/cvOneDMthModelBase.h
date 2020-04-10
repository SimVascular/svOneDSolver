#ifndef CVONEDMTHMODELBASE_H
#define CVONEDMTHMODELBASE_H

//
//  cvOneDMthModelBase.h - Header for a class to handle finite element equations
//  ~~~~~~~~~~~~
//
//  This class is the base class in a finite element
//  computation.  It handles the equation formulation, and
//  construction of the Newton-Rhapson matrix derivatives,
//  and generation of the element dense matrices.
//

# include <vector>

# include "cvOneDModel.h"
# include "cvOneDFEAVector.h"
# include "cvOneDFEAMatrix.h"
# include "cvOneDSubdomain.h"
# include "cvOneDFEAJoint.h"

class cvOneDMthModelBase{

  public:

    static int impedIncr;

    cvOneDMthModelBase(const cvOneDModel* modl);
    cvOneDMthModelBase(const vector<cvOneDSubdomain*>& subdList, const vector<cvOneDFEAJoint*>& jtList,
                       const vector<int>& outletList);
    virtual ~cvOneDMthModelBase();
    virtual long GetTotalNumberOfEquations() const {return numberOfEquations;}
    virtual int  GetNumberOfElementEquations() const {return 4;}
    virtual void TimeUpdate(double pTime, double deltaT);
    virtual void FormNewtonLHS(cvOneDFEAMatrix* lhsMatrix) = 0;
    // forms an approximation to the global consistent tangent
    virtual void FormNewtonRHS(cvOneDFEAVector* rhsVector) = 0;
    // forms minus the global residual vector
    virtual void SetBoundaryConditions();
    virtual double CheckMassBalance();
    virtual void ApplyBoundaryConditions();
    virtual void GetNodalEquationNumbers( long node, long* eqNumbers, long ith);
    virtual void GetEquationNumbers( long element, long* eqNumbers, long ith);
    virtual long GetUpmostEqnNumber(long ele, long ith) =0;
    virtual void EquationInitialize(const cvOneDFEAVector* pSolution, cvOneDFEAVector* cSolution){prevSolution = pSolution; currSolution = cSolution;}
    virtual void SetInflowRate(double *t, double *flow, int size, double cycleT);
    typeOfEquation GetType() const {return type;}
    double GetCycleTime() const {return cycleTime;}

  protected:

    double GetFlowRate();

    typeOfEquation type;
    long numberOfEquations;

    // If all these lists are set up just once and won't be changed in the future,
    // they should be normal array, but I keep them as dynamic arrays because
    // if user wants to input more segments, the program can be restarted
    // without relocate everything, the only thing to be done is to append
    // more segments in the Math class.
    vector<cvOneDSubdomain*> subdomainList;
    vector<cvOneDFEAJoint*> jointList;
    vector<int> outletList;
    // Placement of equations in the globalsystem
    long* equationNumbers;

    // Solution at time t_{n}
    const cvOneDFEAVector* prevSolution;

    // Current approximation to the solution at t_{n+1}
    cvOneDFEAVector* currSolution;

    double previousTime;    // t_{n}
    double deltaTime;        // deltat
    double currentTime;    // t_{n+1}
    double *flrt, *time;
    double cycleTime;
    int  nFlowPts; // added by bnsteel

};

#endif // CVONEDMTHMODELBASE_H
