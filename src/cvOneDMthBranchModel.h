#ifndef CVONEDMTHBRANCHMODEL_H
#define CVONEDMTHBRANCHMODEL_H

//
//  cvOneDMthBranchModel.h
//  ~~~~~~~~~~~~~~~~~~
//
//  Mathematical formulations of joints. Lagrange multipliers are used to have
//  pressure continuity and mass conservation (flux conservation) at joints.
//
//  History:
//  May 1999, J.Wan
//      creation of file,

# include "cvOneDMthModelBase.h"
# include "cvOneDFEAVector.h"
# include "cvOneDFEAMatrix.h"

class cvOneDMthBranchModel: public cvOneDMthModelBase{

  public:

    cvOneDMthBranchModel(const vector<cvOneDSubdomain*> &subdList,
                         const vector<cvOneDFEAJoint*> &jtList, const vector<int>& outletList);
    ~cvOneDMthBranchModel(){}
    int GetNumberOfJoints() {return numOfJoints;}
    void FormNewtonLHS(cvOneDFEAMatrix* lhsMatrix);
    void FormNewtonRHS(cvOneDFEAVector* rhsVector);
    void GetEquationNumbers(long ele, long* eqNumbers, long ithJoint);
    long GetUpmostEqnNumber(long ele, long ithJoint);

  private:

    void FormLagrangeRHSbyP(long ithJoint, cvOneDFEAVector* rhsVector);
    void FormLagrangeRHSbyQ(long ithJoint, cvOneDFEAVector* rhsVector);
    void FormLagrangeLHSbyP(long ithJoint, cvOneDFEAMatrix* lhs);
    void FormLagrangeLHSbyQ(long ithJoint, cvOneDFEAMatrix* lhs);
    int numOfJoints;
};

#endif // CVONEDMTHBRANCHMODEL_H
