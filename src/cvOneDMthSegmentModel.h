#ifndef CVONEDMTHSEGMENTMODEL_H
#define CVONEDMTHSEGMENTMODEL_H

//
//  cvOneDMthSegmentModel.h
//  ~~~~~~~~~~~~
//  This class is a derived class of mthModelBase. It handles the mathematical formulations
//  of segment model. Essentially it creates the stiffness matrix and
//  right-hand side for each element
//


# include <iostream>

# include "cvOneDMthModelBase.h"
# include "cvOneDUtility.h"

class cvOneDMthSegmentModel : public cvOneDMthModelBase{

  public:

    cvOneDMthSegmentModel(const vector<cvOneDSubdomain*> &subdList,
                          const vector<cvOneDFEAJoint*> &jtList, 
                          const vector<int> &outletList,
                          long quadPoints_);
    ~cvOneDMthSegmentModel();

    void FormNewtonLHS( cvOneDFEAMatrix* lhsMatrix);
    // forms an approximation to the global consistent tangent
    void FormNewtonRHS( cvOneDFEAVector* rhsVector);
    // forms minus the global residual vector
    void SetEquationNumbers( long element, cvOneDDenseMatrix* elementMatrix, int ith);
    long GetUpmostEqnNumber(long ele, long ith) { return -2;}
    // 1=Brooke's one, 0=none IV 04-28-03
    static int STABILIZATION;

  private:

    void FormElementLHS(long element, cvOneDDenseMatrix* elementMatrix, long ithSubdomain);
    void FormElementRHS(long element, cvOneDFEAVector* elementVector, long ithSubdomain);
    void FormMixedBCLHS(int ith, cvOneDSubdomain* sub, cvOneDDenseMatrix* elementMatrix){;}
    void FormMixedBCRHS(int ith, cvOneDSubdomain* sub, cvOneDDenseMatrix* elementMatrix){;}
    double N_Stenosis( long ith);
    double N_MinorLoss(long ith);
    double GetInflowRate();

  private:

    long quadPoints;
    double* weight;
    double* xi;
    cvOneDQuadrature quadrature_;
};

#endif // CVONEDMTHSEGMENTMODEL_H
