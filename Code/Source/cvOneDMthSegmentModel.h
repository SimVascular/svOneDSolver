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

    void FormNewton(cvOneDFEAVector* rhsVector, cvOneDFEAMatrix* lhsMatrix);
    void FormNewtonLHS( cvOneDFEAMatrix* lhsMatrix);
    // forms an approximation to the global consistent tangent
    void FormNewtonRHS( cvOneDFEAVector* rhsVector);
    // forms minus the global residual vector
    void SetEquationNumbers( long element, cvOneDDenseMatrix* elementMatrix, int ith);
    long GetUpmostEqnNumber(long ele, long ith) { return -2;}
    // 1=Brooke's one, 0=none IV 04-28-03
    static int STABILIZATION;

  private:

    void FormElement(long element, long ith, cvOneDFEAVector* elementVector, cvOneDDenseMatrix* elementMatrix);
    void FormElementLHS(long element, cvOneDDenseMatrix* elementMatrix, long ithSubdomain);
    void FormElementRHS(long element, cvOneDFEAVector* elementVector, long ithSubdomain);
    void FormMixedBCLHS(int ith, cvOneDSubdomain* sub, cvOneDDenseMatrix* elementMatrix){;}
    void FormMixedBCRHS(int ith, cvOneDSubdomain* sub, cvOneDDenseMatrix* elementMatrix){;}
    double N_Stenosis( long ith);
    void N_MinorLoss(long ith, double* N_vec);
    double GetInflowRate();

  private:

    long quadPoints;
    double* weight;
    double* xi;
    cvOneDQuadrature quadrature_;
};

#endif // CVONEDMTHSEGMENTMODEL_H
