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

#ifndef CVONEDMTHBRANCHMODEL_H
#define CVONEDMTHBRANCHMODEL_H

//
//  cvOneDMthBranchModel.h
//  ~~~~~~~~~~~~~~~~~~
//
//  Mathematical formulations of joints. Lagrange multipliers are used to have
//  pressure continuity and mass conservation (flux conservation) at joints.
//


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
