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

#ifndef CVONEDBFSOLVER_H
#define CVONEDBFSOLVER_H

//
//  cvOneDBFSolver.h - Header for a One-Dimensional Network Blood Flow Solver
//  ~~~~~~~~~~~~
//  This is a static class which abstracts the blood flow solver.
//  Essentially, this is a driver for 1D finite element blood flow solver
//  without any interfaces.
//

# include <vector>
# include <iostream>
# include <ostream>
# include <cstdlib>
# include <cstdio>
# include <math.h>

# include "cvOneDTypes.h"
# include "cvOneDModel.h"
# include "cvOneDSubdomain.h"
# include "cvOneDMthModelBase.h"
# include "cvOneDFEAJoint.h"

using namespace std;

class cvOneDBFSolver{

 public:

    // Solver Initialization
    static void SetDeltaTime(double dt);
    static void SetStepSize(long size);
    static void SetInletBCType(BoundCondType bc);
    static void SetMaxStep(long maxs);
    static void SetQuadPoints(long quadPoints_);
    static void SetConvergenceCriteria(double convCriteria);

    // Set the Model Pointer
    static void SetModelPtr(cvOneDModel *mdl);
    static cvOneDModel* GetModelPtr(){return model;}

    // Solve the blood flow problem
    static void Solve(void);

    // Get the solution;
    static double GetSolution(int i, int j){return TotalSolution[i][j];}//IV 082103

    // Cleanup
    static void Cleanup(void);

    static void DefineInletFlow(double* time, double* flrt, int num);

    static void QuerryModelInformation(void);
    static double currentTime;
    static double deltaTime;
    static BoundCondType inletBCtype;
    static bool useFiniteDifferencesTangent;
    static int ASCII;

    // Result Output
    static void postprocess_Text();
    static void postprocess_VTK();
    static void postprocess_VTK_XML3D_ONEFILE();
    static void postprocess_VTK_XML3D_MULTIPLEFILES();

    // Find Segment index given the ID
    static int getSegmentIndex(int segID);

 private:

    // Query Model, Allocate Memory, Set Initial Conditions
    static void CreateGlobalArrays(void);

    //initialize the solution, flow rate and area
    static void CalcInitProps(long subdomainID);
    //the main solve part
    static void GenerateSolution(void);
    //create MthSegmentModel and MthBranchModel if exists. Also specify inflow profile
    static void DefineMthModels(void);
    static void AddOneModel(cvOneDMthModelBase* model);

    static bool wasSet;

    // Pointer to the Model
    static cvOneDModel *model;

    static vector<cvOneDSubdomain*> subdomainList;
    static vector<cvOneDFEAJoint*> jointList;
    static vector<int> outletList;
    static cvOneDFEAVector *currentSolution;
    static cvOneDFEAVector *previousSolution;
    static cvOneDFEAVector *increment;
    static cvOneDFEAVector *rhs;
    static cvOneDFEAVector *relLength; //for pressure calculation
    static cvOneDMatrix<double> TotalSolution;

    // Generic matrix, can be skyline or sparse
    static cvOneDFEAMatrix *lhs;

    static long stepSize;
    static long maxStep;
    static long quadPoints;
    static vector<cvOneDMthModelBase*> mathModels;

    static long numFlowPts;
    static double *flowRate;
    static double *flowTime;
    static double Period;
    static double convCriteria;
};

#endif //CVONEDBFSOLVER_H
