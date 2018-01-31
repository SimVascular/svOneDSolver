#ifndef CVONEDBFSOLVER_H
#define CVONEDBFSOLVER_H

//
//  cvOneDBFSolver.h - Header for a One-Dimensional Network Blood Flow Solver
//  ~~~~~~~~~~~~
//  This is a static class which abstracts the blood flow solver.
//  Essentially, this is a driver for 1D finite element blood flow solver
//  without any interfaces.
//
//  History:
//
//  ? 2002  B.Steele
//      Added impedance boundary conditions
//  Dec., 2000 J.Wan, B.Steele
//      Added resistence boundary condition and more complicated models.
//  May 1999, J.Wan, B.Steele, G.R.Feijoo, S.A.Spicer and S.Strohband
//      Creation of file, class project of ME234C of T.J.R. Hughes and C.Taylor
//
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
    //static double GetSolution(void)
    static double GetSolution(int i, int j){return TotalSolution[i][j];}//IV 082103

    // Cleanup
    static void Cleanup(void);

    static void DefineInletFlow(double* time, double* flrt, int num);

	// was private, had to make public to add Pressure Wave BC
	static void QuerryModelInformation(void);
    static double currentTime;
    static double deltaTime;
	static BoundCondType inletBCtype;
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

    //static double currentTime;
    //static double deltaTime;
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
