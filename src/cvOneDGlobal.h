#ifndef CVONEDGLOBAL_H
#define CVONEDGLOBAL_H

# include <vector>

# include "cvOneDModel.h"
# include "cvOneDMaterialManager.h"
# include "cvOneDMthSegmentModel.h"
# include "cvOneDDataTable.h"
# include "cvOneDLinearSolver.h"
# include "cvOneDBFSolver.h"

using namespace std;

class cvOneDGlobal{
  
  public:

    // GLOBAL FLAGS
    static bool isCreating;
    static bool isSolving;
    static int outputType;
    static int vtkOutputType;
    static int CONSERVATION_FORM;

    // DEBUG MODE
    static bool debugMode;

    // CURRENT MODEL INDEX
    static long currentModel;

    // VECTOR OF CREATED MODELS
    static vector<cvOneDModel*> gModelList;

    // GLOBAL MATERIAL MANAGER OBJECT
    static cvOneDMaterialManager* gMaterialManager;

    // GLOBAL Mth SEGMENT MODEL
    static cvOneDMthSegmentModel* gMthSegmentModel;

    // GLOBAL SOLVER INSTANCE
    static cvOneDBFSolver* gBFSolver;

    // GLOBAL VECTOR OF DATATABLES
    static vector<cvOneDDataTable*> gDataTables;

    // Generic Solver instance
    static cvOneDLinearSolver *solver;

};

#endif // CVONEDGLOBAL_H
