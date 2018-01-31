#include "cvOneDGlobal.h"

// GLOBAL FLAGS
bool cvOneDGlobal::isCreating = false;
bool cvOneDGlobal::isSolving = false;
int  cvOneDGlobal::outputType = 0; // Default Text Output
int  cvOneDGlobal::vtkOutputType = 0; // Default Multiple Files
int  cvOneDGlobal::CONSERVATION_FORM = 0;

// DEBUG MODE
bool cvOneDGlobal::debugMode = false;

// CURRENT MODEL INDEX
long cvOneDGlobal::currentModel = -1;

// VECTOR OF CREATED MODELS
vector<cvOneDModel*> cvOneDGlobal::gModelList;

// GLOBAL MATERIAL MANAGER OBJECT
cvOneDMaterialManager* cvOneDGlobal::gMaterialManager = NULL;

// GLOBAL Mth SEGMENT MODEL
cvOneDMthSegmentModel* cvOneDGlobal::gMthSegmentModel = NULL;

// GLOBAL SOLVER INSTANCE
cvOneDBFSolver* cvOneDGlobal::gBFSolver = NULL;

// GLOBAL VECTOR OF DATATABLES
vector<cvOneDDataTable*> cvOneDGlobal::gDataTables;

// GENERIC SOLVER INSTANCE
cvOneDLinearSolver* cvOneDGlobal::solver = NULL;

