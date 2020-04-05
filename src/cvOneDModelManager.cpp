# include "cvOneDModelManager.h"

cvOneDModelManager::cvOneDModelManager(char *mdlName){
  // We're creating a model
  cvOneDGlobal::isCreating = true;

  cvOneDModel* newModel = new cvOneDModel;
  newModel->setModelName(mdlName);
  cvOneDGlobal::currentModel = cvOneDGlobal::gModelList.size();
  newModel->setModelID(cvOneDGlobal::currentModel);
  cvOneDGlobal::gModelList.push_back(newModel);
}

cvOneDModelManager::~cvOneDModelManager(){
}

int cvOneDModelManager::CreateMaterial(char *matName, char *MaterialTypeString,
                                       double density, double dynamicViscosity,
                                       double profile_exponent, double pRef,
                                       int numParams, double *params, int *matID){

  if (cvOneDGlobal::gMaterialManager == NULL){
    cvOneDGlobal::gMaterialManager = new cvOneDMaterialManager();
  }
  if(!strcmp (MaterialTypeString, "MATERIAL_OLUFSEN")){
    *matID = cvOneDGlobal::gMaterialManager->AddNewMaterialOlufsen(density,dynamicViscosity,
                             profile_exponent,pRef,params);
    return CV_OK;
  }else if(!strcmp (MaterialTypeString, "MATERIAL_LINEAR")){
    double EHR = params[0];
    *matID = cvOneDGlobal::gMaterialManager->AddNewMaterialLinear(density,dynamicViscosity,
                                                                  profile_exponent,pRef,EHR);
    return CV_OK;
  }else{
    return CV_ERROR;
  }
}

int cvOneDModelManager::CreateSegment(char   *segName,long segID, double  segLen,
                                      long    numEls,long    inNode,long    outNode,
                                      double  InitialInletArea,double  InitialOutletArea,
                                      double  InitialFlow,int matID,char* lossType,
                                      double branchAngle,int upstreamSegment,int branchSegment,
                                      char* boundType,double* value, double* time, int num){

  MinorLossScope::MinorLoss loss;
  BoundCondTypeScope::BoundCondType boundT;

  // convert char string to boundary condition type
  if(!strcmp(boundType, "NOBOUND")){
    boundT = BoundCondTypeScope::NOBOUND;
  }else if(!strcmp(boundType, "PRESSURE")){
    boundT = BoundCondTypeScope::PRESSURE;
  }else if(!strcmp(boundType, "RESISTANCE")){
    boundT = BoundCondTypeScope::RESISTANCE;
  }else if(!strcmp(boundType, "RCR")){
    boundT = BoundCondTypeScope::RCR;
  }else{
    return CV_ERROR;
  }

  // convert char string to boundary condition type
  if (!strcmp(lossType, "NONE")){
    loss = MinorLossScope::NONE;
  }else{
    return CV_ERROR;
  }

  // Create a new Segment
  bool IsOutlet = false;
  if (boundT != BoundCondTypeScope::NOBOUND){
    // fprintf(stdout,"  isOutlet on seg %s\n",segName);
    IsOutlet = true;
  }

  cvOneDSegment *seg = new cvOneDSegment(InitialInletArea,InitialOutletArea,
                                         InitialFlow,IsOutlet);

  //seg -> setSegmentID(ModelList[currentModel]->getNumberOfSegments());
  seg -> setSegmentID(segID);
  seg -> setSegmentName(segName);
  seg -> setParentModel((void *)&cvOneDGlobal::gModelList[cvOneDGlobal::currentModel]);
  seg -> setSegmentLength(segLen);
  seg -> setNumElements(numEls);
  seg -> setInOutJoints(inNode, outNode);
  seg -> setMaterialID(matID);
  // Set minor Loss coefficient type
  seg->SetMinorLossType(loss);

  // Set the Boundary Conditions
  seg -> setBoundCondition(boundT);
  switch(boundT) {
  case BoundCondTypeScope::RESISTANCE_TIME:
      seg->setBoundResistanceValue(value,time,num);
      break;
  case BoundCondTypeScope::RCR:
      seg->setBoundRCRValue(value,num);
      break;

  case BoundCondTypeScope::RESISTANCE:
     seg->setBoundRCRValue(value,num);//using setBoundRCRValue to set resistance and Pd.
     break;

  default:
      seg -> setBoundValue(value[0]);
      break;
  }

  seg -> setMeshType(MeshTypeScope::UNIFORM);

  cvOneDGlobal::gModelList[cvOneDGlobal::currentModel]->addSegment(seg);

  return CV_OK;
}

int cvOneDModelManager::CreateNode(char * nodeName,double x,double y,double z){

  cvOneDNode *node = new cvOneDNode();
  (node ->Name)[0] = '\0';
  strcpy(node -> Name,nodeName);
  node -> x = x;
  node -> y = y;
  node -> z = z;

  cvOneDGlobal::gModelList[cvOneDGlobal::currentModel]->addNode(node);

  return CV_OK;
}


int cvOneDModelManager::CreateJoint(char * jointName,double x,double y,double z,
                                    int numInSegs,int numOutSegs,
                                    int *InSegs,int *OutSegs){

  cvOneDJoint *joint = new cvOneDJoint();
  (joint ->Name)[0] = '\0';
  strcpy(joint -> Name,jointName);
  joint -> x = x;
  joint -> y = y;
  joint -> z = z;

  int i;
  for(i=0; i<numInSegs; i++){
    joint -> InletSegments.push_back(InSegs[i]);
  }

  for(i=0; i<numOutSegs; i++){
    joint -> OutletSegments.push_back(OutSegs[i]);
  }

  cvOneDGlobal::gModelList[cvOneDGlobal::currentModel]->addJoint(joint);

  return CV_OK;
}

// ===========
// SOLVE MODEL
// ===========
int cvOneDModelManager::SolveModel(double dt, long stepSize,
                                   long maxStep, long quadPoints,
                                   int len, char* boundType, double* values,
                                   double* times, double conv, int useIV, int usestab){

  BoundCondTypeScope::BoundCondType boundT;

  // set the creation flag to off.
  cvOneDGlobal::isCreating = false;

  // convert char string to boundary condition type
  if(!strcmp( boundType, "NOBOUND")){
    boundT = BoundCondTypeScope::NOBOUND;
    printf("Inlet Condition Type: NOBOUND\n");
  }else if(!strcmp( boundType, "PRESSURE")){
    boundT = BoundCondTypeScope::PRESSURE;
    printf("Inlet Condition Type: PRESSURE\n");
  }else if(!strcmp( boundType, "FLOW")){
    boundT = BoundCondTypeScope::FLOW;
    printf("Inlet Condition Type: FLOW\n");
  }else if(!strcmp( boundType, "RESISTANCE")){
    boundT = BoundCondTypeScope::RESISTANCE;
    printf("Inlet Condition Type: RESISTANCE\n");
  }else if(!strcmp( boundType, "RESISTANCE_TIME")){
    boundT = BoundCondTypeScope::RESISTANCE_TIME;
    printf("Inlet Condition Type: RESISTANCE_TIME\n");
  }else if(!strcmp( boundType, "RCR")){
    boundT = BoundCondTypeScope::RCR;
    printf("Inlet Condition Type: RCR\n");
  }else{
    return CV_ERROR;
  }

  // Set Solver Options
  cvOneDMthSegmentModel::STABILIZATION = usestab; // 1=stabilization, 0=none
  cvOneDGlobal::CONSERVATION_FORM = useIV;
  cvOneDBFSolver::ASCII = 1;

  cvOneDBFSolver::SetModelPtr(cvOneDGlobal::gModelList[cvOneDGlobal::currentModel]);

  // We need to get these from the solver
  cvOneDBFSolver::SetDeltaTime(dt);
  cvOneDBFSolver::SetStepSize(stepSize);
  cvOneDBFSolver::SetMaxStep(maxStep);
  cvOneDBFSolver::SetQuadPoints(quadPoints);
  cvOneDBFSolver::SetInletBCType(boundT);
  cvOneDBFSolver::DefineInletFlow(times, values, len);
  cvOneDBFSolver::SetConvergenceCriteria(conv);

  cvOneDGlobal::isSolving = true;

  cvOneDBFSolver::Solve();

  cvOneDGlobal::isSolving = false;

  return CV_OK;
}

// ================
// CREATE DATATABLE
// ================
int cvOneDModelManager::CreateDataTable(char* dtName,char* dtType, cvDoubleVec values){

  cvDoubleVec tempTime;
  cvDoubleVec tempValues;

  cvOneDDataTable* table = new cvOneDDataTable();
  table->setName(string(dtName));
  table->setType(string(dtType));
  // Get Values based on Type
  if(boost::to_upper_copy(string(dtType)) == "LIST"){
    // Values are expressed as couples of time and value
    if(values.size() % 2){
      throw cvException("ERROR: Values in data table are not in time value format.\n");
    }
    tempTime.clear();
    tempValues.clear();
    for(int loopA=0;loopA<values.size()/2;loopA++){
      tempTime.push_back(values[loopA * 2]);
      tempValues.push_back(values[loopA * 2 + 1]);
    }
    // Set Values in Table
    table->setTime(tempTime);
    table->setValues(tempValues);
    // If Debug: Show Admittance Values
    if(cvOneDGlobal::debugMode){
      printf("--- Debug\n");
      printf("%15s %15s\n","Time","Value");
      for(int loopA=0;loopA<table->getSize();loopA++){
        printf("%15e %15e\n",table->getTime(loopA),table->getValues(loopA));
      }
    }
  }else{
    throw cvException("ERROR: Invalid data table type.\n");
  }
  // ADD Data Table to the Global List
  cvOneDGlobal::gDataTables.push_back(table);
  return CV_OK;
}
