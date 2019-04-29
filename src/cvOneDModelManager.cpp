# include "cvOneDModelManager.h"

cvOneDModelManager::cvOneDModelManager(char *mdlName){
  // We're creating a model
  cvOneDGlobal::isCreating = true;

  cvOneDModel* newModel = new cvOneDModel;
  newModel->setModelName(mdlName);
  cvOneDGlobal::currentModel = cvOneDGlobal::gModelList.size();
  newModel->setModelID(cvOneDGlobal::currentModel);
  cvOneDGlobal::gModelList.push_back(newModel);
  // modelID_ = currentModel;
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
  }else if(!strcmp(boundType, "AREA")){
    boundT = BoundCondTypeScope::AREA;
  }else if(!strcmp(boundType, "FLOW")){
    boundT = BoundCondTypeScope::FLOW;
  }else if(!strcmp(boundType, "RESISTANCE")){
    boundT = BoundCondTypeScope::RESISTANCE;
  }else if(!strcmp(boundType, "RESISTANCE_TIME")){
    boundT = BoundCondTypeScope::RESISTANCE_TIME;
  }else if(!strcmp(boundType, "PRESSURE_WAVE")){
    boundT = BoundCondTypeScope::PRESSURE_WAVE;
  }else if(!strcmp(boundType, "WAVE")){
    boundT = BoundCondTypeScope::WAVE;
  }else if(!strcmp(boundType, "LPN")){
    boundT = BoundCondTypeScope::LPN;
  }else if(!strcmp(boundType, "RCR")){
    boundT = BoundCondTypeScope::RCR;
  }else if(!strcmp(boundType, "CORONARY")){
    boundT = BoundCondTypeScope::CORONARY;
  }else if(!strcmp(boundType, "ADMITTANCE")){
    boundT = BoundCondTypeScope::ADMITTANCE;
  }else if(!strcmp(boundType, "PULMONARY")){
    boundT = BoundCondTypeScope::ADMITTANCE; //Why is PULMONARY boundary condition the same as ADMITTANCE? (MLD 012518)
  }else{
    return CV_ERROR;
  }

  // convert char string to boundary condition type
  if (!strcmp(lossType, "NONE")){
    loss = MinorLossScope::NONE;
  }else if(!strcmp(lossType, "STENOSIS")){
    loss = MinorLossScope::STENOSIS;
  }else if(!strcmp(lossType, "BRANCH_THROUGH_DIVIDING")){
    loss = MinorLossScope::BRANCH_THROUGH_DIVIDING;
  }else if(!strcmp(lossType, "BRANCH_SIDE_DIVIDING")){
    loss = MinorLossScope::BRANCH_SIDE_DIVIDING;
  }else if(!strcmp(lossType, "BRANCH_THROUGH_CONVERGING")){
    loss = MinorLossScope::BRANCH_THROUGH_CONVERGING;
  }else if(!strcmp(lossType, "BRANCH_SIDE_CONVERGING")){
    loss = MinorLossScope::BRANCH_SIDE_CONVERGING;
  }else if(!strcmp(lossType, "BIFURCATION_BRANCH")){
    loss = MinorLossScope::BIFURCATION_BRANCH;
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
  // Define Branch Angle even with no pressure losses to plot VTK // DES
  seg->SetBranchAngle(branchAngle);
  if(loss != MinorLossScope::NONE){
     seg->SetUpstreamSeg(upstreamSegment);
     seg->SetBranchSeg(branchSegment);
  }

  // Set the Boundary Conditions
  seg -> setBoundCondition(boundT);
  switch(boundT) {
  case BoundCondTypeScope::PRESSURE_WAVE:
      seg->setBoundPressureValue(value,time,num);
      break;
  case BoundCondTypeScope::RESISTANCE_TIME:
      seg->setBoundResistanceValue(value,time,num);
      break;
  case BoundCondTypeScope::ADMITTANCE:
      //fprintf(stdout,"Setting ADMITTANCE boundary values (%i)\n",num);
      //printf("SEGMENT VALUES\n");
      //for(int loopA=0;loopA<num;loopA++){
      //  printf("%e\n",value[loopA]);
      //}
      //getchar();
      seg->setBoundImpedanceValue(value,num);
      break;
  case BoundCondTypeScope::RCR:
      seg->setBoundRCRValue(value,num);
      break;
  case BoundCondTypeScope::CORONARY:
      seg->SetBoundCoronaryValues(value,time,num);
      break;
  case BoundCondTypeScope::WAVE:
      seg->setBoundWaveValue(value,num);
      break;
  case BoundCondTypeScope::LPN:
      // COMPLETE DES !!!
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
  }else if(!strcmp( boundType, "AREA")){
    boundT = BoundCondTypeScope::AREA;
    printf("Inlet Condition Type: AREA\n");
  }else if(!strcmp( boundType, "FLOW")){
    boundT = BoundCondTypeScope::FLOW;
    printf("Inlet Condition Type: FLOW\n");
  }else if(!strcmp( boundType, "RESISTANCE")){
    boundT = BoundCondTypeScope::RESISTANCE;
    printf("Inlet Condition Type: RESISTANCE\n");
  }else if(!strcmp( boundType, "RESISTANCE_TIME")){
    boundT = BoundCondTypeScope::RESISTANCE_TIME;
    printf("Inlet Condition Type: RESISTANCE_TIME\n");
  }else if(!strcmp( boundType, "PRESSURE_WAVE")){
    boundT = BoundCondTypeScope::PRESSURE_WAVE;
    printf("Inlet Condition Type: PRESSURE_WAVE\n");
  }else if(!strcmp( boundType, "WAVE")){
    boundT = BoundCondTypeScope::WAVE;
    printf("Inlet Condition Type: WAVE\n");
  }else if(!strcmp( boundType, "LPN")){
    boundT = BoundCondTypeScope::LPN;
    printf("Inlet Condition Type: LPN\n");
  }else if(!strcmp( boundType, "RCR")){
    boundT = BoundCondTypeScope::RCR;
    printf("Inlet Condition Type: RCR\n");
  }else if(!strcmp( boundType, "CORONARY")){
    boundT = BoundCondTypeScope::CORONARY;
    printf("Inlet Condition Type: CORONARY\n");
  }else if(!strcmp( boundType, "ADMITTANCE")){
    boundT = BoundCondTypeScope::ADMITTANCE;
    printf("Inlet Condition Type: ADMITTANCE");
  }else if(!strcmp( boundType, "PULMONARY")){
    boundT = BoundCondTypeScope::ADMITTANCE;
    printf("Inlet Condition Type: PULMONARY\n");
  }else if(!strcmp( boundType, "CLOSEDLOOP")){    // MLD 012518
    boundT = BoundCondTypeScope::CLOSEDLOOP;
    printf("Inlet Condition Type: CLOSEDLOOP\n");
  }else{
    return CV_ERROR;
  }

  // Set Solver Options
  cvOneDMthSegmentModel::STABILIZATION = usestab; // 1=Brooke's stabilization, 0=none
  cvOneDGlobal::CONSERVATION_FORM = useIV;
  cvOneDBFSolver::ASCII = 1;

  cvOneDBFSolver::SetModelPtr(cvOneDGlobal::gModelList[cvOneDGlobal::currentModel]);

  // We need to get these from Java
  cvOneDBFSolver::SetDeltaTime(dt);
  cvOneDBFSolver::SetStepSize(stepSize);
  cvOneDBFSolver::SetMaxStep(maxStep);
  cvOneDBFSolver::SetQuadPoints(quadPoints);
  cvOneDBFSolver::SetInletBCType(boundT);
  cvOneDBFSolver::DefineInletFlow(times, values, len);
  cvOneDBFSolver::SetConvergenceCriteria(conv);

  //printf("SOLVER FLOW RATE\n");
  //for(int loopA=0;loopA<len;loopA++){
  //  printf("%e %e\n",times[loopA],values[loopA]);
  //}
  //getchar();
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
  int numTimeSteps;
  double lengthRadiusRatio;
  double rootRadius;
  double period;
  double scaleFactor;
  int rtnFourier;
  double Rd;
  double Ru;
  double Cap;
  int min_order;
  double rootRad;
  double ratio;

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
  }else if(boost::to_upper_copy(string(dtType)) == "IMPEDANCE"){
    // Check Size
    if(values.size() != 6){
      throw cvException("ERROR: Invalid number of parameters for Impedance.\n");
    }
    numTimeSteps = (int)values[0];
    lengthRadiusRatio = values[1];
    rootRadius = values[2];
    period = values[3];
    scaleFactor = values[4];
    rtnFourier = (int)values[5];
    if(CalcImpedance(numTimeSteps,lengthRadiusRatio,rootRadius,period,scaleFactor,rtnFourier,tempTime,tempValues) == CV_OK){
      // Set Values in Table
      table->setTime(tempTime);
      table->setValues(tempValues);
      // If Debug: Show Admittance Values
      if(cvOneDGlobal::debugMode){
        printf("--- Debug\n");
        printf("%15s %15s\n","Time","Impedance");
        for(int loopA=0;loopA<table->getSize();loopA++){
          printf("%15e %15e\n",table->getTime(loopA),table->getValues(loopA));
        }
      }
    }else{
      throw cvException("ERROR: Cannot Compute Impedance in CreateDataTable.\n");
    }
  }else if(boost::to_upper_copy(string(dtType)) == "RCRIMPEDANCE"){
    // Check Size
    if(values.size() != 7){
      throw cvException("ERROR: Invalid number of parameters for RCRImpedance.\n");
    }
    numTimeSteps = (int)values[0];
    Rd = values[1];
    Ru = values[2];
    period = values[3];
    Cap = values[4];
    scaleFactor = values[5];
    rtnFourier = (int)values[6];
    if(CalcImpedanceRCR(numTimeSteps,Rd,Ru,period,Cap,scaleFactor,rtnFourier,tempTime,tempValues) == CV_OK){
      // Set Values in Table
      table->setTime(tempTime);
      table->setValues(tempValues);
    }else{
      throw cvException("ERROR: Cannot Compute RCRImpedance in CreateDataTable.\n");
    }
  }else if(boost::to_upper_copy(string(dtType)) == "MORPHIMPEDANCE"){
    // Check Size
    if(values.size() != 5){
      throw cvException("ERROR: Invalid number of parameters for MORPHImpedance.\n");
    }
    numTimeSteps = (int)values[0];
    min_order = (int)values[1];
    rootRad = values[2];
    period = values[3];
    rtnFourier = (int)values[4];
    if(CalcMorphImpedance(numTimeSteps,min_order,rootRad,period,rtnFourier,tempTime,tempValues) == CV_OK){
      // Set Values in Table
      table->setTime(tempTime);
      table->setValues(tempValues);
    }else{
      throw cvException("ERROR: Cannot Compute MORPHImpedance in CreateDataTable.\n");
    }
  }else if(boost::to_upper_copy(string(dtType)) == "ADMITTANCE"){
    // Check Size
    if(values.size() != 6){
      throw cvException("ERROR: Invalid number of parameters for Admittance.\n");
    }
    numTimeSteps = (int)values[0];
    ratio = values[1];
    rootRad = values[2];
    period = values[3];
    scaleFactor = values[4];
    rtnFourier = (int)values[5];
    if(CalcAdmittance(numTimeSteps,ratio,rootRad,period,scaleFactor,rtnFourier,tempTime,tempValues) == CV_OK){
      // Set Values in Table
      table->setTime(tempTime);
      table->setValues(tempValues);
      // If Debug: Show Admittance Values
      if(cvOneDGlobal::debugMode){
        printf("--- Debug\n");
        printf("%15s %15s\n","Time","Admittance");
        for(int loopA=0;loopA<table->getSize();loopA++){
          printf("%15e %15e\n",table->getTime(loopA),table->getValues(loopA));
        }
      }
    }else{
      throw cvException("ERROR: Cannot Compute Admittance in CreateDataTable.\n");
    }
  }else if(boost::to_upper_copy(string(dtType)) == "RCRADMITTANCE"){
    // Check Size
    if(values.size() != 7){
      throw cvException("ERROR: Invalid number of parameters for RCRAdmittance.\n");
    }
    numTimeSteps = (int)values[0];
    Rd = values[1];
    Ru = values[2];
    period = values[3];
    Cap = values[4];
    scaleFactor = values[5];
    rtnFourier = (int)values[6];
    if(CalcRCRAdmittance(numTimeSteps,Rd,Ru,period,Cap,scaleFactor,rtnFourier,tempTime,tempValues) == CV_OK){
      // Set Values in Table
      table->setTime(tempTime);
      table->setValues(tempValues);
    }else{
      throw cvException("ERROR: Cannot Compute RCRAdmittance in CreateDataTable.\n");
    }
  }else if(boost::to_upper_copy(string(dtType)) == "MORPHADMITTANCE"){
    // Check Size
    if(values.size() != 5){
      throw cvException("ERROR: Invalid number of parameters for MORPHAdmittance.\n");
    }
    numTimeSteps = (int)values[0];
    min_order = (int)values[1];
    rootRad = values[2];
    period = values[3];
    rtnFourier = (int)values[4];
    if(CalcMorphAdmittance(numTimeSteps,min_order,rootRad,period,rtnFourier,tempTime,tempValues) == CV_OK){
      // Set Values in Table
      table->setTime(tempTime);
      table->setValues(tempValues);
    }else{
      throw cvException("ERROR: Cannot Compute MORPHAdmittance in CreateDataTable.\n");
    }
  }else{
    throw cvException("ERROR: Invalid data table type.\n");
  }
  // ADD Data Table to the Global List
  cvOneDGlobal::gDataTables.push_back(table);
  return CV_OK;
}

// =================
// COMPUTE IMPEDANCE
// =================
int cvOneDModelManager::CalcImpedance(int numTimeSteps, double lengthRadiusRatio, double rootRadius, double period, double scaleFactor, int rtnFourier, cvDoubleVec& time, cvDoubleVec& values){

  // Create Mette Impedance
  cvOneDMetteImpedance* mette = new cvOneDMetteImpedance();

  // Set Exercise Conditions
  bool exerciseFlag;
  if (scaleFactor > 0) {
    exerciseFlag = true;
  } else {
    exerciseFlag = false;
    scaleFactor = 1.0;
  }
  mette->SetExerciseFlag(exerciseFlag);
  mette->SetExerciseFactor(scaleFactor);

  // Compute Impedance
  double* results;
  time.clear();
  values.clear();
  if(rtnFourier == 0){
    results = mette->calculateImpedanceTime(rootRadius,lengthRadiusRatio,period,numTimeSteps);
    // NATHAN SUGGESTION !!!
    double dt = period / (numTimeSteps - 1);
    for (int i = 0; i < numTimeSteps; i++) {
      time.push_back(dt*i);
      values.push_back(results[i]);
    }
    time.push_back(period);
    values.push_back(results[0]);
  }else{
    results = mette->calculateImpedanceFourier(rootRadius,lengthRadiusRatio,period,numTimeSteps);
    for (int i = 0; i < numTimeSteps; i++) {
      values.push_back(results[i]);
    }
  }
  // Clean Memory
  delete results;
  // Return
  return CV_OK;
}

// =====================
// COMPUTE RCR IMPEDANCE
// =====================
int cvOneDModelManager::CalcImpedanceRCR(int numTimeSteps, double Rd, double Ru, double period, double Cap, double scaleFactor, int rtnFourier,
                                         cvDoubleVec& time, cvDoubleVec& values){
  // CREATE METTE IMPEDANCE
  cvOneDMetteImpedance* mette = new cvOneDMetteImpedance();

  // SET EXCERCISE FLAG
  bool exerciseFlag;
  if (scaleFactor > 0) {
    exerciseFlag = true;
  } else {
    exerciseFlag = false;
    scaleFactor = 1.0;
  }
  mette->SetExerciseFlag(exerciseFlag);
  mette->SetExerciseFactor(scaleFactor);

  // COMPUTE IMPEDANCE
  double* results;
  time.clear();
  values.clear();
  if(rtnFourier == 0){
    results = mette->calculateRCRImpedanceTime(Rd,Ru,Cap,period,numTimeSteps);
    // NATHAN SUGGESTION !!!
    double dt = period / (numTimeSteps - 1);
    for (int i = 0; i < numTimeSteps; i++) {
      time.push_back(dt*i);
      values.push_back(results[i]);
    }
    time.push_back(period);
    values.push_back(results[0]);
  }else{
    results = mette->calculateRCRImpedanceTime(Rd,Ru,Cap,period,numTimeSteps);
    for (int i = 0; i < numTimeSteps; i++) {
      values.push_back(results[i]);
    }
  }
  // Clean Memory
  delete results;
  // Return
  return CV_OK;
}

// ==============================
// COMPUTE MORPHOMETRIC IMPEDANCE
// ==============================
int cvOneDModelManager::CalcMorphImpedance(int numTimeSteps, int min_order, double rootRad, double period, int rtnFourier,
                                           cvDoubleVec& time, cvDoubleVec& values){

  // CREATE MORPHOMETRIC IMPEDANCE
  cvOneDMorphImpedance* morph = new cvOneDMorphImpedance();

  // COMPUTE IMPEDANCE
  double* results;
  time.clear();
  values.clear();
  if(rtnFourier == 0){
    results = morph->calculateImpedanceTime(rootRad,min_order,period,numTimeSteps);
    // NATHAN SUGGESTION !!!
    double dt = period / (numTimeSteps - 1);
    for(int i = 0; i < numTimeSteps; i++){
      time.push_back(dt*i);
      values.push_back(results[i]);
    }
    time.push_back(period);
    values.push_back(results[0]);
  }else{
    results = morph->calculateImpedanceFourier(rootRad,min_order,period,numTimeSteps);
    for(int i = 0; i < numTimeSteps; i++){
      values.push_back(results[i]);
    }
  }
  // Clean Memory
  delete results;
  // Return
  return CV_OK;
}

// ==================
// COMPUTE ADMITTANCE
// ==================
int cvOneDModelManager::CalcAdmittance(int numTimeSteps,double ratio,double rootRad,double period,double scaleFactor,int rtnFourier,
                                       cvDoubleVec& time, cvDoubleVec& values){

  // CREATE METTE IMPEDANCE
  cvOneDMetteImpedance* mette = new cvOneDMetteImpedance();

  // SET EXERCISE FLAG
  bool exerciseFlag;
  if (scaleFactor > 0) {
    exerciseFlag = true;
  } else {
    exerciseFlag = false;
    scaleFactor = 1.0;
  }
  mette->SetExerciseFlag(exerciseFlag);
  mette->SetExerciseFactor(scaleFactor);

  // COMPUTE ADMITTANCE
  double* results;
  time.clear();
  values.clear();
  if(rtnFourier == 0){
    results = mette->calculateAdmittanceTime(rootRad,ratio,period,numTimeSteps);
    // NATHAN SUGGESTION !!!
    double dt = period / (numTimeSteps - 1);
    for(int i = 0; i < numTimeSteps - 1; i++){
      time.push_back(dt*i);
      values.push_back(results[i]);
    }
    time.push_back(period);
    values.push_back(results[0]);
  }else{
    results = mette->calculateAdmittanceFourier(rootRad,ratio,period,numTimeSteps);
    for(int i = 0; i < numTimeSteps; i++){
      values.push_back(results[i]);
    }
  }
  // Clean Memory
  delete results;
  // Return
  return CV_OK;
}

// ======================
// COMPUTE RCR ADMITTANCE
// ======================
int cvOneDModelManager::CalcRCRAdmittance(int numTimeSteps,double Rd,double Ru,double period,double Cap,double scaleFactor,int rtnFourier,
                                          cvDoubleVec& time, cvDoubleVec& values){
  // CREATE METTE IMPEDANCE
  cvOneDMetteImpedance* mette = new cvOneDMetteImpedance();

  // SET EXERCISE FLAG
  bool exerciseFlag;
  if (scaleFactor > 0) {
    exerciseFlag = true;
  } else {
    exerciseFlag = false;
    scaleFactor = 1.0;
  }
  mette->SetExerciseFlag(exerciseFlag);
  mette->SetExerciseFactor(scaleFactor);

  // COMPUTE IMPEDANCE
  double* results;
  time.clear();
  values.clear();
  if(rtnFourier == 0){
    results = mette->calculateRCRAdmittanceTime(Rd,Ru,Cap,period,numTimeSteps);
    // NATHAN SUGGESTION !!!
    double dt = period / (numTimeSteps - 1);
    for(int i = 0; i < numTimeSteps; i++){
      time.push_back(dt*i);
      values.push_back(results[i]);
    }
    time.push_back(period);
    values.push_back(results[0]);
  }else{
    results = mette->calculateRCRAdmittanceFourier(Rd,Ru,Cap,period,numTimeSteps);
    for(int i = 0; i < numTimeSteps; i++){
      values.push_back(results[i]);
    }
  }
  // Clean Memory
  delete results;
  // Return
  return CV_OK;
}

// ===============================
// COMPUTE MORPHOMETRIC ADMITTANCE
// ===============================
int cvOneDModelManager::CalcMorphAdmittance(int numTimeSteps,int min_order,double rootRad,double period,int rtnFourier,
                                            cvDoubleVec& time, cvDoubleVec& values){

  // Create Morphometric Impedance
  cvOneDMorphImpedance* morph = new cvOneDMorphImpedance();

  // COMPUTE IMPEDANCE
  double* results;
  time.clear();
  values.clear();
  if(rtnFourier == 0){
    results = morph->calculateAdmittanceTime(rootRad,min_order,period,numTimeSteps);
    // NATHAN SUGGESTION !!!
    double dt = period / (numTimeSteps - 1);
    for(int i = 0; i < numTimeSteps; i++){
      time.push_back(dt*i);
      values.push_back(results[i]);
    }
    time.push_back(period);
    values.push_back(results[0]);
  }else{
    results = morph->calculateAdmittanceFourier(rootRad,min_order,period,numTimeSteps);
    for(int i = 0; i < numTimeSteps; i++){
      values.push_back(results[i]);
    }
  }
  // Clean Memory
  delete results;
  // Return
  return CV_OK;
}

// ===========
// COMPUTE LRR
// ===========
int cvOneDModelManager::CalcLRR(int numTimeSteps,double resistance,double rootRad,double period,double scaleFactor,
                                double& lrr){

  // CREATE METTE IMPEDANCE
  cvOneDMetteImpedance* mette = new cvOneDMetteImpedance();

  // SET EXERCISE FLAG
  bool exerciseFlag;
  if (scaleFactor > 0) {
    exerciseFlag = true;
  } else {
    exerciseFlag = false;
    scaleFactor = 1.0;
  }
  mette->SetExerciseFlag(exerciseFlag);
  mette->SetExerciseFactor(scaleFactor);

  lrr = 0.0;
  lrr = mette->resistBasedTree(rootRad,resistance,period,numTimeSteps);

  // Return
  return CV_OK;
}
