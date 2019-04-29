# include "main.h"

using namespace std;

// ====================
// WRITE PROGRAM HEADER
// ====================
void WriteHeader(){
  printf("---------------------------------\n");
  printf(" oneDSolver \n");
  printf(" 1D Finite Element Hemodynamics  \n");
  printf("---------------------------------\n");
}

// ====================================
// GET DATA TABLE ENTRY FROM STRING KEY
// ====================================
int getDataTableIDFromStringKey(string key){
  bool found = false;
  int count = 0;
  while((!found)&&(count<cvOneDGlobal::gDataTables.size())){
    found = (boost::to_upper_copy(key) == boost::to_upper_copy(cvOneDGlobal::gDataTables[count]->getName()));
    // Update Counter
    if(!found){
      count++;
    }
  }
  if(!found){
    throw cvException(string("ERROR: Cannot find data table entry from string key: " + key + ".\n").c_str());
    return -1;
  }else{
    return count;
  }
}

// ===============================
// CREATE MODEL AND RUN SIMULATION
// ===============================
void createAndRunModel(cvOneDOptions* opts){

  // MESSAGE
  printf("\n");
  printf("Creating and Running Model ...\n");

  // CREATE MODEL MANAGER
  cvOneDModelManager* oned = new cvOneDModelManager((char*)opts->modelName.c_str());

  // CREATE NODES
  printf("Creating Nodes ... \n");
  int totNodes = opts->nodeName.size();
  int nodeError = CV_OK;
  for(int loopA=0;loopA<totNodes;loopA++){
    // Finally Create Joint
    nodeError = oned->CreateNode((char*)opts->nodeName[loopA].c_str(),
                                   opts->nodeXcoord[loopA], opts->nodeYcoord[loopA], opts->nodeZcoord[loopA]);
    if(nodeError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating NODE " + to_string(loopA) + "\n").c_str());
    }
  }

  // CREATE JOINTS
  printf("Creating Joints ... \n");
  int totJoints = opts->jointName.size();
  int jointError = CV_OK;
  int* asInlets = NULL;
  int* asOutlets = NULL;
  string currInletName;
  string currOutletName;
  int jointInletID = 0;
  int jointOutletID = 0;
  int totJointInlets = 0;
  int totJointOutlets = 0;
  for(int loopA=0;loopA<totJoints;loopA++){
    // GET NAMES FOR INLET AND OUTLET
    currInletName = opts->jointInletName[loopA];
    currOutletName = opts->jointOutletName[loopA];
    // FIND JOINTINLET INDEX
    jointInletID = getListIDWithStringKey(currInletName,opts->jointInletListNames);
    if(jointInletID < 0){
      throw cvException(string("ERROR: Cannot Find JOINTINLET for key " + currInletName).c_str());
    }
    totJointInlets = opts->jointInletListNumber[jointInletID];
    // FIND JOINTOUTLET INDEX
    jointOutletID = getListIDWithStringKey(currOutletName,opts->jointOutletListNames);
    if(jointInletID < 0){
      throw cvException(string("ERROR: Cannot Find JOINTOUTLET for key " + currOutletName).c_str());
    }
    // GET TOTALS
    totJointInlets = opts->jointInletListNumber[jointInletID];
    totJointOutlets = opts->jointOutletListNumber[jointOutletID];
    // ALLOCATE INLETS AND OUTLET LIST
    asInlets = NULL;
    asOutlets = NULL;
    if(totJointInlets > 0){
      asInlets = new int[totJointInlets];
      for(int loopB=0;loopB<totJointInlets;loopB++){
        asInlets[loopB] = opts->jointInletList[jointInletID][loopB];
      }
    }
    if(totJointOutlets > 0){
      asOutlets = new int[totJointOutlets];
      for(int loopB=0;loopB<totJointOutlets;loopB++){
        asOutlets[loopB] = opts->jointOutletList[jointOutletID][loopB];
      }
    }
    // Finally Create Joint
    jointError = oned->CreateJoint((char*)opts->jointName[loopA].c_str(),
                                   opts->nodeXcoord[loopA], opts->nodeYcoord[loopA], opts->nodeZcoord[loopA],
                                   totJointInlets, totJointOutlets,asInlets,asOutlets);
    if(jointError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating JOINT " + to_string(loopA) + "\n").c_str());
    }
    // Deallocate
    delete [] asInlets;
    delete [] asOutlets;
    asInlets = NULL;
    asOutlets = NULL;
  }

  // CREATE MATERIAL
  printf("Creating Materials ... \n");
  int totMaterials = opts->materialName.size();
  int matError = CV_OK;
  double doubleParams[3];
  int matID = 0;
  string currMatType = "MATERIAL_OLUFSEN";
  int numParams = 0;
  for(int loopA=0;loopA<totMaterials;loopA++){
    if(boost::to_upper_copy(opts->materialType[loopA]) == "OLUFSEN"){
      currMatType = "MATERIAL_OLUFSEN";
      numParams = 3;
    }else{
      currMatType = "MATERIAL_LINEAR";
      numParams = 1;
    }
    doubleParams[0] = opts->materialParam1[loopA];
    doubleParams[1] = opts->materialParam2[loopA];
    doubleParams[2] = opts->materialParam3[loopA];
    // CREATE MATERIAL
    matError = oned->CreateMaterial((char*)opts->materialName[loopA].c_str(),
                                    (char*)currMatType.c_str(),
                                    opts->materialDensity[loopA],
                                    opts->materialViscosity[loopA],
                                    opts->materialExponent[loopA],
                                    opts->materialPRef[loopA],
                                    numParams, doubleParams,
                                    &matID);
    if(matError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating MATERIAL " + to_string(loopA) + "\n").c_str());
    }

  }

  // CREATE DATATABLES
  printf("Creating Data Tables ... \n");
  int totCurves = opts->dataTableName.size();
  int curveError = CV_OK;
  for(int loopA=0;loopA<totCurves;loopA++){
    curveError = oned->CreateDataTable((char*)opts->dataTableName[loopA].c_str(),(char*)opts->dataTableType[loopA].c_str(),opts->dataTableVals[loopA]);
    if(curveError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating DATATABLE " + to_string(loopA) + "\n").c_str());
    }
  }

  // SEGMENT DATA
  printf("Creating Segments ... \n");
  int segmentError = CV_OK;
  int totalSegments = opts->segmentName.size();
  int curveTotals = 0;
  double* curveTime = NULL;
  double* curveValue = NULL;
  string matName;
  string curveName;
  int currMatID = 0;
  int dtIDX = 0;
  for(int loopA=0;loopA<totalSegments;loopA++){

    // GET MATERIAL
    matName = opts->segmentMatName[loopA];
    currMatID = getListIDWithStringKey(matName,opts->materialName);
    if(currMatID < 0){
      throw cvException(string("ERROR: Cannot Find Material for key " + matName).c_str());
    }

    // GET CURVE DATA
    curveName = opts->segmentDataTableName[loopA];
    if(boost::to_upper_copy(curveName) != "NONE"){
      dtIDX = getDataTableIDFromStringKey(curveName);
      curveTotals = cvOneDGlobal::gDataTables[dtIDX]->getSize();
      curveTime = new double[curveTotals];
      curveValue = new double[curveTotals];
      for(int loopA=0;loopA<curveTotals;loopA++){
        curveTime[loopA] = cvOneDGlobal::gDataTables[dtIDX]->getTime(loopA);
        curveValue[loopA] = cvOneDGlobal::gDataTables[dtIDX]->getValues(loopA);
      }
    }else{
      curveTotals = 1;
      curveTime = new double[curveTotals];
      curveValue = new double[curveTotals];
      curveTime[0] = 0.0;
      curveValue[0] = 0.0;
    }
    segmentError = oned->CreateSegment((char*)opts->segmentName[loopA].c_str(),
                                       (long)opts->segmentID[loopA],
                                       opts->segmentLength[loopA],
                                       (long)opts->segmentTotEls[loopA],
                                       (long)opts->segmentInNode[loopA],
                                       (long)opts->segmentOutNode[loopA],
                                       opts->segmentInInletArea[loopA],
                                       opts->segmentInOutletArea[loopA],
                                       opts->segmentInFlow[loopA],
                                       currMatID,
                                       (char*)opts->segmentLossType[loopA].c_str(),
                                       opts->segmentBranchAngle[loopA],
                                       opts->segmentUpstreamSegment[loopA],
                                       opts->segmentBranchSegment[loopA],
                                       (char*)opts->segmentBoundType[loopA].c_str(),
                                       curveValue,
                                       curveTime,
                                       curveTotals);
    if(segmentError == CV_ERROR){
      throw cvException(string("ERROR: Error Creating SEGMENT " + to_string(loopA) + "\n").c_str());
    }
    // Deallocate
    delete [] curveTime;
    curveTime = NULL;
    delete [] curveValue;
    curveValue = NULL;
  }

  // PRINT MODEL FOR DEBUG
  //cvOneDGlobal::gModelList[0]->PrintModel();

  double* vals;
  int tot;
  cvOneDGlobal::gModelList[cvOneDGlobal::currentModel]->getSegment(0)->getBoundImpedanceValues(&vals,&tot);
  //printf("TEST UP HERE 1 ...");
  //for(int loopA=0;loopA<tot;loopA++){
  //  printf("%e\n",vals[loopA]);
  //}
  //getchar();

  // SOLVE MODEL
  printf("Solving Model ... \n");
  int solveError = CV_OK;
  string inletCurveName = opts->inletDataTableName;
  int inletCurveIDX = getDataTableIDFromStringKey(inletCurveName);
  int inletCurveTotals = cvOneDGlobal::gDataTables[inletCurveIDX]->getSize();
  double* inletCurveTime = new double[inletCurveTotals];
  double* inletCurveValue = new double[inletCurveTotals];
  for(int loopB=0;loopB<inletCurveTotals;loopB++){
    inletCurveTime[loopB] = cvOneDGlobal::gDataTables[inletCurveIDX]->getTime(loopB);
    inletCurveValue[loopB] = cvOneDGlobal::gDataTables[inletCurveIDX]->getValues(loopB);
  }
  // Solve Model
  solveError = oned->SolveModel(opts->timeStep,
                                opts->stepSize,
                                opts->maxStep,
                                opts->quadPoints,
                                inletCurveTotals,
                                (char*)opts->boundaryType.c_str(),
                                inletCurveValue,
                                inletCurveTime,
                                opts->convergenceTolerance,
                                // Formulation Type
                                opts->useIV,
                                // Stabilization
                                opts->useStab,
                                //shockcapture viscosity
                                opts->useShockcap,
                                //shockcapture smooth parameter beta
                                opts->smoothbeta,
                                //shochcapture S,Q ref values
                                opts->Sref,
                                opts->Qref
                                );
  if(solveError == CV_ERROR){
    throw cvException(string("ERROR: Error Solving Model\n").c_str());
  }
  delete [] inletCurveTime;
  delete [] inletCurveValue;
}

// ======================
// READ SINGLE MODEL FILE
// ======================
void readModelFile(string inputFile, cvOneDOptions* opts, cvStringVec includedFiles){

  // Message
  printf("\n");
  printf("Reading file: %s ... \n",inputFile.c_str());

  // Declare
  cvStringVec tokenizedString;
  cvLongVec   tempIntVec;
  string      matType;
  cvDoubleVec temp;
  bool doInclude = false;

  // Declare input File
  ifstream infile;
  infile.open(inputFile);
  if(infile.fail()){
    throw cvException("ERROR: Input file does not exist.\n");
  }
  int lineCount = 1;
  int reminder = 0;
  int totSegments = 0;

  // Read Data From File
  std::string buffer;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    // Check for Empty buffer
    if(!buffer.empty()){
      // CHECK THE ELEMENT TYPE
      if(boost::to_upper_copy(tokenizedString[0]) == string("MODEL")){
        //printf("Found Model.\n");
        if(opts->modelNameDefined){
          throw cvException("ERROR: Model name already defined\n");
        }
        if(tokenizedString.size() > 2){
          throw cvException(string("ERROR: Too many parameters for MODEL token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 2){
          throw cvException(string("ERROR: Not enough parameters for MODEL token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          opts->modelName = tokenizedString[1];
        }catch(...){
          throw cvException(string("ERROR: Invalid Model Name. Line " + to_string(lineCount) + "\n").c_str());
        }
        opts->modelNameDefined = true;
      }else if(boost::to_upper_copy(tokenizedString[0]) == string("NODE")){
        // printf("Found Joint.\n");
        if(tokenizedString.size() > 5){
          throw cvException(string("ERROR: Too many parameters for NODE token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 5){
          throw cvException(string("ERROR: Not enough parameters for NODE token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Get Node Name
          opts->nodeName.push_back(tokenizedString[1]);
          opts->nodeXcoord.push_back(atof(tokenizedString[2].c_str()));
          opts->nodeYcoord.push_back(atof(tokenizedString[3].c_str()));
          opts->nodeZcoord.push_back(atof(tokenizedString[4].c_str()));
        }catch(...){
          throw cvException(string("ERROR: Invalid NODE Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(boost::to_upper_copy(tokenizedString[0]) == string("JOINT")){
        // printf("Found Joint.\n");
        if(tokenizedString.size() > 5){
          throw cvException(string("ERROR: Too many parameters for JOINT token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 5){
          throw cvException(string("ERROR: Not enough parameters for JOINT token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Get Joint Name
          opts->jointName.push_back(tokenizedString[1]);
//          opts->jointNode.push_back(atof(tokenizedString[2].c_str()));
          opts->jointNode.push_back(tokenizedString[2]);
          opts->jointInletName.push_back(tokenizedString[3]);
          opts->jointOutletName.push_back(tokenizedString[4]);
        }catch(...){
          throw cvException(string("ERROR: Invalid JOINT Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(boost::to_upper_copy(tokenizedString[0]) == string("JOINTINLET")){
        // printf("Found JointInlet.\n");
        if(tokenizedString.size() < 3){
          throw cvException(string("ERROR: Not enough parameters for JOINTINLET token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Get Joint Name
          opts->jointInletListNames.push_back(tokenizedString[1]);
          totSegments = atoi(tokenizedString[2].c_str());
          opts->jointInletListNumber.push_back(totSegments);
          tempIntVec.clear();
          if(totSegments > 0){
            for(size_t loopA=3;loopA<tokenizedString.size();loopA++){
              tempIntVec.push_back(atoi(tokenizedString[loopA].c_str()));
            }
          }
          opts->jointInletList.push_back(tempIntVec);
        }catch(...){
          throw cvException(string("ERROR: Invalid JOINTINLET Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("JOINTOUTLET")){
        // printf("Found JointOutlet.\n");
        if(tokenizedString.size() < 3){
          throw cvException(string("ERROR: Not enough parameters for JOINTOUTLET token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Get Joint Name
          opts->jointOutletListNames.push_back(tokenizedString[1]);
          totSegments = atoi(tokenizedString[2].c_str());
          opts->jointOutletListNumber.push_back(totSegments);
          tempIntVec.clear();
          if(totSegments > 0){
            for(size_t loopA=3;loopA<tokenizedString.size();loopA++){
              tempIntVec.push_back(atoi(tokenizedString[loopA].c_str()));
            }
          }
          opts->jointOutletList.push_back(tempIntVec);
        }catch(...){
          throw cvException(string("ERROR: Invalid JOINTOUTLET Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SEGMENT")){
        // printf("Found Segment.\n");
        if(tokenizedString.size() > 17){
          throw cvException(string("ERROR: Too many parameters for SEGMENT token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 17){
          throw cvException(string("ERROR: Not enough parameters for SEGMENT token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // char* segName,
          opts->segmentName.push_back(tokenizedString[1]);
          // long segID,
          opts->segmentID.push_back(atoi(tokenizedString[2].c_str()));
          // double  segLen,
          opts->segmentLength.push_back(atof(tokenizedString[3].c_str()));
          // long    numEls,
          opts->segmentTotEls.push_back(atoi(tokenizedString[4].c_str()));
          // long    inNode,
          opts->segmentInNode.push_back(atoi(tokenizedString[5].c_str()));
          // long    outNode,
          opts->segmentOutNode.push_back(atoi(tokenizedString[6].c_str()));
          // double  InitialInletArea,
          opts->segmentInInletArea.push_back(atof(tokenizedString[7].c_str()));
          // double  InitialOutletArea,
          opts->segmentInOutletArea.push_back(atof(tokenizedString[8].c_str()));
          // double  InitialFlow,
          opts->segmentInFlow.push_back(atof(tokenizedString[9].c_str()));
          // int matID,
          opts->segmentMatName.push_back(tokenizedString[10].c_str());
          // char* lossType,
          opts->segmentLossType.push_back(tokenizedString[11]);
          // double branchAngle,
          opts->segmentBranchAngle.push_back(atof(tokenizedString[12].c_str()));
          // int upstreamSegment,
          opts->segmentUpstreamSegment.push_back(atoi(tokenizedString[13].c_str()));
          // int branchSegment,
          opts->segmentBranchSegment.push_back(atoi(tokenizedString[14].c_str()));
          // char* boundType,
          opts->segmentBoundType.push_back(tokenizedString[15]);
          // Curve ID Instead of num,value,time
          // double* value,
          // double* time,
          // int num
          opts->segmentDataTableName.push_back(tokenizedString[16]);
        }catch(...){
          throw cvException(string("ERROR: Invalid SEGMENT Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SOLVEROPTIONS")){
        // printf("Found Solver Options.\n");
        if(opts->solverOptionDefined){
          throw cvException("ERROR: SOLVEROPTIONS already defined\n");
        }
        if(tokenizedString.size() > 14){
          throw cvException(string("ERROR: Too many parameters for SOLVEROPTIONS token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 10){
          throw cvException(string("ERROR: Not enough parameters for SOLVEROPTIONS token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // double dt,
          opts->timeStep = atof(tokenizedString[1].c_str());
          // long stepSize,
          opts->stepSize = atof(tokenizedString[2].c_str());
          // long maxStep,
          opts->maxStep = atof(tokenizedString[3].c_str());
          // long quadPoints,
          opts->quadPoints = atoi(tokenizedString[4].c_str());
          // int CurveID,
          opts->inletDataTableName = tokenizedString[5].c_str();
          // char* boundType,
          opts->boundaryType = tokenizedString[6];
          // double conv,
          opts->convergenceTolerance = atof(tokenizedString[7].c_str());
          // int useIV,
          opts->useIV = atoi(tokenizedString[8].c_str());
          // int usestab
          opts->useStab = atoi(tokenizedString[9].c_str());
          //int useshockcap
          opts->useShockcap=0;
        }catch(...){
          throw cvException(string("ERROR: Invalid SOLVEROPTIONS Format. Line " + to_string(lineCount) + "\n").c_str());
        }
        if (tokenizedString.size() == 14){
          opts->useShockcap= atoi(tokenizedString[10].c_str());
          opts->smoothbeta= atof(tokenizedString[11].c_str());
          opts->Sref= atof(tokenizedString[12].c_str());
          opts->Qref= atof(tokenizedString[13].c_str());
         }
        opts->solverOptionDefined = true;
      }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("OUTPUT")){
        if(tokenizedString.size() > 3){
          throw cvException(string("ERROR: Too many parameters for OUTPUT token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 2){
          throw cvException(string("ERROR: Not enough parameters for OUTPUT token. Line " + to_string(lineCount) + "\n").c_str());
        }
        // Output Type
        if(boost::to_upper_copy(tokenizedString[1]) == "TEXT"){
          cvOneDGlobal::outputType = OutputTypeScope::OUTPUT_TEXT;
        }else if(boost::to_upper_copy(tokenizedString[1]) == "VTK"){
          cvOneDGlobal::outputType = OutputTypeScope::OUTPUT_VTK;
        }else if(boost::to_upper_copy(tokenizedString[1]) == "BOTH"){
          cvOneDGlobal::outputType = OutputTypeScope::OUTPUT_BOTH;
        }else{
          throw cvException("ERROR: Invalid OUTPUT Type.\n");
        }
        if(tokenizedString.size() > 2){
          cvOneDGlobal::vtkOutputType = atoi(tokenizedString[2].c_str());
          if(cvOneDGlobal::vtkOutputType > 1){
            throw cvException("ERROR: Invalid OUTPUT VTK Type.\n");
          }
        }
      }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("DATATABLE")){
        // printf("Found Data Table.\n");
        try{

          // Get Datatable Name
          opts->dataTableName.push_back(tokenizedString[1]);
          // Add the type of the datatable
          opts->dataTableType.push_back(tokenizedString[2]);

          bool foundEnd = false;
          temp.clear();
          while(!foundEnd){
            std::getline(infile,buffer);
            lineCount++;
            // Trim String
            boost::trim(buffer);
            // Tokenize String
            boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
            // Check for Empty buffer
            if(!buffer.empty()){
              if(boost::to_upper_copy(tokenizedString[0]) == std::string("ENDDATATABLE")){
                foundEnd = true;
              }else{
                for(int loopA=0;loopA<tokenizedString.size();loopA++){
                  temp.push_back(atof(tokenizedString[loopA].c_str()));
                }
              }
            }
          }
          // Add all the values to the option array
          opts->dataTableVals.push_back(temp);
        }catch(...){
          throw cvException(string("ERROR: Invalid DATATABLE Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("MATERIAL")){
        // printf("Found Material.\n");
        if(tokenizedString.size() > 10){
          throw cvException(string("ERROR: Too many parameters for MATERIAL token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 8){
          throw cvException(string("ERROR: Not enough parameters for MATERIAL token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Material Name
          opts->materialName.push_back(tokenizedString[1]);
          // Material Type
          matType = tokenizedString[2];
          opts->materialType.push_back(matType);
          // Density
          opts->materialDensity.push_back(atof(tokenizedString[3].c_str()));
          // Dynamic Viscosity
          opts->materialViscosity.push_back(atof(tokenizedString[4].c_str()));
          // Reference Pressure
          opts->materialPRef.push_back(atof(tokenizedString[5].c_str()));
          // Material Exponent
          opts->materialExponent.push_back(atof(tokenizedString[6].c_str()));
          // Extra Material Parameters
          if(boost::to_upper_copy(matType) == "OLUFSEN"){
            opts->materialParam1.push_back(atof(tokenizedString[7].c_str()));
            opts->materialParam2.push_back(atof(tokenizedString[8].c_str()));
            opts->materialParam3.push_back(atof(tokenizedString[9].c_str()));
          }else if(boost::to_upper_copy(matType) == "LINEAR"){
            opts->materialParam1.push_back(atof(tokenizedString[7].c_str()));
            opts->materialParam2.push_back(0.0);
            opts->materialParam3.push_back(0.0);
          }else{
            throw cvException(string("ERROR: Invalid MATERIAL Type. Line " + to_string(lineCount) + "\n").c_str());
          }
        }catch(...){
          throw cvException("ERROR: Invalid MATERIAL Format.\n");
        }
      }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("INCLUDE")){
        // Check if the file is active
        if(boost::to_upper_copy(tokenizedString[2]) == "TRUE"){
          doInclude = true;
        }else if(boost::to_upper_copy(tokenizedString[2]) == "FALSE"){
          doInclude = false;
        }else{
          throw cvException(string("ERROR: Invalid INCLUDE switch format. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // If active include file in list
          if(doInclude){
            includedFiles.push_back(tokenizedString[1]);
          }
        }catch(...){
          throw cvException(string("ERROR: Invalid INCLUDE Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if((tokenizedString.size() == 0)||(tokenizedString[0].at(0) == '#')||(tokenizedString[0].find_first_not_of(' ') == std::string::npos)){
        // printf("Found Blank.\n");
        // COMMENT OR BLANK LINE: DO NOTHING
      }else{
        // TOKEN NOT RECOGNIZED
        throw cvException(string("ERROR: Invalid Token in input file, line: "  + to_string(lineCount) + "\n").c_str());
      }
    }
    // printf("Line: %d, Buffer: %s\n",lineCount,buffer.c_str());
    // getchar();

    // Increment Line Count
    lineCount++;
  }
  // Close File
  infile.close();
}

// ====================
// READ MODEL FROM FILE
// ====================
void readModel(string inputFile, cvOneDOptions* opts){

  // List of included Files
  cvStringVec includedFiles;
  string currentFile;

  // Read First File
  readModelFile(inputFile,opts,includedFiles);

  //Read Nested Files
  while(includedFiles.size() > 0){

    // Get the first file Name
    currentFile = includedFiles[0];
    // Delete the First element
    includedFiles.erase(includedFiles.begin());
    // Read the file and store new included files
    readModelFile(inputFile,opts,includedFiles);
  }
}

// ==============
// RUN ONEDSOLVER
// ==============
void runOneDSolver(string inputFile){

  // Create Options
  cvOneDOptions* opts = new cvOneDOptions();

  // Read Model From File
  readModel(inputFile,opts);

  // Model Checking
  opts->check();

  // Print Input Data Echo
  string fileName("echo.out");
  opts->printToFile(fileName);

  // Create Model and Run Simulation
  createAndRunModel(opts);

  // Delete Options
  delete opts;
}

// =============
// MAIN FUNCTION
// =============
int main(int argc, char** argv){

  // Write Program Header
  WriteHeader();

  try{

    // Run Simulation
    string inputFile(argv[1]);
    runOneDSolver(inputFile);

  }catch(exception& e){
    // Print Exception Message
    printf("%s\n",e.what());
    // Execution Terminated
    printf("Terminated.\n");
    return -1;
  }
  printf("Completed!\n");
  return 0;

}
