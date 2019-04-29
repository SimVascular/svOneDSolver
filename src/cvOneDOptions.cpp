#include "cvOneDOptions.h"

// CONSTRUCTOR
cvOneDOptions::cvOneDOptions(){
  modelNameDefined = false;
  solverOptionDefined = false;
}

// DESTRUCTOR
cvOneDOptions::~cvOneDOptions(){
}

// PRINT MODEL NAME
void cvOneDOptions::printModelName(FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"MODEL NAME: %s\n",modelName.c_str());
}

// PRINT NODE DATA
void cvOneDOptions::printNodeData(FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"NODE DATA\n");
  for(long int loopA=0;loopA<nodeName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"NODE: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",nodeName[loopA].c_str());
    fprintf(f,"X: %f\n",nodeXcoord[loopA]);
    fprintf(f,"Y: %f\n",nodeYcoord[loopA]);
    fprintf(f,"Z: %f\n",nodeZcoord[loopA]);
  }
}

// PRINT JOINT DATA
void cvOneDOptions::printJointData(FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"JOINT DATA\n");
  for(long int loopA=0;loopA<jointName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"JOINT: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",jointName[loopA].c_str());
    fprintf(f,"NODE: %s\n",jointNode[loopA].c_str());
    fprintf(f,"X: %f\n",nodeXcoord[loopA]);
    fprintf(f,"Y: %f\n",nodeYcoord[loopA]);
    fprintf(f,"Z: %f\n",nodeZcoord[loopA]);
    fprintf(f,"INLET NAME: %s\n",jointInletName[loopA].c_str());
    fprintf(f,"OUTLET NAME: %s\n",jointOutletName[loopA].c_str());
  }
}

// PRINT SEGMEMENT DATA
void cvOneDOptions::printSegmentData(FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"SEGMENT DATA\n");
  for(long int loopA=0;loopA<segmentName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"SEGMENT N.: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",segmentName[loopA].c_str());
    fprintf(f,"ID: %ld\n",segmentID[loopA]);
    fprintf(f,"LENGTH: %f\n",segmentLength[loopA]);
    fprintf(f,"NUM ELEMENTS: %ld\n",segmentTotEls[loopA]);
    fprintf(f,"IN NODE: %ld\n",segmentInNode[loopA]);
    fprintf(f,"OUT NODE: %ld\n",segmentOutNode[loopA]);
    fprintf(f,"IN AREA: %f\n",segmentInInletArea[loopA]);
    fprintf(f,"OUT AREA: %f\n",segmentInOutletArea[loopA]);
    fprintf(f,"IN FLOW: %f\n",segmentInFlow[loopA]);
    fprintf(f,"MAT DATATABLE NAME: %s\n",segmentMatName[loopA].c_str());
    fprintf(f,"LOSS TYPE: %s\n",segmentLossType[loopA].c_str());
    fprintf(f,"BRANCH ANGLE: %f\n",segmentBranchAngle[loopA]);
    fprintf(f,"UPSTREAM SEGMENT: %ld\n",segmentUpstreamSegment[loopA]);
    fprintf(f,"BRANCH SEGMENT: %ld\n",segmentBranchSegment[loopA]);
    fprintf(f,"BOUNDARY CONDITIONS TYPE: %s\n",segmentBoundType[loopA].c_str());
    fprintf(f,"INLET DATA TABLE: %s\n",segmentDataTableName[loopA].c_str());
  }
}

// PRINT SOLVER OPTION DATA
void cvOneDOptions::printSolverOptions(FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"PRINT SOLVER OPTION DATA\n");
  fprintf(f,"TIME STEP: %f\n",timeStep);
  fprintf(f,"STEP SIZE: %ld\n",stepSize);
  fprintf(f,"MAX STEP: %ld\n",maxStep);
  fprintf(f,"QUADRATURE POINTS: %ld\n",quadPoints);
  fprintf(f,"INLET DATA TABLE: %s\n",inletDataTableName.c_str());
  fprintf(f,"BOUNDARY TYPE: %s\n",boundaryType.c_str());
  fprintf(f,"CONVERGENCE TOLERANCE: %f\n",convergenceTolerance);
  fprintf(f,"USE IV: %d\n",useIV);
  fprintf(f,"USE STABILIZATION: %d\n",useStab);
  fprintf(f,"USE SHOCKCAPTURE: %d\n",useShockcap);
  if (useShockcap==1){
  fprintf(f,"Smooth parameter beta: %f\n",smoothbeta);
  fprintf(f,"Reference S: %f\n",Sref);
  fprintf(f,"Reference Q: %f\n",Qref);
   }

}

// PRINT MATERIAL DATA
void cvOneDOptions::printMaterialData(FILE* f){
  fprintf(f,"--- \n");
  // MATERIAL
  for(long int loopA=0;loopA<materialName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"MATERIAL NUM: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",materialName[loopA].c_str());
    fprintf(f,"TYPE: %s\n",materialType[loopA].c_str());
    fprintf(f,"DENSITY: %f\n",materialDensity[loopA]);
    fprintf(f,"VISCOSITY: %f\n",materialViscosity[loopA]);
    fprintf(f,"REF PRESSURE: %f\n",materialPRef[loopA]);
    fprintf(f,"EXPONENT: %f\n",materialExponent[loopA]);
    fprintf(f,"PARAM 1: %f\n",materialParam1[loopA]);
    fprintf(f,"PARAM 2: %f\n",materialParam2[loopA]);
    fprintf(f,"PARAM 3: %f\n",materialParam3[loopA]);
  }
}

// PRINT JOINTINLET DATA
void cvOneDOptions::printJointInletData(FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"JOINTINLET DATA\n");
  for(long int loopA=0;loopA<jointInletListNames.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"JOINTINLET NUM: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",jointInletListNames[loopA].c_str());
    fprintf(f,"TOTAL : %ld\n",jointInletListNumber[loopA]);
    if(jointInletListNumber[loopA] > 0){
      fprintf(f,"LIST ITEMS\n");
      for(size_t loopB=0;loopB<jointInletList[loopA].size();loopB++){
        fprintf(f,"%ld \n",jointInletList[loopA][loopB]);
      }
    }
    fprintf(f,"\n");
  }
}

// PRINT JOINTOPUTLET DATA
void cvOneDOptions::printJointOutletData(FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"JOINTOUTLET DATA\n");
  for(long int loopA=0;loopA<jointOutletListNames.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"JOINTOUTLET NUM: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",jointOutletListNames[loopA].c_str());
    fprintf(f,"TOTAL : %ld\n",jointOutletListNumber[loopA]);
    if(jointOutletListNumber[loopA] > 0){
      fprintf(f,"LIST ITEMS\n");
      for(size_t loopB=0;loopB<jointOutletList[loopA].size();loopB++){
        fprintf(f,"%ld \n",jointOutletList[loopA][loopB]);
      }
    }
    fprintf(f,"\n");
  }
}

// PRINT DATA TABLES
void cvOneDOptions::printDataTables(FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"DATA TABLES\n");
  for(long int loopA=0;loopA<dataTableName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"DATA TABLE ID: %ld\n",loopA);
    fprintf(f,"DATA TABLE NAME: %s\n",dataTableName[loopA].c_str());
    fprintf(f,"DATA TABLE TYPE: %s\n",dataTableType[loopA].c_str());
    for(size_t loopB=0;loopB<dataTableVals[loopA].size();loopB++){
      fprintf(f,"VALUE n. %ld: %f\n",loopB+1,dataTableVals[loopA][loopB]);
    }
  }
}

void cvOneDOptions::printToFile(string fileName){
  FILE* f;
  f = fopen(fileName.c_str(),"w");

  printf("\n");
  printf("Printing Model Echo ... \n");
  fflush(stdout);
  printModelName(f);
  //printf("1 ... ");
  //fflush(stdout);
  printNodeData(f);
  //printf("1 ... ");
  //fflush(stdout);
  printJointData(f);
  //printf("2... ");
  //fflush(stdout);
  printSegmentData(f);
  //printf("3 ... ");
  //fflush(stdout);
  printSolverOptions(f);
  //printf("4 ... ");
  //fflush(stdout);
  printMaterialData(f);
  //printf("5 ... ");
  //fflush(stdout);
  printJointInletData(f);
  //printf("6 ... ");
  //fflush(stdout);
  printJointOutletData(f);
  //printf("7 ... ");
  //fflush(stdout);
  printDataTables(f);
  //printf("8 ... ");
  //fflush(f);
  fclose(f);
}

void cvOneDOptions::checkSegmentLengthConsistency(){
  bool inconsistencyFound = false;
  // Get total number of segments
  int totSegs = segmentLength.size();
  double segLength = 0.0;
  int inNode = 0;
  int outNode = 0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  double nodeDist = 0.0;
  for(int loopA=0;loopA<totSegs;loopA++){
    // Get Current Segment Length
    segLength = segmentLength[loopA];
    // Get end nodes
    inNode = segmentInNode[loopA];
    outNode = segmentOutNode[loopA];
    // Get node spatial distance
    dx = nodeXcoord[outNode] - nodeXcoord[inNode];
    dy = nodeYcoord[outNode] - nodeYcoord[inNode];
    dz = nodeZcoord[outNode] - nodeZcoord[inNode];

    nodeDist = sqrt(dx*dx + dy*dy + dz*dz);
    if(fabs(nodeDist-segLength)>1.0e-8){
      inconsistencyFound = true;
   //   segmentLength[loopA] = nodeDist;
    }
  }
  if(inconsistencyFound){
    printf("WARNING: Inconsistency detected between segment length and end node distance.\n");
  //  printf("Changing the segment lengths.\n");
  }
}


// PERFORM DATA CHECKING
template <class T>
int checkForDoubleEntry(vector<T> vec){
  for(int loopA=0;loopA<vec.size();loopA++){
    for(int loopB=loopA+1;loopB<vec.size();loopB++){
      if(vec[loopA] == vec[loopB]){
        return loopA;
      }
    }
  }
  return -1;
}

// PERFORM DATA CHECKING
template <class T>
bool checkContains(T value, vector<T> vec){
  for(int loopA=0;loopA<vec.size();loopA++){
    if(vec[loopA] == value){
      return true;
    }
  }
  return false;
}

// CHECK FOR POSITIVE VALUES
int checkForPositiveVal(cvDoubleVec vec){
  for(int loopA=0;loopA<vec.size();loopA++){
    if(vec[loopA] < 0.0){
      return loopA;
    }
  }
  return -1;
}

// =====================
// PERFORM DATA CHECKING
// =====================
void cvOneDOptions::check(){

  int dblIDX = 0;
  int chkIDX = 0;

  // Check for double node name
  dblIDX = checkForDoubleEntry(nodeName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Node Name: " + nodeName[dblIDX] + "\n").c_str());
  }
  // Check for double joint name
  dblIDX = checkForDoubleEntry(jointName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Joint Name: " + jointName[dblIDX] + "\n").c_str());
  }
  // Check for double jointInlet name
  dblIDX = checkForDoubleEntry(jointInletListNames);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double JointInlet Name: " + jointInletListNames[dblIDX] + "\n").c_str());
  }
   // Check for double jointOutlet name
  dblIDX = checkForDoubleEntry(jointOutletListNames);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double JointOutlet Name: " + jointOutletListNames[dblIDX] + "\n").c_str());
  }
  // Check for double material name
  dblIDX = checkForDoubleEntry(materialName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Material Name: " + materialName[dblIDX] + "\n").c_str());
  }
  // Check for double data table name
  dblIDX = checkForDoubleEntry(dataTableName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Data Table Name: " + dataTableName[dblIDX] + "\n").c_str());
  }
  // Check for double segment Name
  dblIDX = checkForDoubleEntry(segmentName);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Segment Name: " + segmentName[dblIDX] + "\n").c_str());
  }
  // Check for double segment ID
  dblIDX = checkForDoubleEntry(segmentID);
  if(dblIDX >= 0){
    throw cvException(string("ERROR: Double Segment ID: " + to_string(segmentID[dblIDX]) + "\n").c_str());
  }
  // Check for negative area in input
  chkIDX = checkForPositiveVal(segmentInInletArea);
  if(chkIDX >= 0){
    throw cvException(string("ERROR: Negative Inlet area in segment : " + segmentName[chkIDX] + "\n").c_str());
  }
  chkIDX = checkForPositiveVal(segmentInOutletArea);
  if(chkIDX >= 0){
    throw cvException(string("ERROR: Negative Outlet area in segment : " + segmentName[chkIDX] + "\n").c_str());
  }

  // Check the consistency between node coords and segment lengths
  checkSegmentLengthConsistency();

  // Check if the segments refer to a node that is not there
  checkSegmentHasNodes();

  // Check if the joints refer to a node that is not there
  checkJointHasNodes();

}

void cvOneDOptions::checkSegmentHasNodes(){
  int inNode = 0;
  int outNode = 0;
  for(int loopA=0;loopA<segmentName.size();loopA++){
    // Get end nodes
    inNode = segmentInNode[loopA];
    outNode = segmentOutNode[loopA];
    // Check
    if((inNode < 0)||(inNode >= nodeName.size())){
      throw cvException(string("ERROR: Missing Node in Segment: " + segmentName[loopA] + "\n").c_str());
    }
    if((outNode < 0)||(outNode >= nodeName.size())){
      throw cvException(string("ERROR: Missing Node in Segment: " + segmentName[loopA] + "\n").c_str());
    }
  }
}

void cvOneDOptions::checkJointHasNodes(){
  string currNodeName;
  for(int loopA=0;loopA<jointName.size();loopA++){
    // Get end nodes
    currNodeName = jointNode[loopA];
    // Check If
    if(!checkContains(currNodeName,nodeName)){
      throw cvException(string("ERROR: Missing Node " + currNodeName + " in Joint: " + jointName[loopA] + "\n").c_str());
    }
  }
}
