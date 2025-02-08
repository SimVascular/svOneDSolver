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

#include "main.h"
#include <string.h>

#include "cvOneDOptionsLegacySerializer.h"

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
    found = upper_string(key) == upper_string(cvOneDGlobal::gDataTables[count]->getName());
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
    if(upper_string(opts->materialType[loopA]) == "OLUFSEN"){
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

    if(upper_string(curveName) != "NONE") {
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

  double* vals;
  int tot;

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
                                opts->useStab);
  if(solveError == CV_ERROR){
    throw cvException(string("ERROR: Error Solving Model\n").c_str());
  }
  delete [] inletCurveTime;
  delete [] inletCurveValue;
}

// ==============
// RUN ONEDSOLVER
// ==============
void runOneDSolver(string inputFile){

  // Create options from legacy format file
  cvOneDOptions* opts = new cvOneDOptions();
  cvOneD::readOptionsLegacyFormat(inputFile,opts);

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
    return 1;
  }
  printf("Completed!\n");
  return 0;

}
