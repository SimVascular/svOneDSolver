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

#include "cvOneDOptionsLegacySerializer.h"

#include <iostream>
#include <fstream>
#include <string.h>

#include "cvOneDUtility.h"

namespace cvOneD{

// ======================
// READ SINGLE MODEL FILE
// ======================
void readOptionsLegacyFormat(string inputFile, options* opts){

  // Message
  printf("\n");
  printf("Reading file: %s ... \n",inputFile.c_str());

  // Declare
  cvStringVec tokenizedString;
  cvLongVec   tempIntVec;
  string      matType;
  cvDoubleVec temp;
  
  bool solverOptionsDefined = false;
  bool modelNameDefined = false;

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
    buffer = trim_string(buffer);

    // Tokenize String
    tokenizedString = split_string(buffer, " ,\t");
    if (tokenizedString.size() == 0) { 
      continue;
    }

    // Check for Empty buffer
    if(!buffer.empty()){
      // CHECK THE ELEMENT TYPE
      if(upper_string(tokenizedString[0]) == "MODEL"){
        //printf("Found Model.\n");
        if(modelNameDefined){
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
        modelNameDefined = true;
      }else if(upper_string(tokenizedString[0]) == "NODE"){
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
      }else if(upper_string(tokenizedString[0]) == "JOINT"){
        // printf("Found Joint.\n");
        if(tokenizedString.size() > 5){
          throw cvException(string("ERROR: Too many parameters for JOINT token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 5){
          throw cvException(string("ERROR: Not enough parameters for JOINT token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // Get Joint Name
          opts->jointName.push_back(tokenizedString[1]);

          // The legacy input files included a "joint node", but it is
          // never utilized in simulation. The legacy behavior was, instead
          // to simply associate each joint with the sequentially listed
          // nodeXcoord, nodeYcoord, and nodeZcoord.
          //
          // To preserve this legacy behavior, we're not going to deserialize
          // the jointNode. Instead, we're going to replace "jointNode" with
          // the sequential node item name. That's **wrong** -- we should be 
          // using "jointNode"!
          // But! it preserves the legacy behavior of the existing legacy 
          // options files when we refactor the actual utilization of the 
          // options to use the jointNode.
          //
          // The actual population of the jointNode option is done in a 
          // post-processing step after we finish processing the file.
          opts->jointInletName.push_back(tokenizedString[3]);
          opts->jointOutletName.push_back(tokenizedString[4]);

        }catch(...){
          throw cvException(string("ERROR: Invalid JOINT Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(upper_string(tokenizedString[0]) == "JOINTINLET"){
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
      }else if(upper_string(tokenizedString[0]) == "JOINTOUTLET"){
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
      }else if(upper_string(tokenizedString[0]) == "SEGMENT"){
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
      }else if(upper_string(tokenizedString[0]) == "SOLVEROPTIONS"){
        // printf("Found Solver Options.\n");
        if(solverOptionsDefined){
          throw cvException("ERROR: SOLVEROPTIONS already defined\n");
        }
        if(tokenizedString.size() > 10){
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
        }catch(...){
          throw cvException(string("ERROR: Invalid SOLVEROPTIONS Format. Line " + to_string(lineCount) + "\n").c_str());
        }
        solverOptionsDefined = true;
      }else if(upper_string(tokenizedString[0]) == std::string("COUPLING")){
        // printf("Found Coupling Options.\n");
        if(tokenizedString.size() > 4){
          throw cvException(string("ERROR: Too many parameters for COUPLING token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 4){
          throw cvException(string("ERROR: Not enough parameters for COUPLING token. Line " + to_string(lineCount) + "\n").c_str());
        }
        try{
          // string couplingStatus (ON/OFF)
          opts->couplingStatus = tokenizedString[1];
          // string couplingType (DIR/NEU)
          opts->couplingType = tokenizedString[2];
          // int couplingSubsteps
          opts->couplingSubsteps = atoi(tokenizedString[3].c_str());
        }catch(...){
          throw cvException(string("ERROR: Invalid COUPLING Format. Line " + to_string(lineCount) + "\n").c_str());
        }
      }else if(upper_string(tokenizedString[0]) == std::string("OUTPUT")){
        if(tokenizedString.size() > 3){
          throw cvException(string("ERROR: Too many parameters for OUTPUT token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 2){
          throw cvException(string("ERROR: Not enough parameters for OUTPUT token. Line " + to_string(lineCount) + "\n").c_str());
        }

        // Collect and store the output type data
        opts->outputType = tokenizedString[1];
        if(tokenizedString.size() > 2){
          opts->vtkOutputType = atoi(tokenizedString[2].c_str());
        }

      }else if(upper_string(tokenizedString[0]) == std::string("DATATABLE")){
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
            buffer = trim_string(buffer);
            // Tokenize String
            tokenizedString = split_string(buffer, " ,\t");
            if (tokenizedString.size() == 0) { 
              break;
            }
            // Check for Empty buffer
            if(!buffer.empty()){
              if(upper_string(tokenizedString[0]) == std::string("ENDDATATABLE")){
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
      }else if(upper_string(tokenizedString[0]) == std::string("MATERIAL")){
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
          if(upper_string(matType) == "OLUFSEN"){
            opts->materialParam1.push_back(atof(tokenizedString[7].c_str()));
            opts->materialParam2.push_back(atof(tokenizedString[8].c_str()));
            opts->materialParam3.push_back(atof(tokenizedString[9].c_str()));
          }else if(upper_string(matType) == "LINEAR"){
            opts->materialParam1.push_back(atof(tokenizedString[7].c_str()));
            opts->materialParam2.push_back(0.0);
            opts->materialParam3.push_back(0.0);
          }else{
            throw cvException(string("ERROR: Invalid MATERIAL Type. Line " + to_string(lineCount) + "\n").c_str());
          }
        }catch(...){
          throw cvException("ERROR: Invalid MATERIAL Format.\n");
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

  // Post-process: to preserve legacy behavior, we're
  // making each joint node the corresponding sequential 
  // node from the node array. See comment above for additional details.
  if(opts->jointName.size() > opts->nodeName.size()){
    throw cvException("ERROR: there must be at least as many nodes as there are joints.\n");
  }
  for(size_t k = 0; k < opts->jointName.size(); ++k){
    opts->jointNode.push_back(opts->nodeName.at(k));
  }

  // Close File
  infile.close();
}

namespace{

// PRINT MODEL NAME
void printModelName(options const& opts,FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"MODEL NAME: %s\n",opts.modelName.c_str());
}

// PRINT NODE DATA
void printNodeData(options const& opts,FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"NODE DATA\n");
  for(long int loopA=0;loopA<opts.nodeName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"NODE: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",opts.nodeName[loopA].c_str());
    fprintf(f,"X: %f\n",opts.nodeXcoord[loopA]);
    fprintf(f,"Y: %f\n",opts.nodeYcoord[loopA]);
    fprintf(f,"Z: %f\n",opts.nodeZcoord[loopA]);
  }
}

// PRINT JOINT DATA
void printJointData(options const& opts,FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"JOINT DATA\n");
  for(long int loopA=0;loopA<opts.jointName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"JOINT: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",opts.jointName[loopA].c_str());
    fprintf(f,"NODE: %s\n",opts.jointNode[loopA].c_str());
    fprintf(f,"X: %f\n",opts.nodeXcoord[loopA]);
    fprintf(f,"Y: %f\n",opts.nodeYcoord[loopA]);
    fprintf(f,"Z: %f\n",opts.nodeZcoord[loopA]);
    fprintf(f,"INLET NAME: %s\n",opts.jointInletName[loopA].c_str());
    fprintf(f,"OUTLET NAME: %s\n",opts.jointOutletName[loopA].c_str());
  }
}

// PRINT SEGMEMENT DATA
void printSegmentData(options const& opts,FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"SEGMENT DATA\n");
  for(long int loopA=0;loopA<opts.segmentName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"SEGMENT N.: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",opts.segmentName[loopA].c_str());
    fprintf(f,"ID: %ld\n",opts.segmentID[loopA]);
    fprintf(f,"LENGTH: %f\n",opts.segmentLength[loopA]);
    fprintf(f,"NUM ELEMENTS: %ld\n",opts.segmentTotEls[loopA]);
    fprintf(f,"IN NODE: %ld\n",opts.segmentInNode[loopA]);
    fprintf(f,"OUT NODE: %ld\n",opts.segmentOutNode[loopA]);
    fprintf(f,"IN AREA: %f\n",opts.segmentInInletArea[loopA]);
    fprintf(f,"OUT AREA: %f\n",opts.segmentInOutletArea[loopA]);
    fprintf(f,"IN FLOW: %f\n",opts.segmentInFlow[loopA]); //Take out
    fprintf(f,"MAT DATATABLE NAME: %s\n",opts.segmentMatName[loopA].c_str());
    fprintf(f,"LOSS TYPE: %s\n",opts.segmentLossType[loopA].c_str()); //take out
    fprintf(f,"BRANCH ANGLE: %f\n",opts.segmentBranchAngle[loopA]); //take out
    fprintf(f,"UPSTREAM SEGMENT: %ld\n",opts.segmentUpstreamSegment[loopA]); //take out
    fprintf(f,"BRANCH SEGMENT: %ld\n",opts.segmentBranchSegment[loopA]); //take out
    fprintf(f,"BOUNDARY CONDITIONS TYPE: %s\n",opts.segmentBoundType[loopA].c_str());
    fprintf(f,"INLET DATA TABLE: %s\n",opts.segmentDataTableName[loopA].c_str());
  }
}

// PRINT SOLVER OPTION DATA
void printSolverOptions(options const& opts,FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"PRINT SOLVER OPTION DATA\n");
  fprintf(f,"TIME STEP: %f\n",opts.timeStep);
  fprintf(f,"STEP SIZE: %ld\n",opts.stepSize);
  fprintf(f,"MAX STEP: %ld\n",opts.maxStep);
  fprintf(f,"QUADRATURE POINTS: %ld\n",opts.quadPoints);
  fprintf(f,"INLET DATA TABLE: %s\n",opts.inletDataTableName.c_str());
  fprintf(f,"BOUNDARY TYPE: %s\n",opts.boundaryType.c_str());
  fprintf(f,"CONVERGENCE TOLERANCE: %f\n",opts.convergenceTolerance);
  fprintf(f,"USE IV: %d\n",opts.useIV);
  fprintf(f,"USE STABILIZATION: %d\n",opts.useStab);
}

// PRINT MATERIAL DATA
void printMaterialData(options const& opts,FILE* f){
  fprintf(f,"--- \n");
  // MATERIAL
  for(long int loopA=0;loopA<opts.materialName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"MATERIAL NUM: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",opts.materialName[loopA].c_str());
    fprintf(f,"TYPE: %s\n",opts.materialType[loopA].c_str());
    fprintf(f,"DENSITY: %f\n",opts.materialDensity[loopA]);
    fprintf(f,"VISCOSITY: %f\n",opts.materialViscosity[loopA]);
    fprintf(f,"REF PRESSURE: %f\n",opts.materialPRef[loopA]);
    fprintf(f,"EXPONENT: %f\n",opts.materialExponent[loopA]);
    fprintf(f,"PARAM 1: %f\n",opts.materialParam1[loopA]);
    fprintf(f,"PARAM 2: %f\n",opts.materialParam2[loopA]);
    fprintf(f,"PARAM 3: %f\n",opts.materialParam3[loopA]);
  }
}

// PRINT JOINTINLET DATA
void printJointInletData(options const& opts,FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"JOINTINLET DATA\n");
  for(long int loopA=0;loopA<opts.jointInletListNames.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"JOINTINLET NUM: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",opts.jointInletListNames[loopA].c_str());
    fprintf(f,"TOTAL : %ld\n",opts.jointInletListNumber[loopA]);
    if(opts.jointInletListNumber[loopA] > 0){
      fprintf(f,"LIST ITEMS\n");
      for(size_t loopB=0;loopB<opts.jointInletList[loopA].size();loopB++){
        fprintf(f,"%ld \n",opts.jointInletList[loopA][loopB]);
      }
    }
    fprintf(f,"\n");
  }
}

// PRINT JOINTOPUTLET DATA
void printJointOutletData(options const& opts,FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"JOINTOUTLET DATA\n");
  for(long int loopA=0;loopA<opts.jointOutletListNames.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"JOINTOUTLET NUM: %ld\n",loopA);
    fprintf(f,"NAME: %s\n",opts.jointOutletListNames[loopA].c_str());
    fprintf(f,"TOTAL : %ld\n",opts.jointOutletListNumber[loopA]);
    if(opts.jointOutletListNumber[loopA] > 0){
      fprintf(f,"LIST ITEMS\n");
      for(size_t loopB=0;loopB<opts.jointOutletList[loopA].size();loopB++){
        fprintf(f,"%ld \n",opts.jointOutletList[loopA][loopB]);
      }
    }
    fprintf(f,"\n");
  }
}

// PRINT DATA TABLES
void printDataTables(options const& opts,FILE* f){
  fprintf(f,"--- \n");
  fprintf(f,"DATA TABLES\n");
  for(long int loopA=0;loopA<opts.dataTableName.size();loopA++){
    fprintf(f,"--- \n");
    fprintf(f,"DATA TABLE ID: %ld\n",loopA);
    fprintf(f,"DATA TABLE NAME: %s\n",opts.dataTableName[loopA].c_str());
    fprintf(f,"DATA TABLE TYPE: %s\n",opts.dataTableType[loopA].c_str());
    for(size_t loopB=0;loopB<opts.dataTableVals[loopA].size();loopB++){
      fprintf(f,"VALUE n. %ld: %f\n",loopB+1,opts.dataTableVals[loopA][loopB]);
    }
  }
}

} // namespace

void printToLegacyFile(options const& opts, string const& fileName){
  FILE* f;
  f = fopen(fileName.c_str(),"w");

  printf("\n");
  printf("Printing Model Echo ... \n");
  fflush(stdout);
  printModelName(opts,f);
  //printf("1 ... ");
  //fflush(stdout);
  printNodeData(opts,f);
  //printf("1 ... ");
  //fflush(stdout);
  printJointData(opts,f);
  //printf("2... ");
  //fflush(stdout);
  printSegmentData(opts,f);
  //printf("3 ... ");
  //fflush(stdout);
  printSolverOptions(opts,f);
  //printf("4 ... ");
  //fflush(stdout);
  printMaterialData(opts,f);
  //printf("5 ... ");
  //fflush(stdout);
  printJointInletData(opts,f);
  //printf("6 ... ");
  //fflush(stdout);
  printJointOutletData(opts,f);
  //printf("7 ... ");
  //fflush(stdout);
  printDataTables(opts,f);
  //printf("8 ... ");
  //fflush(f);
  fclose(f);
}

} // namespace cvOneD