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
#include "cvOneDGlobal.h"

namespace cvOneD{

// ======================
// READ SINGLE MODEL FILE
// ======================
void readOptionsLegacyFormat(string inputFile, cvOneDOptions* opts){

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
//          opts->jointNode.push_back(atof(tokenizedString[2].c_str()));
          opts->jointNode.push_back(tokenizedString[2]);
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
        if(opts->solverOptionDefined){
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
        opts->solverOptionDefined = true;
      }else if(upper_string(tokenizedString[0]) == std::string("OUTPUT")){
        if(tokenizedString.size() > 3){
          throw cvException(string("ERROR: Too many parameters for OUTPUT token. Line " + to_string(lineCount) + "\n").c_str());
        }else if(tokenizedString.size() < 2){
          throw cvException(string("ERROR: Not enough parameters for OUTPUT token. Line " + to_string(lineCount) + "\n").c_str());
        }
        // Output Type
        if(upper_string(tokenizedString[1]) == "TEXT"){
          cvOneDGlobal::outputType = OutputTypeScope::OUTPUT_TEXT;
        }else if(upper_string(tokenizedString[1]) == "VTK"){
          cvOneDGlobal::outputType = OutputTypeScope::OUTPUT_VTK;
        }else if(upper_string(tokenizedString[1]) == "BOTH"){
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
  // Close File
  infile.close();
}

} // namespace cvOneD