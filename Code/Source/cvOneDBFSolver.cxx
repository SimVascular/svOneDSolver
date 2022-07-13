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

# include <time.h>

# include "cvOneDGlobal.h"
# include "cvOneDString.h"
# include "cvOneDException.h"
# include "cvOneDBFSolver.h"
# include "cvOneDMaterial.h"
# include "cvOneDMthSegmentModel.h"
# include "cvOneDMthBranchModel.h"

#ifndef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h> 


# define USE_SKYLINE

# ifdef USE_SKYLINE
  # include "cvOneDSkylineMatrix.h"
  # include "cvOneDSkylineLinearSolver.h"
# endif

# ifdef USE_SUPERLU
  # include "sparse/cvOneDSparseMatrix.h"
  # include "sparse/cvOneDSparseLinearSolver.h"
# endif

# ifdef USE_CSPARSE
  # include "sparse/cvOneDSparseMatrix.h"
  # include "sparse/cvOneDSparseLinearSolver.h"
# endif

# define baryeTommHg 0.0007500615613026439

//
//  cvOneDBFSolver.cpp - Source for a One-Dimensional Network Blood Flow Solver
//  ~~~~~~~~~~~~
//  This is a static class which abstracts the blood flow solver.
//  Essentially, this is a driver for 1D finite element blood flow solver
//  without any interfaces.
//

// Static Declarations...
bool                          cvOneDBFSolver::wasSet = false;
cvOneDModel*                  cvOneDBFSolver::model = NULL;
cvOneDFEAVector*              cvOneDBFSolver::currentSolution = NULL;
cvOneDFEAVector*              cvOneDBFSolver::previousSolution = NULL;
cvOneDFEAVector*              cvOneDBFSolver::increment = NULL;
cvOneDFEAVector*              cvOneDBFSolver::rhs = NULL;
cvOneDMatrix<double>          cvOneDBFSolver::TotalSolution;
cvOneDFEAMatrix*              cvOneDBFSolver::lhs = NULL;
vector<cvOneDMthModelBase*>   cvOneDBFSolver::mathModels;
vector<cvOneDSubdomain*>      cvOneDBFSolver::subdomainList;
vector<cvOneDFEAJoint*>       cvOneDBFSolver::jointList;
vector<int>                   cvOneDBFSolver::outletList;
double                        cvOneDBFSolver::currentTime = 0;
double                        cvOneDBFSolver::deltaTime = 0;
double                        cvOneDBFSolver::Period = 0;
long                          cvOneDBFSolver::stepSize = 0;
long                          cvOneDBFSolver::maxStep = 0;
long                          cvOneDBFSolver::quadPoints = 0;
double*                       cvOneDBFSolver::flowRate = NULL;
double*                       cvOneDBFSolver::flowTime = NULL;
long                          cvOneDBFSolver::numFlowPts = 0;
double                        cvOneDBFSolver::convCriteria = 0;
BoundCondType                 cvOneDBFSolver::inletBCtype;
int                           cvOneDBFSolver::ASCII = 1;

// SET MODE PTR
void cvOneDBFSolver::SetModelPtr(cvOneDModel *mdl){
  model = mdl;
}

// ==================
// WRITE TEXT RESULTS
// ==================
void cvOneDBFSolver::postprocess_Text(){
  int j;
  int fileIter;
  int elCount = 0;
  for(fileIter = 0; fileIter < model -> getNumberOfSegments(); fileIter++){

    cvOneDSegment *curSeg = model -> getSegment(fileIter);
    cvOneDMaterial* curMat = subdomainList[fileIter]->GetMaterial();

    long numEls = curSeg -> getNumElements();
    double segLength = curSeg->getSegmentLength();

    long startOut = elCount;
    long finishOut = elCount + ((numEls+1)*2);

    char *tmp1 = curSeg -> getSegmentName();

    char tmp2[512];
    char tmp3[512];
    char tmp4[512];
    char tmp5[512];
    char tmp6[512]; // WSS

    char *btemp= model-> getModelName(); // add to write out binary files for java

    strcpy(tmp2, btemp);
    strcpy(tmp3, btemp);
    strcpy(tmp4, btemp);
    strcpy(tmp5, btemp);
    strcpy(tmp6, btemp);

    strcat(tmp2, tmp1);
    strcat(tmp2, "_flow.dat");
    cout << tmp2 << endl;
    strcat(tmp3, tmp1);
    strcat(tmp3, "_area.dat");
    cout << tmp3 << endl;
    strcat(tmp4, tmp1);
    strcat(tmp4, "_pressure.dat");
    cout << tmp4 << endl;
    strcat(tmp5, tmp1);
    strcat(tmp5,"_Re.dat");
    cout << tmp5 << endl;
    strcat(tmp6, tmp1);
    strcat(tmp6,"_wss.dat");
    cout << tmp6 << endl;

    FILE *fp1,*fp2,*fp3,*fp5,*fp6; //for binary
    ofstream flow,area,pressure,reynolds,wss; // for ASCII

    if(ASCII){

        // Text Files
        flow.open(tmp2, ios::out);
        area.open(tmp3, ios::out);
        pressure.open(tmp4, ios::out);
        reynolds.open(tmp5, ios::out);
        wss.open(tmp6, ios::out);

        flow.precision(OUTPUT_PRECISION);
        area.precision(OUTPUT_PRECISION);
        pressure.precision(OUTPUT_PRECISION);
        reynolds.precision(OUTPUT_PRECISION);
        wss.precision(OUTPUT_PRECISION);

    }else{

      // Binary Files
      fp1 = fopen(tmp2,"wb");
      fp2 = fopen(tmp3,"wb");
      fp3 = fopen(tmp4,"wb");
      fp5 = fopen(tmp5,"wb");
      fp6 = fopen(tmp6,"wb");

    }
    double val;

    // Output the flow file
    for(j=startOut+1;j<finishOut;j+=2){
      for(int i=0;i<TotalSolution.Rows();i++){
        val = (double)TotalSolution[i][j] ;
        if(ASCII){
          flow << TotalSolution[i][j] << " ";
        }else{
          fwrite(&val,sizeof(double),1,fp1);
        }
      }
      if(ASCII) flow << endl;
    }

    // Output the Area/Pressure/Reynolds/WSS file
    int ii;
    double Re, r, flo, wssVal;

    for(ii=0,j=startOut;ii<numEls+1 && j<finishOut;ii++,j+=2){
      double z = (ii/double(numEls))*segLength;
      for(int i=0;i<TotalSolution.Rows();i++){
        //write area
        val = (double)TotalSolution[i][j];
        if(ASCII){
          area << TotalSolution[i][j] << " ";
        }else{
          fwrite(&val,sizeof(double),1,fp2);
        }
        r = sqrt(val/M_PI);
        flo = (double)TotalSolution[i][j+1];
        //usual Re=rho/mu*D*velocity=rho/mu*Q*sqrt(4/Pi/Area)
        Re = curMat->GetDensity()/curMat->GetDynamicViscosity()*flo/sqrt(val)*sqrt(4.0/M_PI);

        // Write pressure - Initial (CGS) Units
        val = (double) curMat->GetPressure(TotalSolution[i][j],z);
        if(ASCII){
          pressure << curMat->GetPressure(TotalSolution[i][j],z) << " ";
        }else{
          fwrite(&val, sizeof(double),1,fp3);
        }

        if(ASCII){
          reynolds << Re << " ";
        }else{
          fwrite(&Re,sizeof(double),1,fp5);
        }

        // Compute Wall Shear Stresses for Poiseuille flow
        wssVal = (4.0*curMat->GetDynamicViscosity()*flo)/(M_PI*r*r*r);
        if(ASCII){
          wss << wssVal << " ";
        }else{
          fwrite(&wssVal,sizeof(double),1,fp6);
        }
      }
      if(ASCII){
        area << endl;
        pressure << endl;
        reynolds << endl;
        wss << endl;
      }
    }

    elCount += 2*(numEls+1);

    if(ASCII){
      area.close();
      pressure.close();
      reynolds.close();
      flow.close();
      wss.close();
    }else{
      fclose(fp1); //flow.dat
      fclose(fp2); //area.dat
      fclose(fp3); //pressure.dat
      fclose(fp5); //Re.dat
      fclose(fp6); //wss.dat
    }
  }
}

// Get the index of a segment given its ID
int cvOneDBFSolver::getSegmentIndex(int segID){
  cvOneDSegment* curSeg;
  bool found = false;
  int count = 0;
  while((!found)&&(count<model->getNumberOfSegments())){
    curSeg = model->getSegment(count);
    found = (curSeg->getSegmentID() == segID);
    if(!found){
      count++;
    }
  }
  if(found){
    return count;
  }else{
    throw cvException("ERROR:Cannot Find Segment ID in cvOneDBFSolver::getSegmentIndex\n");
  }
}

// EVAL SEGMENT AXIS SYSTEM
void evalSegmentLocalAxis(double axis[][3]){
  // Get Reference Axis Y or Z
  double ref[3];
  double mod = 0.0;
  double yDist = sqrt((axis[0][0] - 0.0)*(axis[0][0] - 0.0) +
                      (axis[1][0] - 1.0)*(axis[1][0] - 1.0) +
                      (axis[2][0] - 0.0)*(axis[2][0] - 0.0));
  if(yDist < 1.0e-5){
    // Reference Axis Z
    ref[0] = 0.0;
    ref[1] = 0.0;
    ref[2] = 1.0;
  }else{
    // Reference Axis Y
    ref[0] = 0.0;
    ref[1] = 1.0;
    ref[2] = 0.0;
  }
  // Make the First Outer Product
  axis[0][1] =  axis[1][0]*ref[2] - ref[1]*axis[2][0];
  axis[1][1] = -axis[0][0]*ref[2] + ref[0]*axis[2][0];
  axis[2][1] =  axis[0][0]*ref[1] - ref[0]*axis[1][0];
  mod = sqrt(axis[0][1]*axis[0][1] + axis[1][1]*axis[1][1] + axis[2][1]*axis[2][1]);
  axis[0][1] /= mod;
  axis[1][1] /= mod;
  axis[2][1] /= mod;
  // Make the Second Outer Product
  axis[0][2] =  axis[1][0]*axis[2][1] - axis[2][0]*axis[1][1];
  axis[1][2] = -axis[0][0]*axis[2][1] + axis[2][0]*axis[0][1];
  axis[2][2] =  axis[0][0]*axis[1][1] - axis[1][0]*axis[0][1];
  mod = sqrt(axis[0][2]*axis[0][2] + axis[1][2]*axis[1][2] + axis[2][2]*axis[2][2]);
  axis[0][2] /= mod;
  axis[1][2] /= mod;
  axis[2][2] /= mod;
}

// =======================================
// WRITE 3D XML VTK RESULTS - ALL ONE FILE
// =======================================
void cvOneDBFSolver::postprocess_VTK_XML3D_ONEFILE(){

  // Set Constant Number of Subdivisions on the vessel circumference
  int circSubdiv = 20;

  int currSegID = 0;
  cvOneDSegment* currSeg = NULL;
  cvOneDJoint* currJoint = NULL;
  cvOneDMaterial* curMat = NULL;
  cvOneDNode* currNode = NULL;

  // Set and open VTK file
  char fileName[512];
  char* modelName= model->getModelName();
  strcpy(fileName, modelName);
  strcat(fileName, ".vtp");
  FILE* vtkFile;
  vtkFile=fopen(fileName,"w");

  // Write VTK XML Header
  fprintf(vtkFile,"<?xml version=\"1.0\"?>\n");
  fprintf(vtkFile,"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(vtkFile,"<PolyData>\n");

  // DEFINE INCIDENCE
  std::vector<double> segInlets(model->getNumberOfSegments());
  std::vector<double> segOutlets(model->getNumberOfSegments());
  for(int loopSegment=0;loopSegment<model->getNumberOfSegments();loopSegment++){
    segInlets[loopSegment] = -1;
    segOutlets[loopSegment] = -1;
  }

  // FORM INCIDENCE AND STORE COORDS
  cvDoubleMat nodeList;
  cvDoubleVec temp;
  long* segNodes;
  int inletNodeID = 0;
  int outletNodeID = 0;
  if(model->getNumberOfNodes() > 0){
    for(int loopNode = 0; loopNode < model->getNumberOfNodes(); loopNode++){
      temp.clear();
      // Get joint coordinate
      currNode = model->getNode(loopNode);
      temp.push_back(currNode->x);
      temp.push_back(currNode->y);
      temp.push_back(currNode->z);
      nodeList.push_back(temp);
    }
    // Loop on the segments
    for(int loopSegment = 0; loopSegment < model->getNumberOfSegments(); loopSegment++){
      // Get current segment
      currSeg = model->getSegment(loopSegment);
      // Get End Nodes
      segNodes = currSeg->getInOutJoints();
      // Mark Inlets
      inletNodeID = segNodes[0];
      outletNodeID = segNodes[1];
      // Assign Inlets and Outlets
      segInlets[loopSegment] = inletNodeID;
      segOutlets[loopSegment] = outletNodeID;
    }
  }else{
    // There are not Joints and a single segment
    currSeg = model->getSegment(0);
    // First Node
    temp.clear();
    temp.push_back(0.0);
    temp.push_back(0.0);
    temp.push_back(0.0);
    nodeList.push_back(temp);
    // Second Node
    temp.clear();
    temp.push_back(currSeg->getSegmentLength());
    temp.push_back(0.0);
    temp.push_back(0.0);
    nodeList.push_back(temp);
    // Store the connectivity
    segInlets[0] = 0;
    segOutlets[0] = 1;
  }

  for(int loopSegment=0;loopSegment<model->getNumberOfSegments();loopSegment++){
    if(segInlets[loopSegment] == -1){
      // CHECK INLETS
      printf("ERROR: INLET FOR SEGMENT %d\n",loopSegment);
    }
    if(segOutlets[loopSegment] == -1){
      // CHECK OUTLETS
      printf("ERROR: OUTLET FOR SEGMENT %d\n",loopSegment);
    }
  }

  // Every Segment is a Piece
  int totSegmentPoints = 0;
  int inletSegJoint = 0;
  int outletSegJoint = 0;
  double segVers[3][3];
  double mod = 0.0;
  double currCentre[3] = {0.0};
  double currIniArea = 0.0;
  double currIniRad = 0.0;
  double currTheta = 0.0;
  cvDoubleVec tmp;
  cvDoubleMat segNodeList;
  double lengthByNodes = 0.0;
  double lengthBySegment = 0.0;
  int startOut = 0;
  int finishOut = 0;
  double segLength = 0.0;
  double z = 0.0;
  double flo = 0.0;
  double area = 0.0;
  double radius = 0.0;
  double Re = 0.0;
  double wss = 0.0;
  double iniArea = 0.0;
  double newArea = 0.0;
  double radDisp = 0.0;
  long segOffset = 0;
  for(int loopSegment=0;loopSegment<model->getNumberOfSegments();loopSegment++){

    // Clear the node list for the segment
    segNodeList.clear();

    // Get Current Segment
    currSeg = model->getSegment(loopSegment);

    // Compute the total number of points for this segment
    totSegmentPoints = (currSeg->getNumElements()+1) * circSubdiv;

    // Set the range for the totalsoluton of this segment
    startOut = segOffset;
    finishOut = segOffset + 2*(currSeg->getNumElements()+1);

    // Write Piece Header
    fprintf(vtkFile,"<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"%ld\" NumberOfPolys=\"0\">\n",totSegmentPoints,currSeg->getNumElements());

    // Get inlet and outlet joints
    inletSegJoint = segInlets[loopSegment];
    outletSegJoint = segOutlets[loopSegment];

    // Compute Segment Versor
    segVers[0][0] = nodeList[outletSegJoint][0] - nodeList[inletSegJoint][0];
    segVers[1][0] = nodeList[outletSegJoint][1] - nodeList[inletSegJoint][1];
    segVers[2][0] = nodeList[outletSegJoint][2] - nodeList[inletSegJoint][2];
    mod = sqrt(segVers[0][0]*segVers[0][0] + segVers[1][0]*segVers[1][0] + segVers[2][0]*segVers[2][0]);
    segVers[0][0] /= mod;
    segVers[1][0] /= mod;
    segVers[2][0] /= mod;

    lengthByNodes = mod;
    lengthBySegment = currSeg->getSegmentLength();

    // Compute Segment Local axis system
    evalSegmentLocalAxis(segVers);

    // Loop on the number of elements
    for(int loopEl=0;loopEl<currSeg->getNumElements() + 1;loopEl++){
      // Compress/Elongate solution by length between nodes rather than defined segment length - to prioritize segment length
      // and maintain similar geometry to node definitions, MD 4/2/19
      currCentre[0] = nodeList[inletSegJoint][0] + loopEl*lengthByNodes/double(currSeg->getNumElements())*segVers[0][0];
      currCentre[1] = nodeList[inletSegJoint][1] + loopEl*lengthByNodes/double(currSeg->getNumElements())*segVers[1][0];
      currCentre[2] = nodeList[inletSegJoint][2] + loopEl*lengthByNodes/double(currSeg->getNumElements())*segVers[2][0];

      // Get initial radius at current location
      currIniArea = currSeg->getInitInletS() + (loopEl/double(currSeg->getNumElements()))*(currSeg->getInitOutletS() - currSeg->getInitInletS());
      currIniRad = sqrt(currIniArea/M_PI);

      // Loop on the subdivisions
      for(int loopSubdiv=0;loopSubdiv<circSubdiv;loopSubdiv++){
        currTheta = loopSubdiv*2*M_PI/double(circSubdiv);
        tmp.clear();
        tmp.push_back(currCentre[0] + currIniRad*segVers[0][1]*cos(currTheta) + currIniRad*segVers[0][2]*sin(currTheta));
        tmp.push_back(currCentre[1] + currIniRad*segVers[1][1]*cos(currTheta) + currIniRad*segVers[1][2]*sin(currTheta));
        tmp.push_back(currCentre[2] + currIniRad*segVers[2][1]*cos(currTheta) + currIniRad*segVers[2][2]*sin(currTheta));
        segNodeList.push_back(tmp);
      }
    }

    // List of Node Coordinates ready for export
    fprintf(vtkFile,"<Points>\n");
    fprintf(vtkFile,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(int loopA=0;loopA<segNodeList.size();loopA++){
      fprintf(vtkFile,"%e %e %e\n",segNodeList[loopA][0],segNodeList[loopA][1],segNodeList[loopA][2]);
    }
    fprintf(vtkFile,"</DataArray>\n");
    fprintf(vtkFile,"</Points>\n");

    // Write Strip Incidence and offset
    fprintf(vtkFile,"<Strips>\n");
    // Strip Connectivity
    fprintf(vtkFile,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int loopA=0;loopA<currSeg->getNumElements();loopA++){
      for(int loopB=0;loopB<circSubdiv;loopB++){
        fprintf(vtkFile,"%d %d ",loopA*circSubdiv+loopB,loopA*circSubdiv+loopB+circSubdiv);
      }
      fprintf(vtkFile,"%d %d ",loopA*circSubdiv+0,loopA*circSubdiv+0+circSubdiv);
      fprintf(vtkFile,"\n");
    }
    fprintf(vtkFile,"</DataArray>\n");
    // Strip Offset
    fprintf(vtkFile,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for(int loopA=0;loopA<currSeg->getNumElements();loopA++){
      fprintf(vtkFile,"%d ",(loopA+1)*(circSubdiv*2+2));
    }
    fprintf(vtkFile,"\n");
    fprintf(vtkFile,"</DataArray>\n");
    fprintf(vtkFile,"</Strips>\n");


    // PRINT OUTPUTS
    fprintf(vtkFile,"<PointData Scalars=\"ScalOutputs\" Vectors=\"VecOutputs\">\n");

    for(int loopTime=0;loopTime<TotalSolution.Rows();loopTime++){

      // PRINT FLOW RATES
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Flowrate_INCR_%05ld_TIME_%.5f\" NumberOfComponents=\"1\" format=\"ascii\">\n",loopTime*stepSize,loopTime*deltaTime*stepSize);
      for(int j=startOut+1;j<finishOut;j+=2){
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",(double)TotalSolution[loopTime][j]);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT AREA
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Area_INCR_%05ld_TIME_%.5f\" NumberOfComponents=\"1\" format=\"ascii\">\n",loopTime*stepSize,loopTime*deltaTime*stepSize);
      for(int j=startOut;j<finishOut;j+=2){
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",(double)TotalSolution[loopTime][j]);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT RADIAL DISPLACEMENTS AS VECTORS
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Disps_INCR_%05ld_TIME_%.5f\" NumberOfComponents=\"3\" format=\"ascii\">\n",loopTime*stepSize,loopTime*deltaTime*stepSize);
      for(int j=startOut;j<finishOut;j+=2){

        // Evaluate Initial Area at current location
        iniArea = currSeg->getInitInletS() + (((j-startOut)/2)/double(currSeg->getNumElements()))*(currSeg->getInitOutletS() - currSeg->getInitInletS());
        // Eval Current Area at current location
        newArea = TotalSolution[loopTime][j];
        // Evaluate Radial displacement
        radDisp = sqrt(newArea/M_PI) - sqrt(iniArea/M_PI);
        for(int k=0;k<circSubdiv;k++){
          // Print the three components for every point
          currTheta = k*2*M_PI/double(circSubdiv);
          tmp.clear();
          tmp.push_back(radDisp*segVers[0][1]*cos(currTheta) + radDisp*segVers[0][2]*sin(currTheta));
          tmp.push_back(radDisp*segVers[1][1]*cos(currTheta) + radDisp*segVers[1][2]*sin(currTheta));
          tmp.push_back(radDisp*segVers[2][1]*cos(currTheta) + radDisp*segVers[2][2]*sin(currTheta));
          // Write values
          fprintf(vtkFile,"%e %e %e ",tmp[0],tmp[1],tmp[2]);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT PRESSURE IN MMHG
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Pressure_mmHg_INCR_%05ld_TIME_%.5f\" NumberOfComponents=\"1\" format=\"ascii\">\n",loopTime*stepSize,loopTime*deltaTime*stepSize);
      segLength = currSeg->getSegmentLength();
      curMat = subdomainList[loopSegment]->GetMaterial();
      int section = 0;
      for(int j=startOut;j<finishOut;j+=2){
        z = (section/(double)currSeg->getNumElements())*segLength;
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",curMat->GetPressure(TotalSolution[loopTime][j],z)*baryeTommHg);
        }
        fprintf(vtkFile,"\n");
        section++;
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT REYNOLDS NUMBER
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Reynolds_INCR_%05ld_TIME_%.5f\" NumberOfComponents=\"1\" format=\"ascii\">\n",loopTime*stepSize,loopTime*deltaTime*stepSize);
      curMat = subdomainList[loopSegment]->GetMaterial();
      for(int j=startOut;j<finishOut;j+=2){
        // Get Flow
        flo = (double)TotalSolution[loopTime][j+1];
        area = (double)TotalSolution[loopTime][j];
        // Get Area
        Re = curMat->GetDensity()/curMat->GetDynamicViscosity()*flo/sqrt(area)*sqrt(4.0/M_PI);
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",Re);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT WSS
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"WSS_INCR_%05ld_TIME_%.5f\" NumberOfComponents=\"1\" format=\"ascii\">\n",loopTime*stepSize,loopTime*deltaTime*stepSize);
      curMat = subdomainList[loopSegment]->GetMaterial();
      for(int j=startOut;j<finishOut;j+=2){
        // Get Flow
        flo = (double)TotalSolution[loopTime][j+1];
        radius = sqrt((double)TotalSolution[loopTime][j]/M_PI);
        // Get Area
        wss = 4.0*curMat->GetDynamicViscosity()*flo/(M_PI*radius*radius*radius);
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",wss);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

    }
    fprintf(vtkFile,"</PointData>\n");
    // Close Piece
    fprintf(vtkFile,"</Piece>\n");

    // Increment Segment Offset
    segOffset += 2*(currSeg->getNumElements()+1);

  } // End Segment Loop

  // Close
  fprintf(vtkFile,"</PolyData>\n");
  fprintf(vtkFile,"</VTKFile>\n");
  printf("Results Exported to VTK File: %s\n",fileName);
}

// =======================================================
// WRITE 3D XML VTK RESULTS - MULTIPLE FILES FOR ANIMATION
// =======================================================
void cvOneDBFSolver::postprocess_VTK_XML3D_MULTIPLEFILES(){

  // Set Constant Number of Subdivisions on the vessel circumference
  int circSubdiv = 20;

  int currSegID = 0;
  cvOneDSegment* currSeg = NULL;
  cvOneDJoint* currJoint = NULL;
  cvOneDMaterial* curMat = NULL;
  cvOneDNode* currNode = NULL;

  // DEFINE INCIDENCE
  std::vector<double> segInlets(model->getNumberOfSegments());
  std::vector<double> segOutlets(model->getNumberOfSegments());
  for(int loopSegment=0;loopSegment<model->getNumberOfSegments();loopSegment++){
    segInlets[loopSegment] = -1;
    segOutlets[loopSegment] = -1;
  }

  // FORM INCIDENCE AND STORE COORDS
  cvDoubleMat nodeList;
  cvDoubleVec temp;
  long* segNodes;
  int inletNodeID = 0;
  int outletNodeID = 0;
  if(model->getNumberOfNodes() > 0){
    for(int loopNode = 0; loopNode < model->getNumberOfNodes(); loopNode++){
      temp.clear();
      // Get joint coordinate
      currNode = model->getNode(loopNode);
      temp.push_back(currNode->x);
      temp.push_back(currNode->y);
      temp.push_back(currNode->z);
      nodeList.push_back(temp);
    }
    // Loop on the segments
    for(int loopSegment = 0; loopSegment < model->getNumberOfSegments(); loopSegment++){
      // Get current segment
      currSeg = model->getSegment(loopSegment);
      // Get End Nodes
      segNodes = currSeg->getInOutJoints();
      // Mark Inlets
      inletNodeID = segNodes[0];
      outletNodeID = segNodes[1];
      // Assign Inlets and Outlets
      segInlets[loopSegment] = inletNodeID;
      segOutlets[loopSegment] = outletNodeID;
    }
  }else{
    // There are not Joints and a single segment
    currSeg = model->getSegment(0);
    // First Node
    temp.clear();
    temp.push_back(0.0);
    temp.push_back(0.0);
    temp.push_back(0.0);
    nodeList.push_back(temp);
    // Second Node
    temp.clear();
    temp.push_back(currSeg->getSegmentLength());
    temp.push_back(0.0);
    temp.push_back(0.0);
    nodeList.push_back(temp);
    // Store the connectivity
    segInlets[0] = 0;
    segOutlets[0] = 1;
  }

  for(int loopSegment=0;loopSegment<model->getNumberOfSegments();loopSegment++){
    if(segInlets[loopSegment] == -1){
      // CHECK INLETS
      printf("ERROR: INLET FOR SEGMENT %d\n",loopSegment);
    }
    if(segOutlets[loopSegment] == -1){
      // CHECK OUTLETS
      printf("ERROR: OUTLET FOR SEGMENT %d\n",loopSegment);
    }
  }

  // LOOP OVER TIME
  int totSegmentPoints = 0;
  int totSegmentSolutions = 0;
  long segOffset = 0;
  int inletSegJoint = 0;
  int outletSegJoint = 0;
  double segVers[3][3];
  double mod = 0.0;
  double currCentre[3] = {0.0};
  double currIniArea = 0.0;
  double currIniRad = 0.0;
  double currTheta = 0.0;
  cvDoubleVec tmp;
  cvDoubleMat segNodeList;
  double lengthByNodes = 0.0;
  double lengthBySegment = 0.0;
  int startOut = 0;
  int finishOut = 0;
  double segLength = 0.0;
  double z = 0.0;
  double flo = 0.0;
  double area = 0.0;
  double radius = 0.0;
  double Re = 0.0;
  double wss = 0.0;
  double iniArea = 0.0;
  double newArea = 0.0;
  double radDisp = 0.0;
  cvStringVec fileList;
  string fileName;

  for(int loopTime=0;loopTime<TotalSolution.Rows();loopTime++){

    // Set and open VTK file for current time step
    fileName = model->getModelName();
    char timeString[512];
    char suffix[512];
    sprintf(timeString, "_%05d", loopTime);
    fileName = fileName + string(timeString) + ".vtp";
    FILE* vtkFile;
    vtkFile = fopen(fileName.c_str(),"w");
    // Add to a list of files
    fileList.push_back(fileName);

    // Write VTK XML Header
    fprintf(vtkFile,"<?xml version=\"1.0\"?>\n");
    fprintf(vtkFile,"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(vtkFile,"<PolyData>\n");

    // Init Segment Offset
    segOffset = 0;

    // LOOP OVER THE SEGMENTS
    for(int loopSegment=0;loopSegment<model->getNumberOfSegments();loopSegment++){

      // Clear the node list for the segment
      segNodeList.clear();

      // Get Current Segment
      currSeg = model->getSegment(loopSegment);

      // Compute the total number of points for this segment
      totSegmentSolutions = (currSeg->getNumElements()+1);
      totSegmentPoints = totSegmentSolutions * circSubdiv;

      // Set the range for the totalsoluton of this segment
      startOut = segOffset;
      finishOut = segOffset + 2*(totSegmentSolutions);

      // Get Material
      curMat = subdomainList[loopSegment]->GetMaterial();

      // Write Piece Header
      fprintf(vtkFile,"<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"%ld\" NumberOfPolys=\"0\">\n",totSegmentPoints,currSeg->getNumElements());

      // Get inlet and outlet joints
      inletSegJoint = segInlets[loopSegment];
      outletSegJoint = segOutlets[loopSegment];

      // Compute Segment Versor
      segVers[0][0] = nodeList[outletSegJoint][0] - nodeList[inletSegJoint][0];
      segVers[1][0] = nodeList[outletSegJoint][1] - nodeList[inletSegJoint][1];
      segVers[2][0] = nodeList[outletSegJoint][2] - nodeList[inletSegJoint][2];
      mod = sqrt(segVers[0][0]*segVers[0][0] + segVers[1][0]*segVers[1][0] + segVers[2][0]*segVers[2][0]);
      segVers[0][0] /= mod;
      segVers[1][0] /= mod;
      segVers[2][0] /= mod;

      lengthByNodes = mod;
      lengthBySegment = currSeg->getSegmentLength();

      // Compute Segment Local axis system
      evalSegmentLocalAxis(segVers);

      // Loop on the number of elements
      for(int loopEl=0;loopEl<currSeg->getNumElements() + 1;loopEl++){
        // Compress/Elongate solution by length between nodes rather than defined segment length - to prioritize segment length
        // and maintain similar geometry to node definitions, MD 4/2/19
        currCentre[0] = nodeList[inletSegJoint][0] + loopEl*lengthByNodes/double(currSeg->getNumElements())*segVers[0][0];
        currCentre[1] = nodeList[inletSegJoint][1] + loopEl*lengthByNodes/double(currSeg->getNumElements())*segVers[1][0];
        currCentre[2] = nodeList[inletSegJoint][2] + loopEl*lengthByNodes/double(currSeg->getNumElements())*segVers[2][0];

        // Get initial radius at current location
        currIniArea = currSeg->getInitInletS() + (loopEl/double(currSeg->getNumElements()))*(currSeg->getInitOutletS() - currSeg->getInitInletS());
        currIniRad = sqrt(currIniArea/M_PI);

        // Loop on the subdivisions
        for(int loopSubdiv=0;loopSubdiv<circSubdiv;loopSubdiv++){
          currTheta = loopSubdiv*2*M_PI/double(circSubdiv);
          tmp.clear();
          tmp.push_back(currCentre[0] + currIniRad*segVers[0][1]*cos(currTheta) + currIniRad*segVers[0][2]*sin(currTheta));
          tmp.push_back(currCentre[1] + currIniRad*segVers[1][1]*cos(currTheta) + currIniRad*segVers[1][2]*sin(currTheta));
          tmp.push_back(currCentre[2] + currIniRad*segVers[2][1]*cos(currTheta) + currIniRad*segVers[2][2]*sin(currTheta));
          segNodeList.push_back(tmp);
        }
      }

      // List of Node Coordinates ready for export
      fprintf(vtkFile,"<Points>\n");
      fprintf(vtkFile,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
      for(int loopA=0;loopA<segNodeList.size();loopA++){
        fprintf(vtkFile,"%e %e %e\n",segNodeList[loopA][0],segNodeList[loopA][1],segNodeList[loopA][2]);
      }
      fprintf(vtkFile,"</DataArray>\n");
      fprintf(vtkFile,"</Points>\n");

      // Write Strip Incidence and offset
      fprintf(vtkFile,"<Strips>\n");
      // Strip Connectivity
      fprintf(vtkFile,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
      for(int loopA=0;loopA<currSeg->getNumElements();loopA++){
        for(int loopB=0;loopB<circSubdiv;loopB++){
          fprintf(vtkFile,"%d %d ",loopA*circSubdiv+loopB,loopA*circSubdiv+loopB+circSubdiv);
        }
        fprintf(vtkFile,"%d %d ",loopA*circSubdiv+0,loopA*circSubdiv+0+circSubdiv);
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");
      // Strip Offset
      fprintf(vtkFile,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
      for(int loopA=0;loopA<currSeg->getNumElements();loopA++){
        fprintf(vtkFile,"%d ",(loopA+1)*(circSubdiv*2+2));
      }
      fprintf(vtkFile,"\n");
      fprintf(vtkFile,"</DataArray>\n");
      fprintf(vtkFile,"</Strips>\n");


      // PRINT OUTPUTS
      fprintf(vtkFile,"<PointData Scalars=\"ScalOutputs\" Vectors=\"VecOutputs\">\n");

      // PRINT FLOW RATES
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Flowrate\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      for(int j=startOut+1;j<finishOut;j+=2){
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",(double)TotalSolution[loopTime][j]);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT AREA
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Area\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      for(int j=startOut;j<finishOut;j+=2){
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",(double)TotalSolution[loopTime][j]);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT RADIAL DISPLACEMENTS AS VECTORS
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Disps\" NumberOfComponents=\"3\" format=\"ascii\">\n");
      for(int j=startOut;j<finishOut;j+=2){

        // Evaluate Initial Area at current location
        iniArea = currSeg->getInitInletS() + (((j-startOut)/2)/double(currSeg->getNumElements()))*(currSeg->getInitOutletS() - currSeg->getInitInletS());
        // Eval Current Area at current location
        newArea = TotalSolution[loopTime][j];
        // Evaluate Radial displacement
        radDisp = sqrt(newArea/M_PI) - sqrt(iniArea/M_PI);
        for(int k=0;k<circSubdiv;k++){
          // Print the three components for every point
          currTheta = k*2*M_PI/double(circSubdiv);
          tmp.clear();
          tmp.push_back(radDisp*segVers[0][1]*cos(currTheta) + radDisp*segVers[0][2]*sin(currTheta));
          tmp.push_back(radDisp*segVers[1][1]*cos(currTheta) + radDisp*segVers[1][2]*sin(currTheta));
          tmp.push_back(radDisp*segVers[2][1]*cos(currTheta) + radDisp*segVers[2][2]*sin(currTheta));
          // Write values
          fprintf(vtkFile,"%e %e %e ",tmp[0],tmp[1],tmp[2]);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT PRESSURE IN MMHG
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Pressure_mmHg\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      segLength = currSeg->getSegmentLength();
      int section = 0;
      for(int j=startOut;j<finishOut;j+=2){
        z = (section/(double)currSeg->getNumElements())*segLength;
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",curMat->GetPressure(TotalSolution[loopTime][j],z)*baryeTommHg);
        }
        fprintf(vtkFile,"\n");
        section++;
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT REYNOLDS NUMBER
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"Reynolds\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      for(int j=startOut;j<finishOut;j+=2){
        // Get Flow
        flo = (double)TotalSolution[loopTime][j+1];
        area = (double)TotalSolution[loopTime][j];
        // Get Area
        Re = curMat->GetDensity()/curMat->GetDynamicViscosity()*flo/sqrt(area)*sqrt(4.0/M_PI);
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",Re);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // PRINT WSS
      fprintf(vtkFile,"<DataArray type=\"Float32\" Name=\"WSS\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      for(int j=startOut;j<finishOut;j+=2){
        // Get Flow
        flo = (double)TotalSolution[loopTime][j+1];
        // Get Radius
        radius = sqrt((double)TotalSolution[loopTime][j]/M_PI);
        // Get WSS
        wss = 4.0*curMat->GetDynamicViscosity()*flo/(M_PI*radius*radius*radius);
        for(int k=0;k<circSubdiv;k++){
          fprintf(vtkFile,"%e ",wss);
        }
        fprintf(vtkFile,"\n");
      }
      fprintf(vtkFile,"</DataArray>\n");

      // Close Pointdata
      fprintf(vtkFile,"</PointData>\n");
      // Close Piece
      fprintf(vtkFile,"</Piece>\n");

      // Increment Segment Offset
      segOffset += 2*(totSegmentSolutions);

    } // End Segment Loop
    // Close
    fprintf(vtkFile,"</PolyData>\n");
    fprintf(vtkFile,"</VTKFile>\n");
    fclose(vtkFile);
  }

  // Create Final PVD File with summary
  fileName = string(model->getModelName()) + string(".pvd");
  FILE* pvdFile;
  pvdFile = fopen(fileName.c_str(),"w");

  fprintf(pvdFile,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">");
  fprintf(pvdFile,"<Collection>");
  // Loop through the dataset
  double currTime = 0.0;
  for(int loopA=0;loopA<fileList.size();loopA++){
    fprintf(pvdFile,"<DataSet timestep=\"%e\" group=\"\" part=\"0\" file=\"%s\"/>",currTime,fileList[loopA].c_str());
    currTime += deltaTime * stepSize;
  }
  fprintf(pvdFile,"</Collection>");
  fprintf(pvdFile,"</VTKFile>");
  fclose(pvdFile);
  printf("Results Exported to VTK.\n");
}

// ====================
// MAIN SOLUTION DRIVER
// ====================
void cvOneDBFSolver::Solve(void){
  char errStr[256];
  // First check to make sure we've set a model pointer
  // Prior to solution attempt.
  if(model == NULL){
    cvOneDError::setErrorNumber(ErrorTypeScope::BAD_VALUE);
    strcpy(errStr,"In BFSolver::Solve(...), No model pointer was set prior to solution attempt.  Bailing out to avoid a coredump!");
    cvOneDError::setErrorString(errStr);
    cvOneDError::CallErrorHandler();
    exit(0);
  }

  // Query the model for information, this is where
  // the subdomain and material information get passed.
  QuerryModelInformation();

  DefineMthModels();

  CreateGlobalArrays();

  long id;
  id = model -> getNumberOfSegments();

  int i;
  for (i=0; i<id; i++){
    CalcInitProps(i);
  }

  // Start Solving the system.
  GenerateSolution();

  // Some Post Processing
  if(cvOneDGlobal::outputType == OutputTypeScope::OUTPUT_TEXT){
    postprocess_Text();
  }else if(cvOneDGlobal::outputType == OutputTypeScope::OUTPUT_VTK){
    if(cvOneDGlobal::vtkOutputType == 0){
      // Export in multifile format
      postprocess_VTK_XML3D_MULTIPLEFILES();
    }else{
      // All results in a single VTK File
      postprocess_VTK_XML3D_ONEFILE();
    }
  }else if(cvOneDGlobal::outputType == OutputTypeScope::OUTPUT_BOTH){
    postprocess_Text();
    if(cvOneDGlobal::vtkOutputType == 0){
      // Export in multifile format
      postprocess_VTK_XML3D_MULTIPLEFILES();
    }else{
      // All results in a single VTK File
      postprocess_VTK_XML3D_ONEFILE();
    }
  }
}

void cvOneDBFSolver::DefineInletFlow(double* time, double* flrt, int num){
  flowTime = new double[num];
  flowRate = new double[num];
  for (int i=0;i<num;i++){
    flowTime[i] = time[i];
    flowRate[i] = flrt[i];
  }
  numFlowPts = num;
  Period = flowTime[num-1];
}

void cvOneDBFSolver::DefineMthModels(){
  mathModels.clear();

  cout << "Subdomain No. "<<subdomainList.size() << endl;
  cout << "Joint No. "<< jointList.size() << endl;
  cout << "Outlet No. "<< outletList.size() << endl;
  cvOneDMthSegmentModel* segM = new cvOneDMthSegmentModel(subdomainList, jointList, outletList, quadPoints);

  //specify inlet flow rate boundary condition with time
  segM->SetInflowRate(flowTime, flowRate, numFlowPts, flowTime[numFlowPts-1]);
  Period = flowTime[numFlowPts-1];

  cvOneDMthBranchModel* branchM = new cvOneDMthBranchModel(subdomainList, jointList, outletList);
  AddOneModel(segM);
  AddOneModel(branchM);
}

void cvOneDBFSolver::AddOneModel(cvOneDMthModelBase* model){
  mathModels.push_back(model);
}

void cvOneDBFSolver::QuerryModelInformation(void)
{
    // place to create the subdomain and material.
    long is, ij;
    is = model -> getNumberOfSegments();
    ij = model -> getNumberOfJoints();

  jointList.resize(0);
  outletList.resize(0);
  subdomainList.resize(0);

    printf("\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Number of Joints: %ld\n",ij);
    printf("Number of Segments: %ld\n",is);
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("\n");

    long int i, j;

    // Loop on Joints
    for(i=0; i<ij; i++){
      cvOneDFEAJoint *feajoint = new cvOneDFEAJoint();
      cvOneDJoint *joint = model->getJoint(i);
      for(j=0;j < model -> getJoint(i)->InletSegments.size(); j++){
        feajoint->AddInletSubdomains(joint->InletSegments[j]);
      }
      for(j=0;j < joint->OutletSegments.size(); j++){
        feajoint->AddOutletSubdomains(joint->OutletSegments[j]);
      }
      jointList.push_back(feajoint);
    }

    // Loop on Segments
    long temp = 0;
    for (i=0; i<is; i++){
      cvOneDSegment* seg = model->getSegment(i);
      long nels = seg->getNumElements();
      double segLen = seg->getSegmentLength();
      int matID = seg->getMaterialID();
      MeshType mType = seg->getMeshType();
      double zin = seg->getInletZ();
      double zout = seg->getOutletZ();

      cvOneDSubdomain* subdomain = new cvOneDSubdomain;
      assert(subdomain != 0);
      subdomain -> SetNumberOfNodes(nels+1);
      subdomain -> SetNumberOfElements(nels);
      subdomain -> SetMeshType(mType);
      subdomain -> Init(zin, zout);


      // Get the Initial Properties of the subdomain...
      double Qo = 0.0;
      double P0 = 0.0;
      double dQ0_dT = 0.0;
      if(i == 0) {
        Qo = seg->getInitialFlow();
        P0 =  seg->getInitialPressure();
        dQ0_dT = 0.0;
      }
      double So = seg->getInitInletS();
      double Sn = seg->getInitOutletS();
      BoundCondType boundT = seg -> getBoundCondition();
      double  boundV= seg -> getBoundValue();

      // Set these in the subdomain.
      subdomain->SetInitialFlow(Qo);
      subdomain->SetInitialdFlowdT(dQ0_dT);
      subdomain->SetInitInletS(So);
      subdomain->SetInitialPressure(P0);
      subdomain->SetInitOutletS(Sn);
      subdomain->SetGlobal1stNodeID(temp);
      subdomain->SetBoundCondition(boundT);
      if(!seg->IsOutlet){
        subdomain->SetBoundCondition(BoundCondTypeScope::NOBOUND);
      }
      subdomain->SetupMaterial(matID);
      subdomain->GetMaterial()->SetPeriod(Period);

      // Set up Minor Loss
      subdomain->SetMinorLossType(seg->GetMinorLossType());
      if(seg->GetMinorLossType() != MinorLossScope::NONE){
        subdomain->SetBranchAngle(seg->GetBranchAngle());
        subdomain->SetUpstreamSeg(seg->GetUpstreamSeg());
        subdomain->SetBranchSeg(seg->GetBranchSeg());
      }

      // Set up boundary condition
      if (boundT == BoundCondTypeScope::RESISTANCE_TIME){
        double* time;
        double* resist;
        int num; // x-fer info from Segment to subdomain
        seg->getBoundResistanceValues(&resist,&time,&num);
        subdomain->SetBoundResistanceWave(time, resist, num);
      }else if(boundT == BoundCondTypeScope::RCR){
        double* rcr;
        int num;
        seg->getBoundRCRValues(&rcr,&num);
        subdomain -> SetBoundRCRValues(rcr,num);
        cout<<"RCR boundary condition"<<endl;

      }else if(boundT == BoundCondTypeScope::RESISTANCE){ // modified resistance taking Pd wgyang
        double* resistance_pd;
        int num;
        seg->getBoundRCRValues(&resistance_pd,&num);
        subdomain -> SetBoundResistPdValues(resistance_pd,num);
        cout<<"RESISTANCE boundary condition"<<endl;

      }else if (boundT == BoundCondTypeScope::PRESSURE_WAVE){
        double* time;
        double* pres;
        int num;
        seg->getBoundPressureValues(&pres,&time,&num);
        subdomain ->SetBoundPresWave(time, pres, num);

      }else if(boundT == BoundCondTypeScope::CORONARY){// Jongmin & Hyunjin
        double* time;
        double* p_lv;
        int num;
        seg->getBoundCoronaryValues(&p_lv, &time, &num);
        subdomain->SetBoundCoronaryValues(time, p_lv,num);
        cout<<"CORONARY boundary condition"<<endl;

      }else{
        subdomain -> SetBoundValue(boundV);
      }


      temp += nels + 1;
      subdomainList.push_back(subdomain);
      if(seg->getIsOutletInfo()){
        outletList.push_back(i);
      }
      if(seg->getIsOutletInfo() == false){
        subdomain->SetBoundCondition(BoundCondTypeScope::NOBOUND);
      }
    }

    // For branch use.
    temp *= 2;
    for(i=0; i<ij; i++){
      jointList[i]->SetGlobal1stLagNodeID(temp);
      temp += jointList[i]->getNumberOfSegments();
    }
}

void cvOneDBFSolver::SetDeltaTime(double dt){deltaTime = dt;}
void cvOneDBFSolver::SetStepSize(long size){stepSize = size;}
void cvOneDBFSolver::SetInletBCType(BoundCondType bc){inletBCtype = bc;}
void cvOneDBFSolver::SetMaxStep(long maxs){maxStep = maxs;}
void cvOneDBFSolver::SetQuadPoints(long point){quadPoints = point;}
void cvOneDBFSolver::SetConvergenceCriteria(double conv){convCriteria = conv;}

void cvOneDBFSolver::CreateGlobalArrays(void){
    assert( wasSet == false);
    long neq = mathModels[0]->GetTotalNumberOfEquations();

    long* maxa = new long[neq + 1];
    assert( maxa != 0);
    clear( neq + 1, maxa);

    int neeq = mathModels[0]->GetNumberOfElementEquations();
    cout <<"Number of equations " << neq << endl;
    long* eqNumbers = new long[neeq];
    assert( eqNumbers != 0);
    long minEq, total;
    int i;

    // Element nodes
    total = 0;
    for(i = 0; i < subdomainList.size(); i++){
      total += subdomainList[i]->GetNumberOfNodes();
      for( long el = 0; el < subdomainList[i]->GetNumberOfElements(); el++){
        mathModels[0]->GetEquationNumbers(el, eqNumbers, i);
        minEq = min(neeq, eqNumbers);
        for( int k = 0; k < neeq; k++){
          long currHeight = maxa[eqNumbers[k]];
          maxa[eqNumbers[k]] = max(currHeight, eqNumbers[k] - minEq);
        }
      }
    }

    // Joints
    total *= 2;
    for(i = 0; i < jointList.size(); i++){
      for(int j = 0; j < jointList[i]->getNumberOfSegments(); j++){
        minEq = mathModels[1]->GetUpmostEqnNumber(j, i);
        maxa[total] = total - minEq;
        total ++;
      }
    }

    // Now maxa contains the column heights
    // change it to hold the position
    // of the first element of the skyline at each column
    maxa[neq] = sum(neq, maxa);
    for( i = neq - 1; i >= 0; i--)
        maxa[i] = maxa[i+1] - maxa[i];

    // INITIALIZE MATRIX STORAGE SCHEME
    // AND ASSOCIATED SOLVER
# ifdef USE_SKYLINE
    lhs = new cvOneDSkylineMatrix(neq, maxa, "globalMatrix");
    cvOneDGlobal::solver = new cvOneDSkylineLinearSolver();
# endif

# ifdef USE_SUPERLU
    lhs = new cvOneDSparseMatrix(neq, maxa, "globalMatrix");
    cvOneDGlobal::solver = new cvOneDSparseLinearSolver();
# endif

# ifdef USE_CSPARSE
    lhs = new cvOneDSparseMatrix(neq, maxa, "globalMatrix");
    cvOneDGlobal::solver = new cvOneDSparseLinearSolver();
# endif

    assert(lhs != 0);

    rhs = new cvOneDFEAVector( neq);
    assert(rhs != 0);

    cvOneDGlobal::solver->SetLHS(lhs);
    cvOneDGlobal::solver->SetRHS(rhs);

    previousSolution = new cvOneDFEAVector(neq, "previousSolution");
    assert(previousSolution != 0);
    previousSolution->Clear();
    currentSolution = new cvOneDFEAVector(neq, "currentSolution");
    assert(currentSolution != 0);
    currentSolution->Clear();
    increment = new cvOneDFEAVector(neq, "increment");
    assert(increment != 0);
    increment->Clear();
}

// Initialize the solution, that is, area as area input and flow rate as 0 except the inlet
void cvOneDBFSolver::CalcInitProps(long ID){
  double segLen = subdomainList[ID] -> GetLength();
  double Qo, dQ0dT;
  Qo = subdomainList[ID] -> GetInitialFlow();
  dQ0dT=0;

  double So = subdomainList[ID] -> GetInitInletS();
  double Sn = subdomainList[ID] -> GetInitOutletS();
  for( long node = 0; node < subdomainList[ID]->GetNumberOfNodes(); node++){
  double zn = subdomainList[ID]->GetNodalCoordinate( node);
  long eqNumbers[2];
  mathModels[0]->GetNodalEquationNumbers(node, eqNumbers, ID);

  // Linear Interpolation
  double zi = (zn - segLen)/(0.0-segLen);
  double Si = (zi*(So - Sn)) + Sn;
  (*previousSolution)[eqNumbers[0]] = Si;

    if(node == 0){
      (*previousSolution)[eqNumbers[1]] = Qo;
    }else{
      (*previousSolution)[eqNumbers[1]] = 0.0;
    }
  }
}

// =================
// GENERATE SOLUTION
// =================
void cvOneDBFSolver::GenerateSolution(void){
  currentTime = 0.0;
  long i = 0;
  clock_t tstart_iter;
  clock_t tstart_solve;
  clock_t tend_iter;
  clock_t tend_solve;


  // Print the formulation used

  if(cvOneDGlobal::CONSERVATION_FORM){
    cout << "Using Conservative Form ..." << endl;
  }else{
    cout << "Using Advective Form ..." << endl;
  }

  // Allocate the TotalSolution Array.
  cout << "maxStep/stepSize: " << maxStep/stepSize << endl;
  long numSteps = maxStep/stepSize;
  TotalSolution.SetSize(numSteps+1, currentSolution -> GetDimension());
  cout << "Total Solution is: " << numSteps << " x ";
  cout << currentSolution -> GetDimension() << endl;

  cvOneDString String1( "step_");
  char String2[] = "99999";
  cvOneDString title;

  previousSolution->Rename( "step_0");
  *currentSolution = *previousSolution;

  double* tmp = previousSolution -> GetEntries();
  int j;
  for (j=0;j<previousSolution -> GetDimension(); j++){
    TotalSolution[0][j] = tmp[j];
  }

  // Initialize the Equations...
  int numMath = mathModels.size();
  for(i = 0; i < numMath; i++){
    mathModels[i]->EquationInitialize(previousSolution, currentSolution);
  }

  double cycleTime = mathModels[0]->GetCycleTime();

  // Global Solution Loop
  long q=1;
  double checkMass = 0;
  int numberOfCycle = 1;
  long iter_total = 0;

  // Time stepping
  for(long step = 1; step <= maxStep; step++){
    increment->Clear();
    for(i = 0; i < numMath; i++){
      mathModels[i]->TimeUpdate(currentTime, deltaTime);
    }

    // Newton-Raphson Iterations...
    int iter = 0;
    double normf = 1.0;
    double norms = 1.0;

    if(fmod(currentTime, cycleTime) <5.0E-6 || -(fmod(currentTime,cycleTime)-cycleTime)<5.0E-6) {
      checkMass = 0;
      cout << "**** Time cycle " << numberOfCycle++ << endl;
    }
    currentTime += deltaTime;

    while(true){
      tstart_iter=clock();

      for(i = 0; i < numMath; i++){
        mathModels[i]->FormNewton(lhs, rhs);
      }

      // PRINT RHS BEFORE BC APP
      if(cvOneDGlobal::debugMode){
        printf("(Debug) Printing LHS and RHS...\n");
        ofstream ofsRHS;
        ofstream ofsLHS;
        ofsRHS.open("rhs_1.txt");
        ofsLHS.open("lhs_1.txt");
        ofsRHS<<" --- RHS: Before ApplyBoundaryConditions" << endl;
        rhs->print(ofsRHS);
        ofsLHS<<" --- LHS: Before ApplyBoundaryConditions" << endl;
        lhs->print(ofsLHS);
        ofsRHS.close();
        ofsLHS.close();
        printf("ECCOLO\n");
        getchar();
      }

      mathModels[0]->ApplyBoundaryConditions();

      // PRINT RHS AFTER BC APP
      if(cvOneDGlobal::debugMode){
        cout<<" --- RHS: After Application of BC " << endl;
        rhs->print(cout);
      }

      // Do not evaluate residuals of lagrange eqns
      if(jointList.size() != 0){
        normf = rhs->Norm(L2_norm,1,2, jointList[0]->GetGlobal1stLagNodeID());
        norms = rhs->Norm(L2_norm,0,2, jointList[0]->GetGlobal1stLagNodeID());
      }else{
        normf = rhs->Norm(L2_norm,1,2);
        norms = rhs->Norm(L2_norm,0,2);
      }

      if (std::isnan(norms) || std::isnan(normf)) {
          throw cvException("Calculated a NaN for the residual.");
      }

      // Check Newton-Raphson Convergence
      if((currentTime != deltaTime || (currentTime == deltaTime && iter != 0)) && normf < convCriteria && norms < convCriteria){
        cout << "    iter: " << std::to_string(iter) << " ";
        cout << "normf: " << normf << " ";
        cout << "norms: " << norms << " ";
        cout << "time: " << ((float)(tend_iter-tstart_iter))/CLOCKS_PER_SEC << endl;
        break;
      }

      // Add increment
      increment->Clear();

      cvOneDGlobal::solver->Solve(*increment);

      currentSolution->Add(*increment);

      // If the area goes less than zero, it tells in which segment the error occurs.
      // Assumes that all the lagrange multipliers are at the end of the vector.
      int negArea=0;
      if(jointList.size() != 0){
        for (long i= 0; i< jointList[0]->GetGlobal1stLagNodeID();i+=2){
          long elCount = 0;
          int fileIter = 0;
          //check if area <0 or =nan
          if (currentSolution->Get(i) < 0.0 || (currentSolution->Get(i) != currentSolution->Get(i))){
           negArea=1;
            while (fileIter < model -> getNumberOfSegments()){
              cvOneDSegment *curSeg = model -> getSegment(fileIter);
              long numEls = curSeg -> getNumElements();
              long startOut = elCount;
              long finishOut = elCount + ((numEls+1)*2);
              char *modelname;
              char *segname;
              if (startOut <= i && i <= finishOut) {
                 modelname = model-> getModelName();
                 segname = curSeg -> getSegmentName();
                 std::string msg = "ERROR: The area of segment '" + std::string(segname) + "' is negative.";
                 throw cvException(msg.c_str());
              }
              elCount += 2*(numEls+1);
              fileIter++;
             }
           }
         }
        }

        if(negArea==1) {
        postprocess_Text();
        assert(0);
        }

      if(cvOneDGlobal::debugMode){
        printf("(Debug) Printing Solution...\n");
        ofstream ofs("solution.txt");
        for(int loopA=0;loopA<currentSolution->GetDimension();loopA++){
          ofs << to_string(loopA) << " " << currentSolution->Get(i) << endl;
        }
        ofs.close();
        getchar();
      }


      // A flag in case the cross sectional area is negative, but don't want to include the lagrange multipliers
      // Assumes that all the lagrange multipliers are at the end of the vector
      if(jointList.size() != 0){
        currentSolution->CheckPositive(0,2,jointList[0]->GetGlobal1stLagNodeID());
      }else{
        currentSolution->CheckPositive(0,2,currentSolution->GetDimension());
      }

      // Set Boundary Conditions
      mathModels[0]->SetBoundaryConditions();
      tend_iter=clock();

      cout << "    iter: " << std::to_string(iter) << " ";
      cout << "normf: " << normf << " ";
      cout << "norms: " << norms << " ";
      cout << "time: " << ((float)(tend_iter-tstart_iter))/CLOCKS_PER_SEC << endl;


      if(iter > MAX_NONLINEAR_ITERATIONS){
        cout << "Error: Newton not converged, exceed max iterations" << endl;
        cout << "norm of Flow rate:" << normf << ", norm of Area:" << norms << endl;
        break;
      }

    // Increment Iteration Number
    iter++;

  }// End while

  checkMass += mathModels[0]->CheckMassBalance() * deltaTime;
  cout << "  Time = " << currentTime << ", ";
  cout << "Mass = " << checkMass << ", ";
  cout << "Tot iters = " << std::to_string(iter) << endl;

  // Save solution if needed
  if(step % stepSize == 0){
    sprintf( String2, "%ld", (unsigned long)step);
    title = String1 + String2;
    currentSolution->Rename(title.data());

    double * tmp = currentSolution -> GetEntries();
    int j;

    for(j=0;j<currentSolution -> GetDimension(); j++){
      TotalSolution[q][j] = tmp[j];
    }
    q++;
  }
  *previousSolution = *currentSolution;
  iter_total += iter;
  } // End global loop

  cout << "\nAvgerage number of Newton-Raphson iterations per time step = "<<(double)iter_total / (double)maxStep<<"\n"<< endl;
}
