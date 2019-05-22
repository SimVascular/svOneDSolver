//
//  Segment.cxx: a class to handle 1D network segments
//  ~~~~~~~~~
//
//  SYNOPSIS...This the C representation of a model segment.  
//             Again, this corresponds roughly to the Java 
//             version but not exactly.
//

#include "cvOneDSegment.h"
#include "cvOneDMthModelBase.h"

// Static Declarations
long cvOneDSegment::NumSegments;

cvOneDSegment::cvOneDSegment(){

  cvOneDSegment::NumSegments++;
  InitInletS = 3.1415926;
  InitOutletS   = 3.1415926; //the ratio couldn't be too big ?
  InitialFlow = 0;
  InitialdFlowdT=0;
  InitialPressure=0;
  IsOutlet    = false;
  boundType   = BoundCondTypeScope::NOBOUND;
  lossType	= MinorLossScope::NONE;
  impedance = NULL;
  values = NULL;
  times = NULL;
  rpCapRd = NULL;
  waveValue = NULL;
  presLength = 0;
  impedanceLength = 0;
  rpCapRdLength = 0;
  waveValueLength = 0;
}

cvOneDSegment::cvOneDSegment(double IIA, double IOA, double IF, bool IO){
  cvOneDSegment::NumSegments++;
  InitInletS     = IIA; //Initial Inlet Area
  InitOutletS    = IOA; //Initial Outlet Area
  InitialFlow    = IF;
  InitialdFlowdT = 0;
  IsOutlet       = IO;
  boundType      = BoundCondTypeScope::NOBOUND;
  lossType	     = MinorLossScope::NONE;
  impedance = NULL;
  values = NULL;
  times = NULL;
  rpCapRd = NULL;
  waveValue = NULL;
  presLength = 0;
  impedanceLength = 0;
  rpCapRdLength = 0;
  waveValueLength = 0;
}

cvOneDSegment::~cvOneDSegment(){
  cvOneDSegment::NumSegments--;
}

cvOneDSegment * cvOneDSegment::New(double IA,double FA,double IF,bool IO){    
  // return new cvOneDSegment(IA, FA, IF, IO, IS);
  return new cvOneDSegment(IA, FA, IF, IO);
}

void cvOneDSegment::Delete(void){
  delete this;    
}

void cvOneDSegment::setSegmentID(long id){
  SegmentID = id;
}

long cvOneDSegment::getSegmentID(void){
  return SegmentID;
}

void cvOneDSegment::setSegmentName(char *Name){
  SegmentName[0] = '\0';
  strcpy(SegmentName,Name);
}

char *cvOneDSegment::getSegmentName(void){
  return SegmentName;
}

void cvOneDSegment::setParentModel(void *mdl){
  parentModel = mdl;
}

void *cvOneDSegment::getParentModel(void){
  return parentModel;
}

void cvOneDSegment::setSegmentLength(double len){
  SegmentLength = len;
  zin = 0;
  zout = len;
}

void cvOneDSegment::setZOfTwoEnds(double zin_, double zout_){
  zin = zin_;
  zout = zout_;
}

double cvOneDSegment::getInletZ(){
  return zin;
}

double cvOneDSegment::getOutletZ(){
  return zout;
}

double cvOneDSegment::getSegmentLength(void){
  return SegmentLength;
}

void cvOneDSegment::setNumElements(long nels){
  NumElements = nels;
}

long cvOneDSegment::getNumElements(void){
  return NumElements;
}

double cvOneDSegment::getInitInletS(void){
  return InitInletS;
}

double cvOneDSegment::getInitOutletS(void){
  return InitOutletS;
}

double cvOneDSegment::getInitialPressure(void){
  return InitialPressure;
}

double cvOneDSegment::getInitialFlow(void){
  return InitialFlow;
}

double cvOneDSegment::getInitialdFlowdT(void){
  //	return InitialdFlowdT;
	return 0;
}

void cvOneDSegment::setMeshType(MeshType mtype){
   SegmentMeshType = mtype;
}

MeshType cvOneDSegment::getMeshType(void){
  return SegmentMeshType;
}

MeshEvaluatorFunctionPtr cvOneDSegment::setMeshEvaluator(MeshEvaluatorFunctionPtr fcn){
  meshEvaluator = fcn;
  return meshEvaluator;
}

void cvOneDSegment::setInOutJoints(long inJoint, long outJoint){
  ConnJoints[0] = inJoint;
  ConnJoints[1] = outJoint;
}

long *cvOneDSegment::getInOutJoints(void){
  return ConnJoints;
}

void cvOneDSegment::tesselate(void){
  // Unimplimented...but coming
}

void cvOneDSegment::putInMatrix(void){
  // Unimplimented...but coming.
}

double cvOneDSegment::getMeshSpacing(long l){
  return segmentMesh -> ElementList[l].h;
}



