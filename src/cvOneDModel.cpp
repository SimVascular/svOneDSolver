//
//  Model.cxx: Class to handle 1D network Models
//

# include <string.h>

# include "cvOneDModel.h"

using namespace std;

//Static Declarations
long cvOneDModel::NumModels;

cvOneDModel::cvOneDModel(){
  cvOneDModel::NumModels++;
}

cvOneDModel::~cvOneDModel(){
  cvOneDModel::NumModels--;
}

cvOneDModel * cvOneDModel::New(void){
  cvOneDModel::NumModels++;
  return new cvOneDModel;
}

void cvOneDModel::Delete(void){
  cvOneDModel::NumModels--;
  delete this;    
}
 
void cvOneDModel::setModelID(long id){
  modelID = id;
}

long cvOneDModel::getModelID(void){
  return modelID;
}

void cvOneDModel::setModelName(char *name){
  modelName[0] = '\0';
  strcpy(modelName,name);
}

char *cvOneDModel::getModelName(void){
  return modelName;
}

long cvOneDModel::getNumberOfSegments(void){
  return (long)SegmentList.size();
}

long cvOneDModel::getNumberOfNodes(void){
  return (long)NodeList.size();
}

long cvOneDModel::getNumberOfJoints(void){
  return (long)JointList.size();
}

cvOneDSegment * cvOneDModel::getSegment(long id){
  return SegmentList[Mapping[id]];
}

cvOneDNode* cvOneDModel::getNode(long id){
  return NodeList[id];
}

cvOneDJoint* cvOneDModel::getJoint(long id){
  return JointList[id];
}

long cvOneDModel::getNumberOfEquations(void){
  return numEquations;
}

int cvOneDModel::addNode(cvOneDNode * newNode){
  NodeList.push_back(newNode);
  return 1;
}

int cvOneDModel::addJoint(cvOneDJoint * newJoint){
  JointList.push_back(newJoint);
  return 1;
}

int cvOneDModel::addSegment(cvOneDSegment *newSeg){

  SegmentList.resize(SegmentList.size() + 1);

  if(Mapping.size() < newSeg -> getSegmentID()+1){
    Mapping.resize(newSeg -> getSegmentID()+1);
  }
  SegmentList[SegmentList.size()-1] = newSeg;
  Mapping[newSeg -> getSegmentID()] = SegmentList.size() - 1; 
  return 1;
}

long cvOneDModel::getTopJoint(void){
  return (topJoint);
}

void cvOneDModel::setTopJoint(long t){
  topJoint=t;
}

void cvOneDModel::SwapSegs(long seg1, long seg2){
  cvOneDSegment* tmp;
  long tmpMap;
  tmp=SegmentList[seg1];
  tmpMap=Mapping[seg1];
  SegmentList[seg1]=SegmentList[seg2];
  Mapping[seg1]=Mapping[seg2];
  SegmentList[seg2]=tmp;
  Mapping[seg2]=tmpMap;
}

long cvOneDModel::FindJointinSegs(vector<long>& segs, 
                                  vector<long>& nextjoint, 
                                  long startseg, 
                                  long jointnr){
  long i, *c, n;
  n=0;
  for (i=startseg;i<getNumberOfSegments();i++){
    c=SegmentList[i]->getInOutJoints();
    if (c[0] == jointnr){ 
      n++; 
      segs.resize(n); 
      nextjoint.resize(n);
      segs[n-1]=i;
      nextjoint[n-1]=c[1];
    }
    if (c[1] == jointnr){ 
      n++; 
      segs.resize(n); 
      nextjoint.resize(n);
      segs[n-1]=i;
      nextjoint[n-1]=c[0];
    }
  }
  return n;
}



