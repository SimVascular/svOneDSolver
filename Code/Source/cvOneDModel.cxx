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



