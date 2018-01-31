#ifndef CVONEDJOINT_H
#define CVONEDJOINT_H

//
//  cvOneDJoint.h: A data structure for model joints. 
//

#include <vector>

#include "cvOneDEnums.h"
#include "cvOneDSegment.h"

using namespace std;

struct cvOneDBoundCondsJoint{
  cvOneDBoundCondsJoint(){}
  ~cvOneDBoundCondsJoint(){}
  BoundCondType bcType;
  cvOneDBoundCondsJoint& operator=(const cvOneDBoundCondsJoint& in){return *this;}
  double bcParameter1; 
  double bcParameter2;
};

struct cvOneDJoint{
  cvOneDJoint(){}
  ~cvOneDJoint(){}
  cvOneDJoint& operator=(const cvOneDJoint& in){return *this;}    
  void setJointID(int i){id = i;}
  int getNumberOfSegments(){return InletSegments.size()+OutletSegments.size();}
  char Name[2048];
  int id;
  double x, y, z, radius;
  vector<cvOneDBoundCondsJoint*> theBCs;
  vector<int> InletSegments;
  vector<int> OutletSegments;
};

#endif // CVONEDJOINT_H

