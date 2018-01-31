#ifndef CVONEDNODE_H
#define CVONEDNODE_H

//
//  cvOneDNode.h: A data structure for model nodes. 
//

#include <vector>
#include "cvOneDEnums.h"
#include "cvOneDSegment.h"

using namespace std;
struct cvOneDBoundConds{

  cvOneDBoundConds(){}
  ~cvOneDBoundConds(){}
  BoundCondType bcType;
  cvOneDBoundConds& operator=(const cvOneDBoundConds& in){return *this;}
  double bcParameter1; double bcParameter2;
};

struct cvOneDNode{
  cvOneDNode(){}
  ~cvOneDNode(){}
  cvOneDNode& operator=(const cvOneDNode& in){return *this;}    
  void setNodeID(int i){id = i;}
  char  Name[2048];
  int id;
  double x, y, z, radius;
  vector<cvOneDBoundConds*> theBCs;
};

#endif // CVONEDNODE_H

