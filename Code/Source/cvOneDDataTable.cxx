#include "cvOneDDataTable.h"

// Constructor
cvOneDDataTable::cvOneDDataTable(){
}
// Destructor
cvOneDDataTable::~cvOneDDataTable(){
}

void cvOneDDataTable::setName(string name){
  this->name = name;
}

void cvOneDDataTable::setType(string type){
  this->type = type;
}

void cvOneDDataTable::setTime(cvDoubleVec time){
  this->time.clear();
  for(int loopA=0;loopA<time.size();loopA++){
    this->time.push_back(time[loopA]);
  }
}

void cvOneDDataTable::setValues(cvDoubleVec values){
  this->values.clear();
  for(int loopA=0;loopA<values.size();loopA++){
    this->values.push_back(values[loopA]);
  }
}

