#ifndef CVONEDDATATABLE_H
#define CVONEDDATATABLE_H

# include <string>
# include "cvOneDTypes.h"

using namespace std;

class cvOneDDataTable{
  protected:
    string name;
    string type;
    cvDoubleVec time;
    cvDoubleVec values;

  public:
    cvOneDDataTable();
    virtual ~cvOneDDataTable();

    // Getters and Setters
    string getName(){return this->name;}
    string getType(){return this->type;}
    double getTime(int index){return this->time[index];}
    double getValues(int index){return this->values[index];}
    int getSize(){return this->values.size();}

    void setName(string name);
    void setType(string type);
    void setTime(cvDoubleVec time);
    void setValues(cvDoubleVec values);

};

#endif // CVONEDDATATABLE_H
