#ifndef CVONEDMODEL_H
#define CVONEDMODEL_H

//
//  cvOneDModel.h: Class to handle 1D network Models.
//  ~~~~~~~
//
//  This is the C representation of the model
//  It corresponds roughly to the java version, 
//  but not exactly.
//
//  History:
//  Jul. 2003, J.Wan
//      *Use STL vector to replace old darray for better memory management and optimization
//      *Deleted bloodpropstruct.h. The struct declaration is moved here instead.
//  May 1999, B.Steele, S.A.Spicer and J.Wan
//      Creation of file

# include <vector>

# include "cvOneDEnums.h"
# include "cvOneDSegment.h"
# include "cvOneDNode.h"
# include "cvOneDJoint.h"
# include "cvOneDMatrix.h"
 
using namespace std;

struct cvOneDBloodPropStruct{
  double viscosity;
  double density;
};

class cvOneDModel{

  public:    
    
    // Default Constructor/Destructor
    cvOneDModel();
    ~cvOneDModel(); 

    // Safe Constructor
    static cvOneDModel * New(void); 

    // Safe Destructor
    void Delete(void);
    

    // Accessors     
    void      setModelName(char *);
    char *    getModelName(void);
    
    void      setModelID(long);
    long      getModelID(void);
    
    long      getNumberOfSegments(void);
    long      getNumberOfNodes(void);
    long      getNumberOfJoints(void);
    cvOneDSegment*  getSegment(long id);
    cvOneDNode*    getNode(long id);
    cvOneDJoint*    getJoint(long id);
    long  getTopJoint(void);
    void  setTopJoint(long t);

    // Finite Element Quantities
    long getNumberOfEquations(void);
    
    // Actions 
    int addSegment(cvOneDSegment *newSeg);
    int addNode(cvOneDNode *Node);
    int addJoint(cvOneDJoint *Joint);

    // just to see what's happening
    void PrintModel(void);

    int Reorder(void);
    
    
    cvOneDBloodPropStruct   BloodProps;
    
  private:

    long FindJointinSegs(vector<long>& results, 
                   vector<long>& nextjoint, 
                   long startseg, 
                   long jointnr);

    void SwapSegs(long seg1, long seg2);

    char modelName[2048];
    long modelID;

    long numEquations;
    long topJoint;
    // double    convCrit;

    vector<cvOneDNode*>     NodeList;
    vector<cvOneDJoint*>    JointList;
    vector<cvOneDSegment*>  SegmentList;

    // Mapping back to Brooke's Java
    vector<long>  Mapping;

    // How many models are there?
    static long NumModels;

};

#endif // CVONEDMODEL_H
