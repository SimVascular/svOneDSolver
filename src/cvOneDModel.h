#ifndef CVONEDMODEL_H
#define CVONEDMODEL_H

//
//  cvOneDModel.h: Class to handle 1D network Models.
//

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

    vector<cvOneDNode*>     NodeList;
    vector<cvOneDJoint*>    JointList;
    vector<cvOneDSegment*>  SegmentList;

    // Mapping of Segments
    vector<long>  Mapping;

    // How many models are there?
    static long NumModels;

};

#endif // CVONEDMODEL_H
