#ifndef CVONEDFEAJOINT_H
#define CVONEDFEAJOINT_H

//  cvOneDFEAJoint.h: A data structure for Lagrange joints. 
//  ~~~~~~~~~~
//
//  History:
//  Jul. 2003, J.Wan
//      Use stl vector for better memory management and optimization
//  May 1999, J.Wan, B.Steele, S.A.Spicer and S.Strohband
//      Creation of file, class project of ME234C of T.J.R. Hughes and C.Taylor

# include <vector>

# include "cvOneDSubdomain.h"

using namespace std;

class cvOneDFEAJoint{
  public:
    cvOneDFEAJoint() {;}
    ~cvOneDFEAJoint(){;}
    void setJointID(int i){id = i;}
    // total number of segments
    int getNumberOfSegments(){return InletSubdomains.size()+OutletSubdomains.size();}
    // inlet segments
    int getNumberOfInletSegments(){return InletSubdomains.size();}
    // outlet segments
    int getNumberOfOutletSegments(){return OutletSubdomains.size();}
    // add an inlet segment by its global segment index
    void AddInletSubdomains(int i){InletSubdomains.push_back(i);}
    // add an inlet segment by its global segment index
    void AddOutletSubdomains(int i){OutletSubdomains.push_back(i);}
    // get the global segment index of ith inlet segment
    int GetInletID(int ith){return InletSubdomains[ith];}
    // get the global segment index of ith outlet segment
    int GetOutletID(int ith){return OutletSubdomains[ith];}
    long GetGlobal1stLagNodeID(void){return global1stLagNodeID;}
    void SetGlobal1stLagNodeID(long ID){global1stLagNodeID = ID;}
  private:
    long global1stLagNodeID;
    int id;
    vector<int> InletSubdomains;
    vector<int> OutletSubdomains;
};

#endif // CVONEDFEAJOINT_H
