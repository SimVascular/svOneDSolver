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
