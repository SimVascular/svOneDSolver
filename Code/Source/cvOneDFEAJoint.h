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

#ifndef CVONEDFEAJOINT_H
#define CVONEDFEAJOINT_H

//  cvOneDFEAJoint.h: A data structure for Lagrange joints. 
//  ~~~~~~~~~~
//

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
