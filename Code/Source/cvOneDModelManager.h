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

#ifndef CVONEDMODELMANAGER_H
#define CVONEDMODELMANAGER_H

# include <string.h>
# include <vector>

# include "cvOneDTypes.h"
# include "cvOneDGlobal.h"
# include "cvOneDException.h"
# include "cvOneDModel.h"

class cvOneDModelManager{
  public:
    // CONSTRUCTOR
    cvOneDModelManager(char *mdlName);
    // DESTUCTOR
    ~cvOneDModelManager();

    // CREATE MATERIAL
    int CreateMaterial(char *matName, char *MaterialTypeString,
                       double density, double dynamicViscosity,
                       double profile_exponent, double pRef,
                       int numParams, double *params, int *matID);

    // CREATE SEGMENT
    int CreateSegment(char* segName,long segID, double  segLen,
                      long numEls,long inNode,long outNode,
                      double InitialInletArea,double InitialOutletArea,
                      double InitialFlow,int matID,char* lossType,
                      double branchAngle,int upstreamSegment,int branchSegment,
                      char* boundType,double* value, double* time, int num );

    // CREATE DATATABLE
    int CreateDataTable(char* dtName,char* dtType, cvDoubleVec values);

    // CREATE NODE
    int CreateNode(char* nodeName,double x,double y,double z);


    // CREATE JOINT
    int CreateJoint(char* jointName,double x,double y,double z,
                    int numInSegs,int numOutSegs,
                    int *InSegs,int *OutSegs);

    // SOLVE MODEL
    int SolveModel(double dt,long stepSize,
                   long maxStep,long quadPoints,
                   int len,char* boundType,double* values,
                   double* times,double conv, int useIV, int usestab);

};

#endif // CVONEDMODELMANAGER_H
