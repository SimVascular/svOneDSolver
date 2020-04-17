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

#ifndef CVONEDOPTIONS_H
#define CVONEDOPTIONS_H

# include <string>
# include <math.h>
# include "cvOneDTypes.h"
# include "cvOneDException.h"

using namespace std;

class cvOneDOptions{
  public:
    // CONSTRUCTOR
    cvOneDOptions();
    // DESTRUCTOR
    ~cvOneDOptions();

    // DATA
    string modelName;
    bool modelNameDefined;

    // NODE DATA
    cvStringVec nodeName;
    cvDoubleVec nodeXcoord;
    cvDoubleVec nodeYcoord;
    cvDoubleVec nodeZcoord;

    // JOINT DATA
    cvStringVec jointName;
    cvStringVec jointNode;
    cvDoubleVec jointXcoord;
    cvDoubleVec jointYcoord;
    cvDoubleVec jointZcoord;
    cvStringVec jointInletName;
    cvStringVec jointOutletName;

    // JOINT INLET AND OUTLET LIST
    cvStringVec jointInletListNames;
    cvLongVec   jointInletListNumber;
    cvLongMat   jointInletList;
    cvStringVec jointOutletListNames;
    cvLongVec   jointOutletListNumber;
    cvLongMat   jointOutletList;

    // MATERIAL
    cvStringVec materialName;
    cvStringVec materialType;
    cvDoubleVec materialDensity;
    cvDoubleVec materialViscosity;
    cvDoubleVec materialPRef;
    cvDoubleVec materialExponent;
    cvDoubleVec materialParam1;
    cvDoubleVec materialParam2;
    cvDoubleVec materialParam3;

    // DATATABLE
    cvStringVec dataTableName;
    cvStringVec dataTableType;
    cvDoubleMat dataTableVals;

    // SEGMENT DATA
    cvStringVec segmentName;
    cvLongVec   segmentID;
    cvDoubleVec segmentLength;
    cvLongVec   segmentTotEls;
    cvLongVec   segmentInNode;
    cvLongVec   segmentOutNode;
    cvDoubleVec segmentInInletArea;
    cvDoubleVec segmentInOutletArea;
    cvDoubleVec segmentInFlow;
    cvStringVec segmentMatName;
    cvStringVec segmentLossType;
    cvDoubleVec segmentBranchAngle;
    cvLongVec   segmentUpstreamSegment;
    cvLongVec   segmentBranchSegment;
    cvStringVec segmentBoundType;
    cvStringVec segmentDataTableName;

    // SOLVER OPTIONS
    double timeStep;
    long   stepSize;
    long   maxStep;
    long   quadPoints;
    string inletDataTableName;
    string boundaryType;
    double convergenceTolerance;
    int    useIV;
    int    useStab;
    string outputType;
    bool solverOptionDefined;

    // PRINTING
    void printToFile(string fileName);
    void printModelName(FILE* f);
    void printNodeData(FILE* f);
    void printJointData(FILE* f);
    void printSegmentData(FILE* f);
    void printSolverOptions(FILE* f);
    void printMaterialData(FILE* f);
    void printJointInletData(FILE* f);
    void printJointOutletData(FILE* f);
    void printDataTables(FILE* f);

    // CHECKING
    void check();
    void checkSegmentLengthConsistency();
    void checkSegmentHasNodes();
    void checkJointHasNodes();
};

#endif // CVONEDOPTIONS_H

