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

#include <string>
#include <math.h>
#include <memory>
#include <optional>
#include "cvOneDTypes.h"
#include "cvOneDException.h"

using namespace std;

namespace cvOneD{

struct options{

    // DATA
    string modelName;

    // NODE DATA
    cvStringVec nodeName;
    cvDoubleVec nodeXcoord;
    cvDoubleVec nodeYcoord;
    cvDoubleVec nodeZcoord;

    // JOINT DATA
    cvStringVec jointName;
    // jointNode: name of the associated node that we will use 
    // to get the XYZ coordinates when creating the joint.
    cvStringVec jointNode; 

    cvStringVec jointInletName;
    cvStringVec jointOutletName;

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

    // COUPLING OPTIONS
    string couplingStatus;      // "ON" or "OFF"
    string couplingType;        // "DIR" or "NEU"
    int couplingSubsteps;       // integer value (e.g., 100)

    // These are to preserve legacy behavior and are 
    // expected to be eventually migrated into a 
    // post-processing step.
    string outputType = "TEXT";
    std::optional<int> vtkOutputType = std::nullopt; 
};

void validateOptions(options const& opts);

} // namespace cvOneD

#endif // CVONEDOPTIONS_H

