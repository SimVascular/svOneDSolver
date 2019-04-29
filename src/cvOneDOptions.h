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
    int    useShockcap;
    double smoothbeta;
    double Sref;
    double Qref;
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
