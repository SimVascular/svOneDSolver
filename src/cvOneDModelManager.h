#ifndef CVONEDMODELMANAGER_H
#define CVONEDMODELMANAGER_H

# include <string.h>
# include <vector>

# include <boost/algorithm/string.hpp>

# include "cvOneDTypes.h"
# include "cvOneDGlobal.h"
# include "cvOneDException.h"
# include "cvOneDModel.h"
# include "cvOneDMetteImpedance.h"
# include "cvOneDMorphImpedance.h"

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
                   double* times,double conv, int useIV, int usestab, int useshockcap, double smoothbeta, double Sref, double Qref);

    // COMPUTE ADMITTANCE AND IMPEDANCE CURVES
    int CalcImpedance(int numTimeSteps, double lengthRadiusRatio, double rootRadius, double period, double scaleFactor, int rtnFourier,
                      cvDoubleVec& time, cvDoubleVec& values);
    int CalcImpedanceRCR(int numTimeSteps, double Rd, double Ru, double period, double Cap, double scaleFactor, int rtnFourier,
                         cvDoubleVec& time, cvDoubleVec& values);
    int CalcMorphImpedance(int numTimeSteps, int min_order, double rootRad, double period, int rtnFourier,
                           cvDoubleVec& time, cvDoubleVec& values);
    int CalcAdmittance(int numTimeSteps,double ratio,double rootRad,double period,double scaleFactor,int rtnFourier,
                       cvDoubleVec& time, cvDoubleVec& values);
    int CalcRCRAdmittance(int numTimeSteps,double Rd,double Ru,double period,double Cap,double scaleFactor,int rtnFourier,
                          cvDoubleVec& time, cvDoubleVec& values);
    int CalcMorphAdmittance(int numTimeSteps,int min_order,double rootRad,double period,int rtnFourier,
                            cvDoubleVec& time, cvDoubleVec& values);
    int CalcLRR(int numTimeSteps,double resistance,double rootRad,double period,double scaleFactor,
                double& lrr);

};

#endif // CVONEDMODELMANAGER_H
