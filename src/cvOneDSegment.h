#ifndef CVONEDSEGMENT_H
#define CVONEDSEGMENT_H

//
//  cvOneDSegment.h: Class to handle 1D network segments
//  ~~~~~~~~~~~
//
//  This the C representation of a model segment.  
//

#include "cvOneDEnums.h"
#include "cvOneDMesh.h"

typedef void (*MeshEvaluatorFunctionPtr)();

class cvOneDSegment{

 public:    
   
    // Default Constructor/Destructor
    cvOneDSegment();
    cvOneDSegment(double IA, double FA, double IF, bool IO);
    ~cvOneDSegment();

    // Safe Constructor
    static cvOneDSegment * New(double IA, double FA, double IF, bool IO);
    
    // Safe Destructor
    void Delete(void);

    // SegmentID/Name/Parent
    void     setSegmentID(long id);
    long     getSegmentID(void);
    void     setSegmentName(char *Name);
    char    *getSegmentName(void);
    void     setParentModel(void *mdl);
    void    *getParentModel(void);

    // Segment Length/Numels/Type (defines mesh spacing);
    void     setSegmentLength(double len);
    double   getSegmentLength(void);
    void     setZOfTwoEnds(double zin, double zout);
    double   getInletZ();
    double   getOutletZ();
    void     setNumElements(long nels);
    long     getNumElements(void);
    double   getInitInletS(void);
    double   getInitOutletS(void);
    double   getInitialFlow(void);
    double   getInitialPressure(void);
    double   getInitialdFlowdT(void);
    void     setIsOutlet(void){IsOutlet = true;} 
    bool     getIsOutletInfo(void){return IsOutlet;}
    
    // Minor Loss Info
    //bool     getIsStenosisInfo(void){return IsStenosis;} 
    MinorLoss GetMinorLossType(void){return lossType;}
    void  SetMinorLossType(MinorLoss loss){ lossType = loss;}
    double  GetBranchAngle(void){return branchAngle;}
    void  SetBranchAngle(double angle){branchAngle=angle;}
    int    GetUpstreamSeg(void){return upstream;}
    void  SetUpstreamSeg(int up){upstream = up;}
    int    GetBranchSeg(void){return branch;}
    void  SetBranchSeg(int up){branch = up;}

    // Mesh Topolgy on Segment (Use Evaluator for user defined mesh)
    void     setMeshType(MeshType mtype);
    MeshType getMeshType(void);

    MeshEvaluatorFunctionPtr setMeshEvaluator(MeshEvaluatorFunctionPtr fcn);
    
    // Material Type

    void setMaterialID(int matID) {matID_ = matID;}
    int  getMaterialID() {return matID_;}
        
    int use_k;

    // Boundary conditions if existent
    void  setBoundCondition(BoundCondType bound){boundType = bound;}
    BoundCondType getBoundCondition(void){return boundType;}

    // SET GENERIC BOUNDARY VALUE
    void  setBoundValue(double value){boundValue = value;}

    // SET BOUNDARY PRESSURE VALUE
    void  setBoundPressureValue(double* pval, double* ptime,int num){
      presLength = num;
      if(values != NULL){
        delete [] values;
        values = NULL;
      }
      if(times != NULL){
        delete [] times;
        times = NULL;
      }
      values = new double[num];
      times = new double[num];
      for(int loopA=0;loopA<num;loopA++){
        values[loopA] = pval[loopA];
        times[loopA] = ptime[loopA];
      }
    }

    // SET BOUNDARY RESISTANCE VALUE
    void  setBoundResistanceValue(double* pval, double* ptime,int num){
      presLength = num;
      if(values != NULL){
        delete [] values;
        values = NULL;
      }
      if(times != NULL){
        delete [] times;
        times = NULL;
      }
      values = new double[num];
      times = new double[num];
      for(int loopA=0;loopA<num;loopA++){
        values[loopA] = pval[loopA];
        times[loopA] = ptime[loopA];
      }
    }

    // SET RCR BOUNDARY CONDITION
    void setBoundRCRValue(double* rcr, int num){//added IV 050803
      rpCapRdLength = num;
      if(rpCapRd != NULL){
        delete [] rpCapRd;
        rpCapRd = NULL;
      }
      rpCapRd = new double[num];
      for(int loopA=0;loopA<num;loopA++){
        rpCapRd[loopA] = rcr[loopA];
      }
    }  

    
    double getBoundValue(){return boundValue;}
    void getBoundPressureValues(double** value, double **time,int* num ){
      *value=values,
      *time=times,
      *num=presLength;
    }
    void getBoundResistanceValues(double** value, double **time,int* num ){
      *value=values,
      *time=times,
      *num=presLength;
    }

    void getBoundRCRValues(double** rcr, int* num){//added IV 050803
      *rcr = rpCapRd;
      *num = rpCapRdLength;
    }

    // Segment Connectivity -- These index into the 
    // global Joint list.
    void     setInOutJoints(long inJoint, long outJoint);
    long *   getInOutJoints(void);

    // Material Properties
    void     setYoungsModulus(double e);
    double   getYoungsModulus(void);
    
    // Generate the segment's mesh
    void     tesselate(void);
    void     putInMatrix(void);

    // return the mesh spacing h for a given
    // mesh element
    double   getMeshSpacing(long l);  
    bool IsOutlet;
        
 private:

    // These define the segment's Identity.
    void* parentModel;  
    long SegmentID;
    char SegmentName[2048];
    
    // These define the segment's mesh properties, the mesh 
    // evaluator to be used by tesselate(), and the mesh itself.
    double SegmentLength;
    double zin, zout;
    long NumElements;
    MeshType SegmentMeshType;
    
    MeshEvaluatorFunctionPtr meshEvaluator;

    cvOneDMesh *    segmentMesh;    

    // This defines the mesh's connectivity ConnNodes[0] is 
    // the inlet Node, ConnNodes[1] is the oulet Node.
    long ConnJoints[2];

    // These are the Material Properties of the Vessel
    //double youngs_modulus;

    // Have we tesselated the segment?
    bool isTesselated;

    // How Many Segments?
    static long NumSegments;

    // Some other info
    double InitialPressure;
    double InitInletS;
    double InitOutletS;
    double InitialFlow;
    double InitialdFlowdT;
    bool IsStenosis;
//  The boundary condition refers to the outlet, either pressure 
//  specified or resistance specified. As for the inlet, the default
//  boundary condition is constant flow rate.
    BoundCondType  boundType;
    double boundValue;  // the value is either pressure, area or resistance,
                        // which is dependent on the boundType.
    double* values;
    double* times;
    double* impedance;
    double* rpCapRd;
    double* waveValue; // *rpCapRd, *wavevalue added by IV
    int presLength;
    int impedanceLength;
    int rpCapRdLength;
    int waveValueLength;// rpCapRdLength, waveValueLength added by IV
    int dummy;

    int upstream;
    int branch;
    double branchAngle;
    MinorLoss lossType;

    int matID_;
    // MaterialType elasticType;
    // double matK[3];
};

#endif // CVONEDSEGMENT_H
    
