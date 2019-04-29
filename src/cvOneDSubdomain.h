#ifndef CVONEDSUBDOMAIN_H
#define CVONEDSUBDOMAIN_H

//
//  cvOneDSubdomain.h - Header for a class to contain the descritization of 
//  ~~~~~~~~~~~~~   of the Geometry.
//


#include <cstring>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <complex> //IV 080703 for wave BC

#include "cvOneDEnums.h"
#include "cvOneDFiniteElement.h"
#include "cvOneDMaterial.h"
#include "cvOneDMaterialOlufsen.h"
#include "cvOneDError.h"

class cvOneDSubdomain{

  public:

    // Constructor/Destructor
    cvOneDSubdomain();
    ~cvOneDSubdomain();

		// Setup the subdomain...
    void SetNumberOfNodes(long);
    void SetNumberOfElements(long);
    void SetMeshType(MeshType mType);

    void Init(double x0, double xL);

    long GetNumberOfNodes() const;
    long GetNumberOfElements() const;
    void GetConnectivity( long element, long* conn) const; //local
    void GetNodes( long element, double* nd) const; 
    double GetNodalCoordinate( long node) const;
	bool GetStenosisInfo() {return isStenosis;}
    
	cvOneDFiniteElement* GetElement( long element) const;

    void SetInitInletS(double So);
    void SetInitOutletS(double Sn);
    void SetFinalArea(double Sn);
    void SetInitialFlow(double Qo);
    void SetInitialPressure(double P0);
    void SetInitialdFlowdT(double dQ0dT);
    void SetStenosisInfo(bool is) {isStenosis = is;}
    double GetInletZ(){return z_in;}
    double GetOutletZ() {return z_out;}
    double GetLength() {return fabs(z_out - z_in);}
    double GetInitInletS(void);
    double GetInitOutletS(void);
    double GetFinalArea();
    double GetInitialFlow(void);
    double GetInitialPressure(void);
    double GetInitialdFlowdT(void);

	long GetGlobal1stNodeID(void) {return global1stNodeID;}
    void SetGlobal1stNodeID(const long id){global1stNodeID = id;}
    
    cvOneDMaterial* GetMaterial(void) const {return mat;}// didn't change
    void SetupMaterial(int matID);
    
    // Boundary conditions if existent
    void  SetBoundCondition(BoundCondType bound){boundType = bound;}
    BoundCondType GetBoundCondition(void){return boundType;}
    void SetBoundValue(double boundV);
    void SetBoundPresWave(double *time, double *pres, int num);
	void SetBoundResistanceWave(double *time, double *resist, int num);
	void SetBoundRCRValues(double *rcr, int num);//added IV 050803
	
	double GetBoundArea(){return boundValue;}
    double GetBoundResistance(){return boundValue;}
    double GetBoundResistance(double currentTime);
	double GetBoundFlowRate(){return boundValue;}
    double GetBoundAreabyPresWave(double currentTime);
	double GetPressure(double currentTime);
	double GetdPressuredt(double currentTime) { return 0; }

	double GetRp(){return proximalResistance;}//added IV 050803
	double GetCap(){return capacitance;}//added IV 050803
	double GetRd(){return distalResistance;}//added IV 050803
	double GetAlphaRCR(){return alphaRCR;}//added IV 050803
	
    //double* GetEigValWave(){return eigValWave;}//added IV 080703

	/*
	* Use the following functions MemIntRCR, MemAdvRCR, dMemIntRCRdP, dMemAdvRCRdP for RCR boundary conditions
	* Mem is for "memory" - time integrals
	* added IV 050803
	*/
	double MemIntRCR(double currP, double previousP, double deltaTime, double currentTime);//compute integral over one time step of the pressure convolution in time int(P(t')*exp(-alphaRCR(t-t')),0,t)
	double MemAdvRCR(double currP, double previousP, double deltaTime, double currentTime);//compute part of the advective term integral over one time step of Q^2, Q=couplingFunction(P)
	double dMemIntRCRdP(double deltaTime);//used for contribution to LHS
	double dMemAdvRCRdP(double currP, double previousP, double deltaTime, double currentTime);//contribution to LHS: compute part of the advective term integral over one time step of Q^2, Q=couplingFunction(P))
	double MemC(double currP, double previousP, double deltaTime, double currentTime);//if RCR essential-see ApplyBoundaryConditions()
	
	// minor loss / stenosis info
	
	// type of minor loss (segment, branch, etc...)
	MinorLoss GetMinorLossType(void){return minorLoss;}
	void SetMinorLossType(MinorLoss loss){minorLoss = loss;}
	
	// angle between loss segment and upstream segment
    double GetBranchAngle(void){return branchAngle;}
	void SetBranchAngle(double angle){branchAngle = angle;}

	// segment upstream of loss segment
	void SetUpstreamSeg(int seg){upstream = seg;}
	int GetUpstreamSeg(void){return upstream;}
	void SetBranchSeg(int seg){branch = seg;}
	int GetBranchSeg(void){return branch;}

	void SaveK(double K, int i);
	double* K;
	double* impedancePressure;
	double  lastImpedancePressure;

  private:
    // The initial state & dimensions.
    int ID;
    double S_initial;
    double S_final;
    double Q_initial;
    double P_initial;
    double dQ_dT_initial;
    double z_in;
    double z_out;

    long numberOfElements;
    // the ID of the first node in the global range
    long global1stNodeID;
    long numberOfNodes;
    long* connectivities;
    double* nodes;
    cvOneDFiniteElement* finiteElement;
    MeshType meshType;

	// minor loss flag / info
    bool isStenosis;
	MinorLoss minorLoss;
	double branchAngle;
	int upstream;
	int branch;

    // Each subdomain can have different materials.
    cvOneDMaterial* mat;

    //  The boundary condition refers to the outlet, either pressure(area) 
    //  specified or resistance specified. As for the inlet, the default
    //  boundary condition is constant flow rate.
     
    BoundCondType  boundType;
    double boundValue;
    double* pressureTime;
    double* pressureWave;
	double* resistanceTime;
	double* resistanceWave;
	double* PressLVWave;
	double* PressLVTime;
	double* presslv;
	int numPressurePts;
	int numPressLVPts;
	
	//RCR BC added IV 050803
	double proximalResistance, capacitance, distalResistance, alphaRCR;
	double rcrTime;
	
    //double MemC(double currP, double previousP, double deltaTime, double currentTime);//for RCR BC -natural, made public to run in Brooke's formulation as well IV 
    double MemD, MemD1, MemD2;//for RCR BC
    double dMemCdP(void);//used for the advective LHS part for RCR BC IV 051303
    double expmDtOne(double deltaTime) {return 1-exp(-alphaRCR*deltaTime);}//for RCR BC IV
    double ConvPressRCR(double previousP, double deltaTime, double currentTime);//used for the pressure convolution in time for RCR BC IV
};

#endif // CVONEDSUBDOMAIN_H
