#ifndef CVONEDSUBDOMAIN_H
#define CVONEDSUBDOMAIN_H

//
//  cvOneDSubdomain.h - Header for a class to contain the descritization of
//  ~~~~~~~~~~~~~   of the Geometry.
//
//  History:
//  Aug. 8, 2003, I. Vignon
//	    Added MemIntWave, MemAdvWave, dMemIntWavedP, dMemAdvWavedP for Wave boundary conditions
//      Added the wave boundary condition
//  May,16, 2003, I. Vignon
//      Added MemIntImp, MemAdvImp, dMemAdvImpdP for impedance boundary conditions
//  May. 8, 2003, I. Vignon
//	    Added MemIntRCR, MemAdvRCR, dMemIntRCRdP, dMemAdvRCRdP for RCR boundary conditions
//  Mar. 26, 2003, I. Vignon
//      Added the RCR boundary condition
//  Oct., 2000 J.Wan
//      Added the Pressure wave boundary condition.
//  May 1999, J.Wan, S.A.Spicer and S.Strohband
//      Creation of file,

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
	void SetBoundImpedanceValues(double *h, int num);
	void SetBoundCoronaryValues(double *time, double *p_lv, int num); //added kimhj 09022005
	void SetBoundRCRValues(double *rcr, int num);//added IV 050803
	void SetBoundWaveValues(double *wave, int num);//added IV 080603
    void SetBoundResistPdValues(double *value, int num); //added wgyang 2019/4

	double GetBoundArea(){return boundValue;}
    double GetBoundResistance(){return boundValue;}
    double GetBoundResistance(double currentTime);
	double getBoundCoronaryValues(double currentTime);
    double GetBoundFlowRate(){return boundValue;}
    double GetBoundAreabyPresWave(double currentTime);
	double GetPressure(double currentTime);
	double GetdPressuredt(double currentTime) { return 0; }

	double* GetImpedance(){return impedance;}
	double* ShiftPressure(double lastP, double currentTime, double currS);
	//double* FlipAndShiftFlow(double lastQ, double currentTime);
	double*	GetDpDs(double lastDpDs);// called AFTER FlipAndShiftPressure
	int		GetNumImpedancePts(){return numImpedancePts;}

	double GetRp(){return proximalResistance;}//added IV 050803
	double GetCap(){return capacitance;}//added IV 050803
	double GetRd(){return distalResistance;}//added IV 050803
	double GetAlphaRCR(){return alphaRCR;}//added IV 050803

    double GetResistanceR(){return resistancevalue;}//added wgyang 2019/4
    double GetResistancePd(){return Pd;}//added wgyang 2019/4

	double GetRa1(){return Ra1;}
	double GetRa2(){return Ra2;}
	double GetCa(){return Ca;}
	double GetCc(){return Cc;}
	double GetRv1(){return Rv1;}
	double GetRv2(){return Rv2;}

    //double* GetEigValWave(){return eigValWave;}//added IV 080703

	/*
	* Use the following functions MemIntRCR, MemAdvRCR, dMemIntRCRdP, dMemAdvRCRdP for RCR boundary conditions
	* Mem is for "memory" - time integrals
	* added IV 050803
	*/
	double MemIntRCR(double currP, double previousP, double deltaTime, double currentTime);//compute integral over one time step of the pressure convolution in time int(P(t')*exp(-alphaRCR(t-t')),0,t)
	double MemIntPexp(double previousP, double deltaTime, double currentTime);//compute integral over one time step of the pressure convolution in time int(P(t')*exp(-alphaRCR(t-t')),0,t) wgyang 2019/4
	double MemIntSexp(double previousS, double deltaTime, double currentTime);//compute integral over one time step of the area convolution in time int(A(t')^(-3/2)*exp(-alphaRCR(t-t')),0,t) wgyang 2019/4
	double MemAdvRCR(double currP, double previousP, double deltaTime, double currentTime);//compute part of the advective term integral over one time step of Q^2, Q=couplingFunction(P)
	double dMemIntRCRdP(double deltaTime);//used for contribution to LHS
	double dMemAdvRCRdP(double currP, double previousP, double deltaTime, double currentTime);//contribution to LHS: compute part of the advective term integral over one time step of Q^2, Q=couplingFunction(P))


	double MemC(double currP, double previousP, double deltaTime, double currentTime);//if RCR essential-see ApplyBoundaryConditions()
	/* Use the following functions MemIntCoronary, TotalMemIntCoronary, MemAdvCoronary, dMenIntCoronarydP,
	* dTotalMemIntCoronarydP, dMemAdvCoronarydP, MemCoronary1, MemCoronary2 for coronary boundary conditions
	* added kimhj 09022005
	*/
	double CORic1(void);
	double CORic2(void);
	double MemIntCoronary(double currP, double previousP, double deltaTime, double currentTime, double exponent);
    double MemAdvCoronary(double currP, double previousP, double deltaTime, double currentTime);
	double dMemIntCoronarydP(double deltaTime, double exponent);
    double dMemAdvCoronarydP(double currP, double previousP, double deltaTime, double currentTime);
	double MemCoronary1(double currP, double previousP, double deltaTime, double currentTime);
	double MemCoronary2(double currP, double previousP, double deltaTime, double currentTime);
    /*
	* Use the following functions MemIntWave, MemAdvWave, dMemIntWavedP, dMemAdvWavedP for wave boundary conditions
	* Mem is for "memory" - time integrals
	* added IV 080603
	*/
	double MemQWave(double currS, double prevS, double initS, double deltaTime,double currentTime);//compute integral of Q
	double MemAdvWave(double currS, double prevS,double initS,double deltaTime,double currentTime);//compute part of the advective term integral over one time step of Q^2, Q=couplingFunction(P)
	double dMemQWavedS(double deltaTime, double currentTime);//used for contribution to LHS
	double dMemAdvWavedS(double currS, double prevS,double initS,double deltaTime, double currentTime);//contribution to LHS: compute part of the advective term integral over one time step of Q^2, Q=couplingFunction(P))

    /*
	* Use the following functions MemIntImp, MemAdvImp, dMemAdvImpdP for impedance boundary conditions
	* Mem is for "memory" - time integrals
	* added IV 051603
	*/
	double MemIntImp(double *press, double deltaTime, double currentTime);
	double MemAdvImp(double *press, double deltaTime, double currentTime);
	double dMemAdvImpdP(double *press, double deltaTime, double currentTime);
	double dMemIntImpdP(double deltaTime);

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
	double* impedance;
	double* resistanceTime;
	double* resistanceWave;
	double* PressLVWave;
	double* PressLVTime;
	double* presslv;
//	double* impedanceFlow;
	double* impedanceDpDs;
	int numPressurePts;
	int numImpedancePts;
	int numPressLVPts;
	double impedanceTime;

	//RCR BC added IV 050803
	double proximalResistance, capacitance, distalResistance, alphaRCR;
	double rcrTime;
    //resistance with Pd wgyang 2019/4
    double resistancevalue;
    double Pd;
    double rcrTime2;
    double rcrTime3;

	//Coronary BC kimhj 08312005
	double p0COR, p1COR, p2COR;
	double q0COR, q1COR, q2COR;
	double b0COR, b1COR;
	double expo1COR, expo2COR, detCOR;
	double CoefZ1, CoefY1, CoefZ2, CoefY2, CoefR;
	double Ra1, Ra2, Ca, Cc, Rv1, Rv2;
	double corTime;

    //wave BC added IV 080603, parameters for the downstream tube
    //reference corss sectional area, number of modes in the "infinite" sum, domain length, end S BC, viscosity coeff, damping factor, wave speed
	double waveSo, numWaveMod, waveLT, waveEndh, waveN, waveAlpha, waveSpeed;
	double waveTime;
    int numWavePts;
    double* eigValWave;

    //double MemC(double currP, double previousP, double deltaTime, double currentTime);//for RCR BC -natural, made public to run in Brooke's formulation as well IV
    double MemD, MemD1, MemD2;//for RCR BC
    double MemConvP;//convolution P(t)*exp(-alpharcr*t) for RCR BC wgyang 2019/4 test!
    double MemConvS;//convolution S(t)^(-3/2)*exp(-alpharcr*t) for RCR BC wgyang 2019/4 test!
	double MemDImp;//for impedance BC
    double MemDWave, MemDWave1, MemDWave2;//for Wave BC
    complex<double> MemEIWaveM;//for Wave BC
    complex<double> MemEIWaveM1;//for Wave BC
    complex<double>* MemEIw;//for Wave BC
    complex<double>* MemEIw1;//for Wave BC
    complex<double>* MemEIw2;//for Wave BC
    complex<double>* dMemEIw;//for Wave BC
    double dMemCdP(void);//used for the advective LHS part for RCR BC IV 051303
    double expmDtOne(double deltaTime) {return 1-exp(-alphaRCR*deltaTime);}//for RCR BC IV
    double ConvPressRCR(double previousP, double deltaTime, double currentTime);//used for the pressure convolution in time for RCR BC IV
    double ConvPressImp(double *press, double currentTime);//used for the pressure convolution with admittance in Impedance BC IV 051603
    double ConvPressCoronary(double previousP, double deltaTime, double currentTime, double exponent); //added kimhj 09022005
    double expmDtOneCoronary(double deltaTime, double exponent); //added kimhj 09022005
    double ConvPressexp(double previousP, double deltaTime, double currentTime);//convolution pressure = int(P(t')exp(-alpharcr(t-t')) in time for RCR BC wgyang test!
    double ConvSexp(double previousS, double deltaTime, double currentTime);//convolution S^(-3/2) = int(S(t')^(-3/2)exp(-alpharcr(t-t')) in time for RCR BC wgyang test!
    double dMemCoronary1dP(void);
    double dMemCoronary2dP(void);
    double WaveEndBC, WaveIni, MemEWave, MemEIWave;//for wave BC
    complex<double>* MemEISpWave;//for wave BC 1st term in double integral convolution
    complex<double>* MemEISmWave;//for wave BC 1st term in double integral convolution
    complex<double>* MemEIS2pWave;//for wave BC 2nd term in double integral convolution
    complex<double>* MemEIS2mWave;//for wave BC 2nd term in double integral convolution
    double previous2S;// need to keep track of S previous to previousS for wave BC
    double memIniWave(double initS,double Time);//compute the IC contribution in Q at time t
    double memEndBCWave(double Time);//compute the endBC contribution in Q at time t
    double convolWave(double currS, double previousS, double deltaTime, double currentTime, double Time);//compute the single integral convolution in Q at time t
    double dblConvolWave(double currS, double previousS, double deltaTime, double currentTime, double Time);//compute the double integral convolution in Q at time t
    double QWave(double currS, double prevS, double initS, double deltaTime, double currentTime,double Time);//compute Q at time t in time slab tn to tn+1
};

#endif // CVONEDSUBDOMAIN_H
