// User Defined Lumped Parameter Model for 0D/1D Coupling
// 	0D Right Heart LPM - modified from Weiguang Yang's LPM
// 	(initial adaptation from Mynard et al's model)
// 	Adapted by Melody Dong 9-11-2017 to Mahdi's GenBC, User2.f
//	files.
// 	mldong@stanford.edu

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "global.h"
#include "userLPN.h"

using namespace std;

// bool pCorr;
// bool qCorr;
int nUnknowns;
int nDirichletSrfs;
int nNeumannSrfs;
// int nXprint;
double pConv;
double qConv;
// int *srfToXPtr;
// int *srfToXdPtr;
double *QNeumann;
double *PDirichlet;
// double *offset;
// double *Xprint;


int userLPN::Initialize(int nTimeStep) {
	bool ierr = false;
	int i;

	//********************************************************************
	// For instance if pressure in 1D solver is in cgs and here mmHg
	// pConv = 1334 = 1.334D3, also if flowrate in 3D solver is mm ^ 3 / s and
	// here is mL / s qConv = 1000 = 1D3. In the case both solver are using same
	// unites you can set these two conversion coefficients to 1D0
	// pConv = 1.333e+3;
	// qConv = -1.0;

	// Only when all the surfaces of you model are coupled with NeumannSrfs
	// you may set this to.TRUE.
	// pCorr = false;
	// qCorr = false;

	//********************************************************************
	// Block for your inputs

	// MELODY MODIFICATIONS: include nDir and Neumann surfaces, unknowns timestep, etc. into external users data file

	// These two value should match to that of defined in solver.inp
	nDirichletSrfs = 1;
	nNeumannSrfs = 0;

	// Number of unknowns that you need inside your lumped parameter network
	nUnknowns = 5;
	// Number of time step between N and N + alpha, for higher values this code
	// would be more accurate and more costly
	nTimeStep = 2000;

	// Number of parameters to be printed in AllData file(the first
	// nUnknowns columns are by default X)
	// nXprint = 4;

	ReadUserData();

	//--------------------------------------------------------------------
	// You don't need to change this part      
	double *tZeroX = new double[nUnknowns];
	//--------------------------------------------------------------------
	
	for (i = 0; i < nUnknowns; i++){
		tZeroX[i] = 0.0;
	}

	//--------------------------------------------------------------------
	// Value of your unknown at time equal to zero(This is going to be used
	// ONLY when you are initiating simulation)

	tZeroX[1] = RVV; //1.23e+2;
	// tZeroX[5] = 

	fstream initialData_file("InitialData.txt");
	
	if (initialData_file.fail()){
		cout << "ERROR:Initial Data file does not exist." << endl;
		cout << "Creating new initData file" << endl;

		initialData_file.open("InitialData.txt", std::ofstream::out | std::ofstream::trunc);

		// Write initial timestep - hardcoded to 0
		initialData_file << tZeroX[0] <<"\n";

		// Write initialized unknowns
		for (i = 0; i < nUnknowns; i++) {
			initialData_file << tZeroX[i] << "\n";
		}

		initialData_file.close();
	}

	delete[] tZeroX;

	return nTimeStep;
}


//####################################################################
// Here you should find the f_i = dx_i / dt, based on the following parameters :
//  current x_i : x(i)
//  Current time : t
//  Flowrates of Neumann faces : Q(i)
//  Pressures of Dirichlet faces : P(i)

double userLPN::FindF(double t, double told, double x[], double *f[], double Q[], double P[], int col){
	double RVP; //actual RV Pressure
	double EDPVR_val = 0; //diastolic RV Pressure
	double ESPVR_val = 0; //systolid RV Pressure
	double a = 0;
	// P[0] = abs(P[0]);

	// Pulmonary Valve Opening (zeta1): 0=closed, 1=open
	if (x[3] > 1.0) {
		x[3] = 1.0;
	}
	if (x[3] < 1.0e-8){
		x[3] = 0;
	}

	double zeta = x[3];

	// Tricuspid Valve Opening (zeta2): 0=closed, 1=open
	if (x[4] > 1.0) {
		x[4] = 1.0;
	}
	if (x[4] < 1.0e-8) {
		x[4] = 0.0;
	}

	double zeta2 = x[4];

	// Calculate Diastolic and Systolic RV pressure based on PV relationship
	Edpvr(x[1], &EDPVR_val);
	Espvr(x[1], &ESPVR_val);

	// cout << "RV Volume: " << x[1] << endl;
	// cout << "	EDPVR_val: " << EDPVR_val << endl;

	// Ventricular Activation, a: 0 to 1 (relaxing vs. contracting)
	Activate(t, Tc, Ts, Tr, &a);

    // calculate RV pressure
	RVP = (1 - a)*fD*EDPVR_val + a*fS*ESPVR_val;

	// cout << "activation: " << a << endl;
	// cout << "RV Pressure: " << RVP << endl;
	// cout << "MPA Pressure: " << P[0] << endl;
	// cout << "zeta1: " << x[3] << endl;

    // for pulmonary valve: effective area
	double Aeff = (Mst*Aann - Mrg*Aann)*zeta + Mrg*Aann;
	double B = rho / (2.0*Aeff*Aeff);
	double L = rho*leff / Aeff;

    // for tricuspid valve: effective area
	double Aeff2 = (Mst2*Aann2 - Mrg2*Aann2)*zeta2 + Mrg2*Aann2;
	double B2 = rho / (2.0*Aeff2*Aeff2);
	double L2 = rho*leff2 / Aeff2;

	// Partial derivatives of X values
	// Flow through tricuspid valve
	if (Aeff2 < 1.0e-6) {
		f[col][0] = 0.0;
	}
	else {
		f[col][0] = (URA - RVP - B2*(x[0])*abs(x[0])) / L2;
	}

	// RV Volume derivative = flow through tricuspid - flow through pulmonary
	f[col][1] = x[0] - x[2];

	// Flow through pulmonary valve
	if (Aeff < 1.0e-6){
		f[col][2] = 0.0;
	}
	else{
		f[col][2] = (RVP - P[0] - B*(x[2])*abs(x[2])) / L;
	}

    // for  dynamics of pulmonary valve
	if (RVP > P[0]) {
		f[col][3] = (1 - zeta)*Kvo*(RVP - P[0]);
	}
	else{
		f[col][3] = zeta*Kvc*(RVP - P[0]);
	}

	// for dynamics of tricuspid valve
	if (URA > RVP) {
		f[col][4] = (1 - zeta2)*Kvo2*(URA - RVP);
		// cout << "f[4]: " << &f[4,col] << endl;
		// cout << "RA Press > RV Press" << endl;
	}
	else{
		f[col][4] = zeta2*Kvc2*(URA - RVP - 0.5); //from 3D USER2.f: - 0.5);
	}


	return RVP;	

}

////////////////////////////////////////////////////////////////////////////////
//
//   end diastolic pressure volume relastionship derived from scaling Levick Figure 7.18a
//

void userLPN::Edpvr(double v, double *p){

	*p = ((1.1435e-5 * v - 2.501e-3)*0.0 * v + EDPVRc1) * v + EDPVRc2;
}


///////////////////////////////////////////////////////////////////////////////
//   end systolic pressure volume relastionship derived from scaling Levick Figure 7.18a
//

void userLPN::Espvr(double v, double *p){

	*p = ((3.157e-6 * v - 1.6534e-3)*0.0 * v + ESPVRc1) * v + ESPVRc2;
}

////////////////////////////////////////////////////////////////////////////
//	Ventricular activation 
//
void userLPN::Activate(double t, double Tc, double Ts, double Tr, double *a){

	double pi;
	double tm;
	double AA;
	pi = 3.14159265358979323846;

	tm = fmod(t,Tc);
	AA = 2.0 / (pi + 2.0);
	if (tm < (Ts / 2.0)){
		*a = (1.0 - AA) * 2.0 * (tm / Ts);
	}
	else if((Ts/2.0) <= tm && tm < Ts) {
		*a = 1.0 - AA + AA*sin(pi*((tm / Ts) - 0.5));
	}
	else if (Ts <= tm && tm < (Ts + Tr)){
		*a = 1.0 - sin((pi / 2.0)*(tm - Ts) / Tr);
	}
	else {
		*a = 0.0;
	}

	
}


////////////////////////////////////////////////////////////////////////////
void userLPN::ReadUserData(){
	int USER_DATA = 30;


	std::fstream user_file;

	user_file.open ("user_param.dat", std::ios_base::in);
	if (user_file.fail()){
		cout << "ERROR:Input file does not exist. user_param.dat" << endl;
		// Put code in to throw an exception
	}
	else{
		std::string buffer;

		for (int i = 0; i < USER_DATA; i++) {
			std::string flag;
			user_file >> flag;

			if (flag == "TC"){
				// double Tc;
				user_file >> Tc;
			} 
			else if (flag == "TS"){
				// double Ts;
				user_file >> Ts;
			}
			else if (flag == "TR"){
				// double Tr;
				user_file >> Tr;
			}
			else if (flag == "URA"){
				// double URA;
				user_file >> URA;
			}
			else if (flag == "LTC"){
				// double LTC;
				user_file >> LTC;
			}
			else if (flag == "RTC"){
				// double RTC;
				user_file >> RTC;
			}
			else if (flag == "R2TC"){
				// double R2TC;
				user_file >> R2TC;
			}
			else if (flag == "FD"){
				// double fD;
				user_file >> fD;
			}
			else if (flag == "FS"){
				// double fS;
				user_file >> fS;
			}
			else if (flag == "LPV"){
				// double LPV;
				user_file >> LPV;
			}
			else if (flag == "RPV"){
				// double RPV;
				user_file >> RPV;
			}
			else if (flag == "R2PV"){
				// double R2PV;
				user_file >> R2PV;
			}
			else if (flag == "RVV") {
				// double RVV;
				user_file >> RVV;
			}
			else if (flag == "MST"){
				// double Mst;
				user_file >> Mst;
			}
			else if (flag == "MRG"){
				// double Mrg;
				user_file >> Mrg;
			}
			else if (flag == "KVO"){
				// double Kvo;
				user_file >> Kvo;
			}
			else if (flag == "KVC"){
				// double Kvc;
				user_file >> Kvc;
			}
			else if (flag == "RHO"){
				// double rho;
				user_file >> rho;
			}
			else if (flag == "LEFF"){
				// double leff;
				user_file >> leff;
			}
			else if (flag == "AANN"){
				// double Aann;
				user_file >> Aann;
			}
			else if (flag == "MST2"){
				// double Mst2;
				user_file >> Mst2;
			}
			else if (flag == "MRG2"){
				// double Mrg2;
				user_file >> Mrg2;
			}
			else if (flag == "KVO2"){
				// double Kvo2;
				user_file >> Kvo2;
			}
			else if (flag == "KVC2"){
				// double Kvc2;
				user_file >> Kvc2;
			}
			else if (flag == "LEFF2"){
				// double leff2;
				user_file >>leff2;
			}
			else if (flag == "AANN2"){
				// double Aann2;
				user_file >> Aann2;
			}
			else if (flag == "EDPVRc1"){
				// double EDPVRc1;
				user_file >> EDPVRc1;
			}
			else if (flag == "EDPVRc2") {
				// double EDPVRc2;
				user_file >> EDPVRc2;
			}
			else if (flag == "ESPVRc1") {
				// double ESPVRc1;
				user_file >> ESPVRc1;
			}
			else if (flag == "ESPVRc2") {
				// double ESPVRc2;
				user_file >> ESPVRc2;
			}
			else {
				break;
			}

		}
		user_file.close();
	}

}
