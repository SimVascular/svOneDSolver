#ifndef USERLPN_H
#define USERLPN_H

#include <string.h>

 
class userLPN{
	public:

		//	double tau;
		//	double HR;
		//	double IC;
		//	double SPV_pre;
		//	double Tpv;

		// extern int USER_DATA = 25;

		int Initialize(int nTimeStep);
		//FindF
		double FindF(double t, double told, double x[], double *f[], double Q[], double P[], int col);
	private:
		// Specify the boundary parameters
		double Tc;
		double Ts;
		double Tr;
		double URA;
		double LTC;
		double RTC;
		double R2TC;
		double fD;
		double fS;
		double LPV;
		double RPV;
		double R2PV;
		double RVV;

		// Parameters for Mynard Valve model. Other local valuables used for local calculation in Mynard Valve modeling are zeta, B, L, Aeff
		double Mst;
		double Mrg;
		double Kvo;
		double Kvc;
		double rho;
		double leff;
		double Aann;

		// Parameters for tricuspid valve. Other local valuables used for local calculation are zeta2, B2, L2, Aeff2
		double Mst2;
		double Mrg2;
		double Kvo2;
		double Kvc2;
		double leff2;
		double Aann2;

		// Parameters of Pressure-Volume Relationships
		double EDPVRc1;
		double EDPVRc2;
		double ESPVRc1;
		double ESPVRc2;


		//Edpvr
		void Edpvr(double v, double *p);
		//Espvr
		void Espvr(double v, double *p);
		//Activate
		void Activate(double t, double Tc, double Ts, double Tr, double *a);
		//ReadUserDat
		void ReadUserData();

};

#endif //USERLPN_H

