// General Boundary Conditions for Closed Loop and Coupling Methods
// 	GenBC for 0D/1D coupling: communicates between OneDSolver and
//	a User Defined LPN
// 	Adapted by Melody Dong 10-17-2017 to Mahdi's Fortran GenBC
// 	mldong@stanford.edu 

#include <istream>
#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include "GenBC0D1D.h"
#include "global.h"
#include "userLPN.h"

using namespace std;

void cGenBC::GenBC0D1D() {

	// *****************************
	//     INITIALIZE VARIABLES
	// *****************************
	//surface index, timestep index, final time, # of Dirichlet Surf, # of Neumann Surf
	int i, n, nTimeStep, nDir, nNeu; 
	double t, dt, t_init, delta_t, curr_t;
	double pConv = 1333.224;

	// GenBC Flag: indicate what step in Closed Loop process
	//	flag = I : Initializing
	// 	flag = T : Iteration Loop
	// 	flag = L : Last Iteration
	//	flag = D : Derivative
	//
	//  Here is the period of time that you should integrate for
	//  each of these flags:
	//
	// Flags              I                                         D&T&L
	//                   ^                                          ^
	// 3D code time step: N.........................................N+1
	// 0D code time step: 1.........................................nTimeStep+1
	// Flowrates:         Qi........................................Qf
	// Time, t:         tInitial....................................tInitial+tFinal
	char flag;

	// Create LPN
	userLPN* userManager = new userLPN();


	// ********************************************************
	//       INITIALIZE TIMESTEP AND INITIAL DATA FROM LPM
	// ********************************************************
	int nUnknowns = 5;
	nTimeStep = 0;
	nTimeStep = userManager->Initialize(nTimeStep);




	// *****************************************
	//          READ Q, P FROM 1D SOLVER
	i = nDirichletSrfs;
	
	double *Pi = new double[i];
	// double *Pf = new double[i];
	// PDirichlet = new double[i, 4]; //4 columns wide for RK4

	// For inlet 0D/1D right heart coupling, Q's are placeholders, no values
	// double *Qi = new double[i];
	// double *Qf = new double[i];
	QNeumann = new double[i, 4];

	fstream genBC_file;
	genBC_file.open ("GenBC.int", std::ios_base::in);
	if (genBC_file.fail()){
		cout << "ERROR:Input file does not exist. GenBC.int" << endl;
		//!!!!Put code in to throw an exception
	} else{
		// Flag step in GenBC Process - MELODY SIMPLIFYING TO 1 DIRICHLET SURFACE FOR CLOSED LOOP
		char temp_flag;
		genBC_file >> temp_flag;
		flag = char(temp_flag);

		// Current Timestep
		genBC_file >> delta_t; //1D delta t
		genBC_file >> curr_t; //1D current t

		double p_temp;
		for (i = 0; i < nDirichletSrfs; i++) {
			genBC_file >> Pi[i]; //1D Pressure in cgs
			Pi[i] = Pi[i] / pConv; //0D Pressure handles in mmHg
			// Pf[i] = Pf[i] / pConv;
		}


		
		genBC_file.close(); 
	}
 



	// ***********************************
	//       INITIALIZE THE UNKNOWNS 
	// ***********************************
	// MELODY: took out check for if initializing
	dt = delta_t / double(nTimeStep); //dt for 0D calculations (time discretized more than 1D)

	double RVP;

	double *X = new double[nUnknowns];
	double *Xo = new double[nUnknowns];
	double **f = new double*[4]; //[col][unknown var]
	for (i = 0; i<4; i++) {
		f[i] = new double[nUnknowns];
	}


	for (i = 0; i < nUnknowns; i++) {
		X[i] = 0;
	}

	// Read initial values from Initial Data txtfile
	fstream initialData_file;
	initialData_file.open ("InitialData.txt", std::ios_base::in);
	if (initialData_file.fail()){
		cout << "ERROR:Input file does not exist. InitialData.txt" << endl;
	} else{
		// Read current time step as t_init (will solve for t_init + 1 solution)
		initialData_file >> t_init; 

		// Read Xo initial data left from previous timestep
		for (i = 0; i < nUnknowns; i++) {
			initialData_file >> Xo[i];
		}

		// if not in GenBC Initialization step and is not the 1st iteration of a new timestep
		if (flag == 'T' && curr_t > delta_t) {
			// Read 1st iteration 1D pressure to initialize all 0D solutions of current timestep
			for (i = 0; i < nDirichletSrfs; i++) {
				initialData_file >> Pi[i];
			}
			cout << "	Pressure to 0D: " << Pi[0] << endl;
		}


		initialData_file.close();
	}

	// If new timestep, use initial data created from last iteration of previous timestep
	//	Use pressure from 1D solution at 1st iteration
	fstream tempInit_file;
	if (flag == 'L' && curr_t > delta_t) {
		// Overwrite initial data with previous initial data if moving onto next 1D timestep
		tempInit_file.open("tempInit.txt", std::ios_base::in);
		if (tempInit_file.fail()) {
			cout << "ERROR: Input file does not exist. Temporary initial data file" << endl;
		} else {
			// Time initial: should be equal to current time
			tempInit_file >> t_init;

			for (i = 0; i < nUnknowns; i++) {
				tempInit_file >> Xo[i];
			}

			tempInit_file.close();

			cout << "	Reading from temp init" << endl;

		}

		// Only call and write to initial data file if last iteration signified (i.e moving onto next iter)
		//	rewrite initial data with previous temp initial data
		initialData_file.open("InitialData.txt", std::ofstream::out | std::ofstream::trunc);

		if (initialData_file.fail()){
			cout << "ERROR:Input file does not exist. InitialData" << endl;
		} else{

			// cout << "WRITING TO INITIAL DATA.TXT" << endl;
			// Write current time
			initialData_file << t_init << "\n";

			for (i = 0; i < nUnknowns; i++) {
				initialData_file << Xo[i] << "\n";
			}

			initialData_file.close();
		}

	}


	// ******************************
	//	   CALL USER LPN TO SOLVE
	// ******************************
	//!!! For the sake of debugging, set nTimeStep to 2000 instead of 0 due to GenBC.int has I in it
	// nTimeStep = 3; //MELODY NEED CHANGE
	int col;
	t = t_init;
	for (n = 0; n < nTimeStep; n++){

		// iteration 1
		col = 0;
		userManager->FindF(t, t, Xo, f, QNeumann, Pi, col);
		int row = 0;
		for (row = 0; row < nUnknowns; row++){
			X[row] = Xo[row] + dt * f[0][row] / 3.0;
		}

		// iteration 2		
		col =1;
		userManager->FindF(t + dt / 3.0, t, X, f, QNeumann, Pi, col);
		for (row = 0; row < nUnknowns; row++){
			X[row] = Xo[row] - dt * f[0][row] / 3.0 + dt * f[1][row];
		}

		// iteration 3
		col = 2;
		userManager->FindF(t + dt * 2.0 / 3.0, t, X, f, QNeumann, Pi, col);
		for (row = 0; row < nUnknowns; row++){
			X[row] = Xo[row] + dt * f[0][row] - dt*f[1][row] + dt*f[2][row];
		}

		// iteration 4
		col = 3;
		RVP = userManager->FindF(t + dt, t, X, f, QNeumann, Pi, col);
		for (row = 0; row < nUnknowns; row++){
			f[0][row] = (f[0][row] + 3.0 * f[1][row] + 3.0 * f[2][row] + f[3][row]) / 8.0;
			Xo[row] = Xo[row] + dt * f[0][row];
		}	

		t = t + dt;



	}



	// ************************************
	// 		   WRITE Q,P TO 1D SOLVER
	// ************************************
	for (i = 0; i < nUnknowns; i++){
		X[i] = Xo[i];
	}


	// MELODY MODIFICATIONS - to test pure pulsatile flow coupling
// 	double QPA [] = {103.8738307861, 
// 181.3854096209, 
// 226.4073983107, 
// 252.2668111621, 
// 263.9044913723, 
// 263.2067109456, 
// 251.1089068276, 
// 229.6576547167, 
// 203.0150392149, 
// 172.4661735049, 
// 134.9362452273, 
// 110.6297968659, 
// 97.7282291611, 
// 92.4031871186, 
// 90.6287892118, 
// 90.1301217815, 
// 90.2975260673, 
// 90.9204284319, 
// 99.6005238695, 
// 113.9055067188};
// 	int max = 20;
// 	double t_div = 0.04285;
// 	int ptr;
// 	ptr = int( floor( fmod((delta_t+curr_t)/t_div,max)));
// 	int ptr2 = ptr+1;
// 	if (ptr == 19) {
// 		ptr2 = 0;
// 	}
// 	X[2] = ((QPA[ptr]-QPA[ptr2])/(-t_div))*(fmod(delta_t+curr_t, t_div*max)-t_div*ptr)+QPA[ptr];
//////////////////////////////////////////////////

	double tempQ;
	genBC_file.open("./GenBC.int", std::ofstream::out | std::ofstream::trunc);
			
	if (genBC_file.fail()){
		cout << "ERROR:Input file does not exist. GenBC.int" << endl;
	} else{

		// MELODY MODIFICATIONS
		genBC_file << curr_t << "\n";

		//Writing nDirichlet flowrates here - might need to convert to char to write
		for (i = 0; i < nDirichletSrfs; i++) {
			
			// tempQ = double(X[2]); //MD: X[2] to get Qpv only for USER2 Right heart model coupled to inlet
			genBC_file << X[2] << "\n";

			cout << "	0D qPV: " << X[2] << endl;

		}

		genBC_file.close();
	}
	

	cout << "	0D t_cur: " << curr_t << "\n 	0D t_delta: " << delta_t << endl;
	cout << "  FLAG: " << flag << endl;
	cout << "	1D Press: " << Pi[0] << endl;
	
	cout << "	0D t: " << t << endl;
	cout << "	t_init: " << t_init << endl;
	// Only write to Initial Data file every iteration with GenBC Initialization
	if (curr_t <= delta_t) {
		initialData_file.open("InitialData.txt", std::ofstream::out | std::ofstream::trunc);

		if (initialData_file.fail()){
			cout << "ERROR:Input file does not exist. InitialData" << endl;
		} else{

			//!!!! NOTE Need error checking
			// clear out the file to zero first - why?
			//initialData_file.write(0, sizeof(double));
			cout << "WRITING TO INITIAL DATA.TXT" << endl;

			initialData_file << t_init << "\n";

			for (i = 0; i < nUnknowns; i++) {
				initialData_file << Xo[i] << "\n";
			}

			initialData_file.close();
		}

	}

	// Write 0D and temporary Initial data if new timestep (flag = L)
	fstream results_file;
	if (flag == 'L' && curr_t > delta_t) {
		// Only call and write to initial data file if last iteration signified (i.e moving onto next iter)
		//	rewrite initial data with previous temp initial data
		initialData_file.open("InitialData.txt", std::ofstream::app | std::ofstream::out);
		if (initialData_file.fail()){
			cout << "ERROR:Input file does not exist. InitialData" << endl;
		} else{
			// Append 1D pressure to Initial data file to be called for next time
			initialData_file << Pi[0] << endl;
			initialData_file.close();
		}

		// Write 0D results
		results_file.open("results.dat", std::ofstream::app | std::ofstream::out);
		if (results_file.fail()) {
			cout << "ERROR: Input file does not exist. results.dat" << endl;
		} else {
			results_file << curr_t << "\t";
			for (i = 0; i<nUnknowns; i++) {
				results_file << X[i] << "\t";
			}

			for (i=0; i<nUnknowns; i++) {
				results_file << f[0][i] << "\t";
			}

			results_file << RVP << "\t"; // RV pressure
			results_file << Pi[0] << "\n"; // PA pressure

			results_file.close();

		}
		
	}






	// Creating a temp initial data file to store data after each iteration
	// fstream tempInit_file;
	tempInit_file.open("tempInit.txt", std::ofstream::out | std::ofstream::trunc);
	if (tempInit_file.fail()) {
		cout << "ERROR: Input file does not exist. Temporary initial data file" << endl;
	} else {
		// Write out marched timestep (i.e next timestep)
		tempInit_file << t << "\n";

		// Write out Initial Data from 0D iter: qTC, Vrv, qPV, zeta1, zeta2
		for (i = 0; i < nUnknowns; i++) {
			tempInit_file << X[i] << "\n";
		}

		tempInit_file.close();

	}



	delete[] Pi;
	// delete[] Pf;
	// delete[] PDirichlet;
	// delete[] Qi;
	// delete[] Qf;
	delete[] QNeumann;
	delete[] X;
	delete[] Xo;
	delete[] f;

	delete userManager;


}