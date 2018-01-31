# include "cvOneDClosedLoop.h"
#include <iostream>
#include <fstream>
#include <ostream>
#include <stdio.h>

// cvOneDClosedLoop: subroutines that facilitate passing Q, P to and from 
// a 0D lumped parameter model

// Initializes simulation with constant steady flow to fill Initial Data for GenBC
void cvOneDClosedLoop::initGenBC(double currentTime, double deltaTime){
  // hardcoded distal PA pressure - MELODY needs to change
  double P = 11997;
  // initialize with a stable PA pressure
  writeClosedLoop_P(P, currentTime, deltaTime); 
  // fill initialdata.txt with known stable initializer
  callGenBC();
}


// Writes pressure from 1D model to external file (GenBC.int) for 0D model to access
void cvOneDClosedLoop::writeClosedLoop_P(double OneD_P, double currentTime, double deltaTime) {

  char flag = 'T';

  // extra info being passed - to be implemented more with neumann coupling and outlet closed loop
  int nNeu = 0; // not used when passing info to 0D model yet - will need to be implemented when using for Neumann coupling
  int nDir = 1; // not used when passing info to 0D model yet - will need to be implemented when using multiple surfaces for coupling
  double t_prev;

  std::fstream genBCInFile;

  // Check if moved onto next iteration - if so, change flag to L
  genBCInFile.open("./GenBC.int", std::ios_base::in);
  genBCInFile >> t_prev;
  cout << "  previous t: " << t_prev << endl;
  cout << "  current t: " << currentTime << endl;

  if (currentTime > double(t_prev)+1e-5 || currentTime < double(t_prev)-1e-5){
    flag = 'L';
    cout << "	last iteration previously: " << t_prev << endl;
  }
  genBCInFile.close();

  // Write 1D values to GenBC
  genBCInFile.open ("./GenBC.int", std::ofstream::out | std::ofstream::trunc);
  if(genBCInFile.is_open()){
		
    // MELODY: for simple case
    genBCInFile << flag << "\n";
    genBCInFile << deltaTime << "\n";
    genBCInFile << currentTime << "\n";
    genBCInFile << OneD_P << "\n";

  }
  genBCInFile.close();
}

// Reads 0D flow from external file (GenBC.int) that 0D model writes to
double cvOneDClosedLoop::getClosedLoop_Q() {

  double Qf;
  double deltaTime;

  std::fstream genBCInFile;
  genBCInFile.open ("./GenBC.int", std::ios_base::in);
  if(genBCInFile.is_open()){
    genBCInFile >> deltaTime;
    genBCInFile >> Qf;
    // std::cout << "Q from GenBC: " << Qf << "\n";
  }
  genBCInFile.close();

  return Qf;
}

// Calls GenBC and UserLPN executable from current working directory
void cvOneDClosedLoop::callGenBC() {

	int a = system("\"./genBC0D1D.exe\"");

}