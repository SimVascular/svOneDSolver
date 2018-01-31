#ifndef CVONEDCLOSEDLOOP_H
#define CVONEDCLOSEDLOOP_H

# include <vector>
# include <iostream>
# include <fstream>

# include "cvOneDMthModelBase.h"


using namespace std;

class cvOneDClosedLoop{
   
 public:

   void initGenBC(double currentTime, double deltaTime);
   void writeClosedLoop_P(double OneD_P, double currentTime, double deltaTime);
   double getClosedLoop_Q();
   void callGenBC();
   void initGenBC();
 };

 #endif //CVONEDCLOSEDLOOP_H