#ifndef CVONEDMETTEIMPEDANCE_H
#define CVONEDMETTEIMPEDANCE_H

//
// Calculate the Impedance for a 1D vascular network.
//

#include <stdio.h>
#include <math.h>
#include <iostream>

#include "cvOneDComplex.h"

using namespace std;

class cvOneDMetteImpedance{    

  public:

    // variables
    bool METTE; // Use mette's structured tree method to compute impedance, will turn on if lrr<1

    bool EXERCISE; // Modify vessel radii wilth this flag
    double ex_fact; // If >1 dilate, if <1 constrict

    double r_min; // Minimum radius of tree
    double r_root; // Root (starting) radius of tree
    double alpha,beta; // Scaling parameters

    bool F90;
    bool F95; // Use dimensionless parameters: true, use f95 version, if false use f90(old)version
    double g; // Used with F95, gravity
    double Lr; // Used with F95, dimensionless length
    double q; // Used with F95, dimensionless flow
    double Fr2; // Used with F95, Froude number?

    double lrr; // Length to radius ratio mette = 50
                //Zamir: average=20 max=70 (for carciac vessels)

    double mu; // viscosity, cp = g/cm/s( mu = 0.035??)
    double rho; // density  g/cm^3

    double k1; // constants for compliance
    double k2; // constants for compliance
    double k3; // constants for compliance

    bool RCR; // use RCR impedance
    double R1; // for RCR curciut computation. (defaults?)
    double R2; // for RCR curciut computation. (defaults?)
    double C1; // for RCR curciut computation. (defaults?)

    int Maxgen; // The number of generations in the structured tree
   
                             // The asymmetry ratio of the structured tree
    double asym ;  //I tried changing this to .85, 34 generations
    double asym1;
    double asym2;

    double expo;     // Exponent in radius relation. radius of daughter
                            // vessels to parent
                            // Miller (1893) deduced this value for the larger vessels
                            // of the lungs. 2.51 was the median value found by 
                            // Papagergiou of measured vessels greater than 2 mm 
                            // in diameter

    double radius_flag1;// 500 micron diameter   500 micron radius
    double expo1;        // Suwa (1963)Miller(1893)Hooper(1977) proposed a 
                                // number of values for the middle vessels. 
                                // This is approximate
    double radius_flag2;// 100 micron diameter   250 micron radius
    double expo2;        // Bu Murray's law radius at which different
                                // exponent kicks in. we want to be able to take min
                                // radius down to the level of the capilaries.  
                                // higher exponents during exercise conditions will 
                                // increase this number, so use  the pre-computed number
                                // of generations to determine new end condition.
    double radius_flag0; // 500 micron diameter
    double ex_expo;      //3.2

    cvOneDComplex** Computed;
    double *z_t;
    double *y_t;

    int localmax;                // max numbver of generations in a tree
    int Z_out;
    int z_out;
    int y_out;
    int writeout;       // flag to writeout Z, Y ,z, or y

    cvOneDComplex Complex_ZERO;
    cvOneDComplex Complex_ONE;
    cvOneDComplex Complex_I;
    int lengthN;
      
    cvOneDMetteImpedance();

    ~cvOneDMetteImpedance();

    //ImpedanceTime returns the impedance in the time domain
    double*  calculateImpedanceTime(double rootR, double lrr,double Period, 
                                       int tmstps);
    //ImpedanceFourier returns the real component of Z in the frequency domain
    double*  calculateImpedanceFourier(double rootR, double lrr,double Period, 
                                       int tmstps);
    //RCRImpedanceTime returns the impedance in the time domain
    double*  calculateRCRImpedanceTime(double Rd, double Ru,double Cap,
                               double Period,int tmstps);
   //RCRImpedanceFourier returns the real part of the impedance in the frequency domain
    double*  calculateRCRImpedanceFourier(double Rd, double Ru,double Cap,
                               double Period,int tmstps);
    //AdmittanceTime returns admittance in the time domain
    double*  calculateAdmittanceTime(double rootR, double lrr,double Period, 
                                       int tmstps);
    //AdmittanceFourier returns the real part of the admittance in the frequency domain
    double*  calculateAdmittanceFourier(double rootR, double lrr,double Period, 
                                       int tmstps);
    //RCRAdmittanceTime returns admittance in the time domain
    double*  calculateRCRAdmittanceTime(double Rd, double Ru,double Cap,
                               double Period,int tmstps);
    //RCRAdmittanceFourier returns the real part of the admittance in the frequency domain
    double*  calculateRCRAdmittanceFourier(double Rd, double Ru,double Cap,
                               double Period,int tmstps);
  //kimhj 08-18-2004 added
    //resist saves the input, DC impedance
    double resist;
    double resistBasedTree( double r_root, double resist, double Period, int tmstps);
    
    /******************************************************************************
     * This will calculate the Impedance                                          *   
     * rootR - radius of parent vessel (we will attach fractal tree to this)      *
     * minR - minimum radius of vessel in structured tree                         *
     * asym - asymmetry variable, default value is 0.4048                         *
     * Period - period of cardiac cycle                                           *
     * tmstps - number of timesteps to solve for (I like to solve for 36)         *
     * resist - resistance at end of structured trees?  Maybe overall resistance  *
     * I excpec structured tree to have                                           *
     *                                                                            *
     * returns admittance or impedance in the time or frequency domain            *
     *                                                                            *
     * change this a little, use consistant r_min = 5                             *
     * microns = 0.0005, vary length to radius ratio only?,                       *
     * resistance always zero at bc                                               *
     ******************************************************************************/

    double*  calculateImpAdmit(double rootR, double lrr, double asym,
                                       double Period, int tmstps, bool returnImp, bool returnFourier);

    cvOneDComplex* comp_imp(int N,double* Omega);

    void SetExerciseFlag(bool flag);
    bool GetExerciseFlag();

    void SetExerciseFactor(double factor);
    double GetExerciseFactor();
    
  private:
  
    cvOneDComplex Z0func (double omega_k, int alpha_pow,int beta_pow);

    cvOneDComplex RCRImpedance(double R, double Rd, double C, double omega);

    cvOneDComplex getF10(double wom);
    
    // Take inverse Fourier transform of array
    cvOneDComplex* InverseFourier(cvOneDComplex* Y);

    // Fourier isn't really used, but figured I'd leave it just in case (nate)
    cvOneDComplex* Fourier(cvOneDComplex* y);

    //kimhj 08-18-2004 added
    //the function "calculate LRR" estimates the lrr when given the resist
    //and root radius using the approximation that asym=asym1 and expo=expo1
    // and doesn't change as the radius changes
    double calculateLRR(double rootR, double resist);

}; // end of class

#endif // CVONEDMETTEIMPEDANCE_H
