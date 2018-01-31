#ifndef CVONEDMORPHIMPEDANCE_H
#define CVONEDMORPHIMPEDANCE_H

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "cvOneDComplex.h"

class cvOneDMorphImpedance{
    
  public:

    double r_root; // root (starting) radius of tree
    
    cvOneDComplex** Computed;
    double *z_t;
    double *y_t;
    
    cvOneDComplex Complex_ZERO;
    cvOneDComplex Complex_ONE;
    cvOneDComplex Complex_I;
    int lengthN;
      
    cvOneDMorphImpedance();

    ~cvOneDMorphImpedance();

    //ImpedanceTime returns the impedance in the time domain
    double*  calculateImpedanceTime(double rootR, int min_order,double Period, 
                                       int tmstps);
    //ImpedanceFourier returns the real component of Z in the frequency domain
    double*  calculateImpedanceFourier(double rootR, int min_order,double Period, 
                                       int tmstps);
    //AdmittanceTime returns admittance in the time domain
    double*  calculateAdmittanceTime(double rootR, int min_order,double Period, 
                                       int tmstps);
    //AdmittanceFourier returns the real part of the admittance in the frequency domain
    double*  calculateAdmittanceFourier(double rootR, int min_order,double Period, 
                                       int tmstps);
    double*  calculateImpAdmit(double rootR, int min_order, double Period, int tmstps,
                               bool returnImp, bool returnFourier);
    
  private:
  
    // Take inverse Fourier transform of array
    cvOneDComplex* InverseFourier(cvOneDComplex* Y);

    // Fourier isn't really used, but figured I'd leave it just in case (nate)
    cvOneDComplex* Fourier(cvOneDComplex* y);
  
}; // end of class

#endif // CVONEDMORPHIMPEDANCE_H
