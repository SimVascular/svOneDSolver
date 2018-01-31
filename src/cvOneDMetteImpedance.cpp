/*This routine computes the admittance and impedance in the Fourier and time domain
* following either:
* 1)Womersley theory applied to a structured tree (see Olufsen's and B. Steele's PhD)
* 2)or the RCR model (electric analog - 0D model)
* It also modifies them to accomodate for exercise if necessary (scaling factor- see B. Steele's thesis)
* 
*
* Author: B.N.Steele
* Converted from java to C by N. Wilson 2004
* I. Vignon added output z(t) 08/06/2004 
*/

# include <ostream>

# include "cvOneDGlobal.h"
# include "cvOneDMetteImpedance.h"

void cvOneDMetteImpedance::SetExerciseFlag(bool flag){
  EXERCISE = true;
}

bool cvOneDMetteImpedance::GetExerciseFlag(){
  return EXERCISE;
}

void cvOneDMetteImpedance::SetExerciseFactor(double factor){
  ex_fact = factor;
}

double cvOneDMetteImpedance::GetExerciseFactor(){
  return ex_fact;
}

cvOneDMetteImpedance::cvOneDMetteImpedance(){

  // Variables Initialization
  METTE = false; // Use mette's structured tree method to compute impedance, will turn on if lrr<1

  EXERCISE = false; // modify vessel radii wilth this flag
  ex_fact = 2.0; // if >1 dilate, if <1 constrict

  r_min = 0.0003; // Minimum radius of tree
  r_root = 0; // Root (starting) radius of tree
  alpha = 0,beta = 0; // Scaling parameters

  F90 = true;
  F95 = false;     // use dimensionless parameters.
                           // true, use f95 version, if false use f90(old)version
  g = 981.0; // Used with F95, gravity in CGS
  Lr = 1.0; // Used with F95, dimensionless length
  q = 10.0*Lr*Lr; // Used with F95, dimensionless flow
  Fr2= q*q/g/Lr/Lr/Lr/Lr/Lr; // Used with F95, Froude number??

  lrr = 20.0; // length to radius ratio mette = 50
              //Zamir: average=20 max=70 (for carciac vessels)

  mu  = 0.049; // Viscosity, cp = g/cm/s( mu = 0.035??)
  rho = 1.055; // Density  g/cm^3

  k1 = 2.0E7; // Constants for compliance
  k2 = -22.5267; // Constants for compliance
  k3 = 8.65E5; // Constants for compliance

  RCR = false; // Use RCR impedance
  R1 = 1200.0; // For RCR circuit computation. (defaults?)
  R2 = 90.0; // For RCR circuit computation. (defaults?)
  C1 = 0.00008; // For RCR circuit computation. (defaults?)

  Maxgen = 60; // The number of generations in the structured tree
               // The asymmetry ratio of the structured tree
  asym  = 0.4048 ;  // I tried changing this to .85, 34 generations
  asym1 = 0.6;
  asym2 = 0.9;

  expo = 2.51; // Exponent in radius relation. radius of daughter vessels to parent
               // Miller (1893) deduced this value for the larger vessels
               // of the lungs. 2.51 was the median value found by
               // Papagergiou of measured vessels greater than 2 mm
               // in diameter.

  radius_flag1 = 0.025;// 500 micron diameter   500 micron radius
  expo1 = 2.76;        // Suwa (1963)Miller(1893)Hooper(1977) proposed a 
                              // number of values for the middle vessels. 
                              // This is approximate
  radius_flag2 = 0.005;// 100 micron diameter   250 micron radius
  expo2  = 2.9;        // Bu Murray's law radius at which different
                              // exponent kicks in. we want to be able to take min
                              // radius down to the level of the capilaries.  
                              // higher exponents during exercise conditions will 
                              // increase this number, so use  the pre-computed number
                              // of generations to determine new end condition.
  radius_flag0 = 0.03; // 500 micron diameter
  //ex_expo  = 2.9;      //3.2, not used

  Computed = NULL;
  z_t = NULL , y_t = NULL;

  localmax = 0;                // max numbver of generations in a tree
  Z_out = 0;
  z_out = 1;
  y_out = 2;
  writeout = 0; // Flag to writeout Z, Y ,z, or y

  Computed = new cvOneDComplex*[Maxgen];
  for (int i = 0; i < Maxgen; i++) {
    Computed[i] = new cvOneDComplex[Maxgen];
  }

  Complex_ZERO = cvOneDComplex(0,0);
  Complex_ONE = cvOneDComplex(1,0);
  Complex_I = cvOneDComplex(0,1);

  lengthN = 0;

}

cvOneDMetteImpedance::~cvOneDMetteImpedance(){
  for (int i = 0; i < Maxgen; i++) {
    delete Computed[i];
  }
  delete [] Computed;
}

double* cvOneDMetteImpedance::calculateImpedanceTime(double rootR, double lrr,double Period, 
                                                     int tmstps){
  bool returnImp = true;
  bool returnFourier = false;
  return calculateImpAdmit( rootR, lrr, asym, Period, tmstps, returnImp, returnFourier);

}

double* cvOneDMetteImpedance::calculateImpedanceFourier(double rootR, double lrr,double Period, 
                                                        int tmstps){
  bool returnImp = true;
  bool returnFourier = true;
  return calculateImpAdmit( rootR, lrr, asym, Period, tmstps, returnImp, returnFourier);

}

double* cvOneDMetteImpedance::calculateRCRImpedanceTime(double Rd, double Ru, 
                                                        double Cap, double Period,int tmstps){
  bool returnImp = true;
  bool returnFourier = false;
  RCR=true;
  R1=Rd;
  R2=Ru;
  C1=Cap;
  return calculateImpAdmit( 0.0, 10.0, asym, Period, tmstps, returnImp, returnFourier );

}
  
double* cvOneDMetteImpedance::calculateRCRImpedanceFourier(double Rd, double Ru, 
                                                           double Cap, double Period,int tmstps){
  bool returnImp = true;
  bool returnFourier = true;
  RCR=true;
  R1=Rd;
  R2=Ru;
  C1=Cap;
  return calculateImpAdmit( 0.0, 10.0, asym, Period, tmstps, returnImp, returnFourier );
}
  
double*  cvOneDMetteImpedance::calculateAdmittanceTime(double rootR, double lrr,double Period, 
                                                       int tmstps){
  bool returnImp = false;
  bool returnFourier = false;
  return calculateImpAdmit( rootR, lrr, asym, Period, tmstps, returnImp, returnFourier);
}
  
double*  cvOneDMetteImpedance::calculateAdmittanceFourier(double rootR, double lrr,double Period, 
                                                          int tmstps){
  bool returnImp = false;
  bool returnFourier = true;
  return calculateImpAdmit( rootR, lrr, asym, Period, tmstps, returnImp, returnFourier);
}

double* cvOneDMetteImpedance::calculateRCRAdmittanceTime(double Rd, double Ru, 
                                                         double Cap, double Period,int tmstps){
  bool returnImp = false;
  bool returnFourier = false;
  RCR=true;
  R1=Rd;
  R2=Ru;
  C1=Cap;
  return calculateImpAdmit( 0.0, 10.0, asym, Period, tmstps, returnImp, returnFourier );
}

double* cvOneDMetteImpedance::calculateRCRAdmittanceFourier(double Rd, double Ru, 
                                                            double Cap, double Period,int tmstps){
  bool returnImp = false;
  bool returnFourier = false;
  RCR=true;
  R1=Rd;
  R2=Ru;
  C1=Cap;
  return calculateImpAdmit( 0.0, 10.0, asym, Period, tmstps, returnImp, returnFourier );
}

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

double* cvOneDMetteImpedance::calculateImpAdmit(double rootR, double ARGlrr, double ARGasym,
                                                double Period, int tmstps, bool returnImp, 
                                                bool returnFourier){

  //int nb_terms=0; // Number of terms to solve for in tree
  double df;
  double trm_rst = 0.00; // Frequency stepsize and terminal resistance
  double * Freq = new double[tmstps+1];
  double * Omega = new double[tmstps+1];
  cvOneDComplex* Z_om; // Root impedance at all frequencies
  cvOneDComplex* Z_hat;
  cvOneDComplex* Y_hat;   
  cvOneDComplex* z_xt;
  cvOneDComplex* y_xt; // Impedance and conductance in time domain
  int N = tmstps;  // Number of timesteps must equal number of frequencies.
  lengthN = N;

  // automatically switch to structured trere of lrr < 1
  if(lrr < 1.0){
    METTE = true;
  }

 // Set mette's default values  (shouldn't these appear in defaults above??)
 if(METTE){
   expo=2.7;
   asym = 0.4048;
   r_min=lrr;
   lrr=50.0;
   fprintf(stdout,"Use Mette's method. rmin=%lf  lrr=%lf \n",r_min,lrr);
 }
      
  // This is from mette's dimension-less parameters.
  if(F95||METTE){
    Period = Period*q/Lr/Lr; // *10
  }

  // Physical parameters
  df  = 1.0/Period; // Frequency interval.

  for(int i=-N/2;i<=N/2;i++){
    Freq[i+N/2] = df*(double)i;  
    Omega[i+N/2]= 2*M_PI*Freq[i+N/2];//    ! Freq.-vector scaled by 2pi.     
  }
  
  // "globals"
  if(!RCR){
    //asym = this.asym;
    beta = std::pow(std::pow(asym,expo/2.0)+1.0,-1.0/expo); // Scaling parameter.
    alpha = std::sqrt(asym)*beta;
    r_root = rootR;  //
    lrr = ARGlrr;  //
    localmax = 0; // initialize
    //nb_terms = 0;
  }else{
    asym = ARGasym;
  }

  // Compute the impedance at the root of the structured tree.
  Z_om = comp_imp(tmstps,Omega);
  Z_hat = Z_om;
  //fprintf(stdout,"Z(0) = %lf\n",Z_om[0].Real());
  //Y_hat = 1/Z_om; find conductance
  Y_hat = new cvOneDComplex[lengthN];
  for(int i=0;i<lengthN;i++){
    Y_hat[i] = Complex_ONE/(Z_hat[i]);
  }

  y_xt = InverseFourier(Y_hat); // Inverse Fourier Transform

/*   //only low frequencies-irene's trials
  if(lengthN/2>10){//line to change for cut frequency
    for(int i=10;i<lengthN/2;i++){//line to change for cut frequency
      Z_hat[i] = Complex_ONE*0.0;
    }
// Apply self-adjointness and rotate
    for(int i=0;i<N/2;i++){
      Z_hat[N-1-i]= conj(Z_hat[i+1]);
    }
    Z_hat[N/2]=Complex_ONE*0.0;
  }
*/
  
  z_xt = InverseFourier(Z_hat); // Fourier transform
  
  y_t=new double[lengthN];// in this section, pick out real part  
  z_t=new double[lengthN];
  double* Y_omR = new double[lengthN];
  double* Z_omR = new double[lengthN];
   
  for(int i=0;i<lengthN;i++){
    y_t[i]=y_xt[i].Real() ;//return this for admittance y(t)
    z_t[i]=z_xt[i].Real() ;//return this for impedance z(t)
    //Y_omR[i]=Y_hat[i].Real();//return this for admittance Y(w)
    Y_omR[i]=sqrt(Y_hat[i].Real()*Y_hat[i].Real()+Y_hat[i].Imag()*Y_hat[i].Imag());//return the modulus for admittance Y(w) - vie 112904
    Z_omR[i]=sqrt(Z_om[i].Real()*Z_om[i].Real()+Z_om[i].Imag()*Z_om[i].Imag());//return the modulus for impedance Z(w) - vie 112904
    //Y_omR[i]=Y_hat[i];//return the whole complex number
    //Z_omR[i]=Z_om[i];//return the whole complex number
    //Z_omR[i]=Z_om[i].Real();//return the real part
    //Z_omR[i]=Z_om[i].Imag();//return the imaginary part

    if(F95){// use non-dimensional terms
      y_t[i]*= 1.0*q/g/rho;// this gives same y as F90
      // y_t[i]*= N/Period;// need to run with Mette's analysis code
    }
  }

  delete [] Freq;
  delete [] Omega;
  delete [] Y_hat; 
  delete [] Z_om;
  delete [] y_xt;
  delete [] z_xt;
  if (returnImp == true) {
  delete [] y_t;
  delete [] Y_omR;
  if(returnFourier == true){
    delete [] z_t;
    return Z_omR;
  }
  else{
    delete [] Z_omR;
    return z_t;
  }
}
else{
  delete [] z_t;
  delete [] Z_omR;
  if(returnFourier == true){
    delete [] y_t;
    return Y_omR;
  }
  else{
    delete [] Y_omR;
    return y_t;
  }
  }

}  // End subroutine CalculateImpAdmit

// =================
// COMPUTE IMPEDANCE
// =================
cvOneDComplex* cvOneDMetteImpedance::comp_imp(int N,double* Omega){

  //printf("COMP_IMP\n");
  //if(cvOneDGlobal::debugMode){
  //  for(int loopA=0;loopA<N+1;loopA++){
  //    printf("Input Omega, comp n. %d: %f\n",loopA,Omega[loopA]);
  //  }
  //}
  //getchar();

  cvOneDComplex  temp;
  cvOneDComplex* Z_om = new cvOneDComplex[N];

  // Omega contains tmstps+1 frequencies because when computing the
  // impedance it is easier to invert it when it is computed for
  // all positive frequencies and since the interval goes from
  // [-N/2/Tper:N/2/Tper-df],
  // we include the frequency N/2/Tper in Omega and from this we get
  // Z(-N/2/Tper) we then end up throwing pushing Z(N/2/Tper) out.

  for(int i=0; i < N;i++){
    Z_om[i] = Complex_ZERO; // Initialize Z_om with zeros
  }

  // For all the positive frequencies compute the impedance.
  // Since we know that the system is self-adjoint we don't
  // have to compute the negative frequencies as well.
  for(int k = N/2; k <= N ;k++){     // this works!!!
    for(int i=0;i<Maxgen;i++){       // For each frequency make sure we don't carry
      for(int j=0;j<Maxgen;j++){     // any values with us from the computations for the                                      
        Computed[i][j]=Complex_ZERO; // previous frequency, so make
      }                              // sure Computed is zero.                      
    }
    // Since Z_om only has N places leave result for the frequency k at
    // Z_om(k-1) we will later make up for this.
    if (RCR){
      Z_om[k-N/2] = RCRImpedance(R2,R1,C1,Omega[k]);
    }else{
      Z_om[k-N/2] = Z0func(Omega[k],0,0);
    }
  }

  // Apply self-adjointness and rotate
  for(int i=0;i<N/2;i++){
    Z_om[N-1-i]= conj(Z_om[i+1]);
  }
  Z_om[N/2]=Complex_ONE*(Z_om[N/2].Real());
   
  return Z_om;
} //end function comp_imp

// ========================================================
// Compute charachteristic ipedance using geometry and root
// impedance of downstream vessels.
// Recursive Function
// ========================================================
cvOneDComplex cvOneDMetteImpedance::Z0func(double omega_k, int alpha_pow,int beta_pow){

  double lrr_d = lrr; // if METTE lrr==50
  double r_d  = std::pow(alpha,alpha_pow)*std::pow(beta,beta_pow)*r_root; // Current radius
  int generations = alpha_pow + beta_pow;  // Current generation.

  if(!METTE){
    expo = 2.51;
    asym = 0.4048;
    if( r_d <= radius_flag1){
      expo = expo1;
      asym = asym1;
      double tmpbeta = std::pow(std::pow(asym,expo/2.0)+1.0,-1.0/expo); // Scaling parameter.
      double tmpalpha = std::sqrt(asym)*tmpbeta;
      r_d  = std::pow(tmpalpha,alpha_pow)*std::pow(tmpbeta,beta_pow)*r_root; // Current radius
      if(r_d <= radius_flag2){
        expo=expo2;
        asym=asym2;
        tmpbeta = std::pow(std::pow(asym,expo/2.0)+1.0,-1.0/expo);// Scaling parameter.
        tmpalpha = std::sqrt(asym)*tmpbeta;
        r_d  = std::pow(tmpalpha,alpha_pow)*std::pow(tmpbeta,beta_pow)*r_root; // Current radius
        lrr_d = lrr;
      }
    }
  }

  double r = r_d/Lr;  // Radius at root, dimension-less. shoud fall in F95 category, but Lr == 1.
  double l = lrr_d*r; // Length of vessel segment.
  if(EXERCISE && r <= radius_flag0 ){
    r *= ex_fact;
  }
  double A  = M_PI*r*r;                    // Cross-sectional area, cm^2.
  double nu = mu/rho;                           // Kinematic blood viscosity, cm^2/s.
  double D  = 3.0*A/2.0/(k1*std::exp(k2*r)+k3); // Compliance.  cm^4/dyne
  if(F95){
    D = D*rho*g/Lr;// Non dimensionalized -? Mette's rf
  }

  // Womersleys parameter
  double wom  = r*std::sqrt(omega_k/nu);

  cvOneDComplex g_omega, // wave speed times compliance,  cm^4/dynes * cm/s = cm^5/(dynes*s)
                c_omega, // wave speed, cm/s
                  kappa, // unitless
        ZL,Zl_0,Zr_0,Z0; // impedance, dynes * s / cm^5

  cvOneDComplex f10 = getF10(wom); // sqrt(1-Fj)
  // printf("f10 Imag Value, wom: %e %e\n",f10.Imag(),wom);
  // Define Wave Speed
  c_omega = Complex_ONE*sqrt(A/D/rho)*(f10); // Wave speed, cm/s.
  if(F95){
    c_omega = Complex_ONE*sqrt(A/D/Fr2)*(f10); // Wave speed, cm/s Fr2=sqrt(q)/g/Lr5.
  }
  g_omega = c_omega*(D); // cm^5/dyne/s.
  // printf("g_omega: %e\n",g_omega.Real());
  // printf("Parames A D rho: %e %e %e\n",A,D,rho);
  //  ! Temporary function of omega_k. Unitless
  if (omega_k != 0)  {
    kappa  = Complex_ONE*(omega_k*l)/(c_omega);
    if(F95) kappa = Complex_ONE*(omega_k* Lr*Lr*Lr/q *(l/Lr))/( c_omega);
  } else {
    kappa  = Complex_ZERO;
  }

  // Determine the resistance of the root of the terminal branches. Recursive call
  // if (generations >= Maxgen)
  // printf("r, rmin: %e %e\n",r,r_min);
  if (r <= r_min || (EXERCISE && r <= r_min * ex_fact) ){
     // cannot use this scheme for exercise where we want to dilate the arterioles and capilaries.
     if (generations > localmax){
       localmax = generations;
     }
     // printf("ZL Set to Zero!\n");
     ZL = Complex_ZERO; // Set default terminal resistance
   } else {
     // Get Z0 recursively at reduced radii.
     // Left
     if (cmodulus(Computed[alpha_pow+1][beta_pow]) != 0.0)
        Zl_0 = Computed[alpha_pow+1][beta_pow];
     else 
        Zl_0 = Z0func(omega_k,alpha_pow+1,beta_pow);
     // Right
     if (cmodulus(Computed[alpha_pow][ beta_pow+1]) != 0.0)
       Zr_0 = Computed[alpha_pow][beta_pow+1];
     else 
       Zr_0 = Z0func(omega_k,alpha_pow,beta_pow+1);

     // Prediction of the resulting impedance from the recursion formula.
     ZL = Complex_ONE/((Complex_ONE/(Zl_0)+(Complex_ONE/(Zr_0))));
   }

   // Finally get Z0(omega) as theoretically derived.
   if (cmodulus(g_omega) != 0.0){
     cvOneDComplex t1  = Complex_I*( sin(kappa) )/(g_omega)+( cos(kappa)*(ZL) );
     cvOneDComplex t2  = cos(kappa)+( Complex_I*(g_omega)*( sin(kappa)*(ZL) ) );
     Z0  = t1/(t2);// Input Impedance, dynes s /cm^5
   } else {
     // Not Zero Area and Radius
     Z0  = Complex_ONE*(8.0*mu*l/(A*r*r))+(ZL);
     if(F95) {
       Z0  = Complex_ONE*(8.0*nu*Lr*l*Fr2/(A*r*r)/q)+(ZL);
     }
   }
   // printf("Z0func Return: %f\n",Z0.Real());
   return Computed[alpha_pow][beta_pow] = Z0;
}

cvOneDComplex cvOneDMetteImpedance::RCRImpedance(double R, double Rd, double C, double omega){
  cvOneDComplex numerator, denominator, Z;
  numerator = cvOneDComplex(R + Rd,omega *R*Rd*C);  // Deleted "new" (nate)
  denominator = cvOneDComplex(1, omega*Rd*C);  // Deleted "new" (nate)
  Z = numerator/(denominator); // Computes impedance
  return Z;
}

cvOneDComplex cvOneDMetteImpedance::getF10(double wom){
  if(fabs(wom)<1.0e-10) {
    return Complex_ZERO;
  }else if(wom > 3.0){ // Used by mette
    return  sqrt(Complex_ONE-(Complex_ONE*(2.0)/(sqrt(Complex_I)*(wom))*(1.0+1.0/2.0/wom)));
  }else if(wom > 2.0){ // Used by mette
    return Complex_ONE*(3.0-wom)*(sqrt(Complex_I*(wom*wom/8.0)+(Complex_ONE*(wom*wom/48.0*wom*wom))))+(Complex_ONE*(wom-2.0)*(sqrt(Complex_ONE-(Complex_ONE*(2.0)/(sqrt(Complex_I)*(wom))*(1.0+1.0/2.0/wom)))));
  }else{
    return sqrt(Complex_I*(wom*wom/8.0)+(Complex_ONE*(wom*wom/48.0*wom*wom))); // Mette's code
  }
}

// Take inverse Fourier transform of array
cvOneDComplex* cvOneDMetteImpedance::InverseFourier(cvOneDComplex* Y){
  // do transform by hand, does not require to be power of 2, but should be even
  cvOneDComplex* y = new cvOneDComplex[lengthN];
  int n=lengthN;
  for(int i=0;i<lengthN;i++){
    cvOneDComplex sum=Complex_ZERO;
    for(int j=0;j<n;j++){
      sum = sum + (Y[j]*(exp(Complex_I*(2.0*M_PI*(double)i*(double)j/(double)n))));
    }
    // y[i]=sum/((double)n);
    y[i]=sum; // modified by vie for a more consistent definition of impedance/admittance in time p(t)=1/T *integral(z(t-tau)q(tau),tau,t-T,t) 
  }
  return y;
}

// Fourier isn't really used, but figured I'd leave it just in case (nate)
cvOneDComplex* cvOneDMetteImpedance::Fourier(cvOneDComplex* y){
  // do transform by hand, does not require to be power of 2, but should be even
  cvOneDComplex* Y=new cvOneDComplex[lengthN];
  for(int i=0;i<lengthN;i++){
    cvOneDComplex sum=Complex_ZERO;
    int n=lengthN;
    for(int j=0;j<n;j++){
      sum=sum+(y[j]*(exp(Complex_I*(-2.0*M_PI*(double)i*(double)j/(double)n))));
    }
    Y[i]=sum/((double)n);
  }
  return Y;
}

// kimhj 08-18-2004 added
// iterate until the root DC impedance matches the resistance within 
// the epsilon. And the maximum iteration is MAX_iter. 
// But it converges really fast(The relationship is linear). 
double cvOneDMetteImpedance::resistBasedTree(double rootR, double resist, double Period, int tmstps){
  
  double * Omega = new double[tmstps+1];
  cvOneDComplex* Z_om;  
  double epsilon=1.0E-7;
  int MAX_iter=20;
  double temp1, temp2,temp;
  double x1,x2,xopt;
  double comp1, comp2,comp;
  int iter=1;

  double df  = 1/Period; //! Frequency interval. 
  for(int i=-tmstps/2;i<=tmstps/2;i++){
    Omega[i+tmstps/2]= 2*M_PI*df*(double)i;    
  }
  
  beta = std::pow(std::pow(asym,expo/2.0)+1.0,-1.0/expo);//    ! Scaling parameter.
  alpha = std::sqrt(asym)*beta;
  r_root = rootR; 
  lrr=calculateLRR(r_root,resist);
  x1=lrr;
  Z_om = comp_imp(tmstps,Omega);
  temp1=Z_om[0].Real();
  comp1=std::abs(temp1-resist)/resist;
  x2=x1+lrr;
  temp=temp1;
  comp=comp1;
  
  while(comp>=epsilon && iter <=MAX_iter){
    lrr=x2;
    Z_om = comp_imp(tmstps,Omega);
    temp2=Z_om[0].Real();
    delete [] Z_om;
    comp2=std::abs(temp2-resist)/resist;
    xopt=(resist-temp1)*(x2-x1)/(temp2-temp1)+x1;
    lrr=xopt;
    Z_om = comp_imp(tmstps,Omega);
    temp=Z_om[0].Real();
    delete [] Z_om;
    comp=std::abs(temp-resist)/resist;
    x2=xopt;
    iter=iter+1;
  } // end while loop
    
  lrr=xopt;
  fprintf(stdout,"Final lrr=%lf\n",lrr);
  delete [] Omega;      
  return lrr;
}
    
// kimhj 08-16-2004 add a function which calculates the length-to-radius ratio
// according to the given resistance value. 
// It receives rootR and resistance as input variables and gives the lrr
// which estimates the DC root impedance adjacent to the resistance. 
// It is approximated using that asym=asym1 and expo=expo1
// and doesn't change as the radius changes
 
double cvOneDMetteImpedance::calculateLRR(double rootR, double resist){
  double testasym=asym1;
  double testexpo=expo1;
  double fraction=r_min/rootR;
  double testbeta = std::pow(std::pow(testasym,testexpo/2.0)+1.0,-1.0/testexpo);
  double testalpha = std::sqrt(testasym)*testbeta;
  int gens=0.5*(floor(log(fraction)/log(testalpha))+floor(log(fraction)/log(testbeta)));
  double part=std::pow(testalpha,3)+std::pow(testbeta,3);
  lrr=resist*M_PI*std::pow(rootR,3)*(part-1)/8/mu/(part-1/std::pow(part,gens));
  return lrr;
}
