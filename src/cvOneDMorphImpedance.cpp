# include "cvOneDMorphImpedance.h"

float* pm_ImpComp (float radius, int min_order, float freq, int num_freq){
  fprintf(stderr,"ERROR:  Vmod code no longer exists!\n");
  fflush(stderr);
  exit(-1);
  return (float*)NULL;
}

cvOneDMorphImpedance::cvOneDMorphImpedance() {

  // variables
  r_root = 0;           // root (starting) radius of tree
  z_t = NULL , y_t = NULL;

  Complex_ZERO=cvOneDComplex(0,0);
  Complex_ONE=cvOneDComplex(1,0);
  Complex_I=cvOneDComplex(0,1);

  lengthN = 0;

}

cvOneDMorphImpedance::~cvOneDMorphImpedance(){
}

double* cvOneDMorphImpedance::calculateImpedanceTime(double rootR, int min_order, 
                                                     double Period, int tmstps){
  bool returnImp = true;
  bool returnFourier = false;
  return calculateImpAdmit( rootR, min_order, Period, tmstps, returnImp, returnFourier);
}
  
double* cvOneDMorphImpedance::calculateImpedanceFourier(double rootR, int min_order,
                                                        double Period, int tmstps){
  bool returnImp = true;
  bool returnFourier = true;
  return calculateImpAdmit( rootR, min_order, Period, tmstps, returnImp, returnFourier);
}

double* cvOneDMorphImpedance::calculateAdmittanceTime(double rootR, int min_order, 
                                                      double Period, int tmstps){
  bool returnImp = false;
  bool returnFourier = false;
  return calculateImpAdmit( rootR, min_order, Period, tmstps, returnImp, returnFourier);
}

double* cvOneDMorphImpedance::calculateAdmittanceFourier(double rootR, int min_order, 
                                                         double Period, int tmstps){
  bool returnImp = false;
  bool returnFourier = true;
  return calculateImpAdmit( rootR, min_order, Period, tmstps, returnImp, returnFourier);
}

double* cvOneDMorphImpedance::calculateImpAdmit(double rootR, int min_order, 
                                                double Period, int tmstps, 
                                                bool returnImp, bool returnFourier){

  int nb_terms = 0; // number of terms to solve for in tree
  double  df, trm_rst = 0.00; // frequency stepsize and terminal resistance
  double* Freq = new double[tmstps+1];
  double* Omega = new double[tmstps+1];
  cvOneDComplex* Z_om; // Root impedance at all frequencies
  cvOneDComplex* Z_hat;
  cvOneDComplex* Y_hat;   
  cvOneDComplex* z_xt;
  cvOneDComplex* y_xt; // impedance and conductance in time domain
  int N = tmstps; // number of timesteps must equal number of frequencies.
  lengthN = N;
  int i;

  // Physical parameters
  df  = 1/Period; // Frequency interval
 
  for(i=-N/2;i<=N/2;i++){
    Freq[i+N/2] = df*(double)i;  
    Omega[i+N/2]= 2*M_PI*Freq[i+N/2];// Freq.-vector scaled by 2pi.     
  }

  double dimp[1000];
  float *imp;
  int num_freq = tmstps/2+1;
  imp = pm_ImpComp ( (float)rootR,  min_order, (float)(2.0*M_PI/Period), num_freq);
  
  for (i = 0; i < num_freq; i++) {
    dimp[2*i] = imp[2*i] / (M_PI*rootR*rootR);
    dimp[2*i+1] = imp[2*i+1] / (M_PI*rootR*rootR);
  }

  // Compute the impedance at the root of the structured tree.
  Z_om = new cvOneDComplex[N];

  for(int i=0; i < N;i++) {
    Z_om[i] = Complex_ZERO; // Initialize Z_om with zeros
  }

  for (int i=0;i<=N/2;i++) {
    Z_om[i]=Complex_ONE*dimp[2*i]+Complex_I*dimp[2*i+1];
  }

  // Apply self-adjointness and rotate
  for(int i=0;i<N/2;i++){
    Z_om[N-1-i]= conj(Z_om[i+1]);
  }
  Z_om[N/2]=Complex_ONE*(Z_om[N/2].Real());

  Z_hat=Z_om;
  // fprintf(stdout,"Z(0) = %lf\n",Z_om[0].Real());
  // Y_hat = 1/Z_om; find conductance
  Y_hat = new cvOneDComplex[lengthN];
  for(int i=0;i<lengthN;i++){
    Y_hat[i] = Complex_ONE/(Z_hat[i]);
  }
  
  y_xt = InverseFourier(Y_hat);// Fourier transform 
  z_xt = InverseFourier(Z_hat);// 
  
  y_t = new double[lengthN];// in this section, pick out real part  
  z_t = new double[lengthN];
  double* Y_omR = new double[lengthN];
  double* Z_omR = new double[lengthN];
   
  for(int i=0;i<lengthN;i++){
    y_t[i]=y_xt[i].Real() ;//return this for admittance y(t)
    z_t[i]=z_xt[i].Real() ;//return this for impedance z(t)
    // Y_omR[i]=Y_hat[i].Real();//return this for admittance Y(w)
    Y_omR[i]=sqrt(Y_hat[i].Real()*Y_hat[i].Real()+Y_hat[i].Imag()*Y_hat[i].Imag());//return the modulus for admittance Y(w) - vie 112904
    // Z_omR[i]=Z_om[i].Real();//return this for impedance Z(w)
    Z_omR[i]=sqrt(Z_om[i].Real()*Z_om[i].Real()+Z_om[i].Imag()*Z_om[i].Imag());//return the modulus for impedance Z(w) - vie 112904
  }

  delete [] Freq;
  delete [] Omega;
  delete [] Y_hat; 
  delete [] Z_om;
  delete [] y_xt;
  delete [] z_xt;

  if(returnImp == true){
    delete [] y_t;
    delete [] Y_omR;
    if(returnFourier == true){
      delete [] z_t;
      return Z_omR;
    }else{
      delete [] Z_omR;
      return z_t;
    }
  }else{
    delete [] z_t;
    delete [] Z_omR;
    if(returnFourier == true){
      delete [] y_t;
      return Y_omR;
    }else{
      delete [] Y_omR;
      return y_t;
    }
  }
}  // end subroutine CalculateImpAdmit

// Take inverse Fourier transform of array
cvOneDComplex* cvOneDMorphImpedance::InverseFourier(cvOneDComplex* Y){
  // do transform by hand, does not require to be power of 2, but should be even
  cvOneDComplex* y = new cvOneDComplex[lengthN];
  int n=lengthN;
  for(int i=0;i<lengthN;i++){
    cvOneDComplex sum=Complex_ZERO;
    for(int j=0;j<n;j++){
      sum=sum+(Y[j]*(exp(Complex_I*(2.0*M_PI*(double)i*(double)j/(double)n))));
    }
    // y[i]=sum/((double)n);
    y[i]=sum;//modified by vie for a more consistent definition of impedance/admittance in time p(t)=1/T *integral(z(t-tau)q(tau),tau,t-T,t) 
  }
  return y;
}

// Fourier isn't really used, but figured I'd leave it just in case (nate)
cvOneDComplex* cvOneDMorphImpedance::Fourier(cvOneDComplex* y){
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
