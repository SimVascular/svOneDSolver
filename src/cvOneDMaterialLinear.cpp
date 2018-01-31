
//  MaterialLinear.cxx - Source for a class to maintain material properties
//  ~~~~~~~~~~~~
//  
//  SYNOPSIS...This class maintains MaterialLinear Properties of the subdomain.
//  08/02/04 VIE : added the three functions that are needed for conservative formulation
//  08/02/04 VIE : changed other functions to be consistent with linear model

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "assert.h"
#include "cvOneDGlobal.h"
#include "cvOneDMaterialLinear.h"
#include "cvOneDUtility.h"

using namespace std;

cvOneDMaterialLinear::cvOneDMaterialLinear(){   
}

cvOneDMaterialLinear::~cvOneDMaterialLinear(){
}

//
//  The following three methods are required by the conservation form
//

//Integral of S from ref P to P(t)
double cvOneDMaterialLinear::GetIntegralpS(double area, double z)const{ 
  double EHR = GetEHR(z);//for this model it is a constant
  double So_ = GetS1(z);
  double IntegralpS = EHR/3.0*So_*(area/So_*sqrt(area/So_)-1.0);
  return IntegralpS;
}

//Integral of dS(p,z,t)dz from ref P to P(t)
double cvOneDMaterialLinear::GetIntegralpD2S (double area, double z)const{
  double EHR = GetEHR(z);//for this model it is a constant
  double So_ = GetS1(z);
  double dSo_dz = GetDS1Dz(z);
  double IntegralpD2S = EHR/3.0*dSo_dz*(area/So_*sqrt(area/So_)-1.0);
  return IntegralpD2S;
}

//required for wave BC only 
double cvOneDMaterialLinear::GetRefWaveSpeed(double area)const{ 
  fprintf(stderr,"ERROR:  GetRefWaveSpeed not implemented!!\n\n");
  assert(0);
  return 0;
  // double EHR   = GetEHR(z);//for this model it is a constant
  // double co = sqrt(EHR/2.0/density);
  // return co;
}

cvOneDMaterialLinear::cvOneDMaterialLinear(const cvOneDMaterialLinear &rhs){
  cvOneDMaterial::operator=(rhs);
  ehr = rhs.ehr;
}
 
cvOneDMaterialLinear& cvOneDMaterialLinear::operator=(const cvOneDMaterialLinear &that){
  if (this != &that) {
    cvOneDMaterial::operator=(that);
    ehr = that.ehr;
  }
  return *this;
}

void cvOneDMaterialLinear::SetEHR(double ehr_val){  
  ehr = ehr_val;
}
  
double cvOneDMaterialLinear::GetProperty(char* what)const{
  if( strcmp( what, "density") == 0){
    return density;
  }else if( strcmp( what, "dynamic viscosity") == 0){
    return dynamicViscosity;
  }else if( strcmp( what, "kinematic viscosity") == 0){
    return kinematicViscosity;
  }else if( strcmp( what, "constant n") == 0){
    return profile_exponent;
  }else if( strcmp( what, "delta") == 0){
    return delta;
  }else if( strcmp( what, "N") == 0){
    return N;
  }else{
    abort();
  }
  // Return
  return 0.0;
}

double cvOneDMaterialLinear::GetEHR(double z)const{
  return ehr;
}

void cvOneDMaterialLinear::SetAreas_and_length(double S_top,double S_bottom,double z){
  Stop = S_top;//inlet
  Sbot = S_bottom;//outlet
  len = z;
  //  cout<<"length "<< len<<" areas "<< Stop<<" "<<Sbot<<endl;
}

double cvOneDMaterialLinear::GetS1(double z)const{
  double r    = Getr1(z);
  double area = r*r*M_PI;
  return area;
}

double cvOneDMaterialLinear::Getr1(double z)const{
  // Linearly interpolated r
  double r_top = sqrt(Stop/M_PI);
  double r_bot = sqrt(Sbot/M_PI);
  double r     = ((z-len)/(-len))*(r_top - r_bot) + r_bot;  
  // mette's method
  //r=r_top*exp(z/len*log(r_bot/r_top)); //RLS: Use the linearly interpolated version.
  return r;  
}

//not used
double cvOneDMaterialLinear::GetDS1Dz(double z)const{
  // Since vessel geometry will always be linear in
  // space, this is just a constant  
  /*
  double der;
  der = ((Sbot - Stop)/(len));
  return der;
  */
  // this shouldn't change solution, but try using dr/dz and compute S  
  double drodz  = GetDr1Dz(z) ;
  double dsodro = 2.0*M_PI*Getr1(z);
  return dsodro*drodz;  // slightly increased pressure/decreased area
}

double cvOneDMaterialLinear::GetDr1Dz(double z) const{
  double r_top = sqrt(Stop/M_PI);
  double r_bot = sqrt(Sbot/M_PI);
  double drodz = ((r_bot - r_top)/len) ; //RLS: These values are the initial radii.

  // mette's method
  // drodz=1/len*log(r_bot/r_top)*Getr1(z); //RLS: Use the linearly interpolated version.
    
  return drodz;
}

double cvOneDMaterialLinear::GetWaveSpeed(double pressure, double z)const{
  // KLUDGE - unimplimented (needed?)
  return 0.0;
}

double cvOneDMaterialLinear::GetArea(double pressure, double z)const{
  // NOTE: o "So_" is the LSA under pressure p1_.
  //         This property comes from the subdomain    
  //
  //       o Po is the the zero transmural pressure 
  double pres = pressure;  
  double So_  = GetS1(z); 
  double EHR  = GetEHR(z);//*4/3
  // double area1 = (16.0*So_*pow(EHR,2))/pow(-3.0*pres + 3.0*p1_ + 4.0*EHR,2);//I moved 4/3 into EHR term, so this will be off now, but did work.
  // double area = So_/pow(1-(pres-p1_)/EHR,2);// this makes more sense to me, bns 8/22/02. gives same answer
  
  // RLS: This is the area computation using the "pressure-strain" modulus, EHR.
  double area = So_*pow(1.0+(pres-p1_)/EHR,2.0);

  if(cvOneDGlobal::debugMode){
    printf("So_: %e\n",So_);
    printf("pres: %e\n",pres);
    printf("p1_: %e\n",p1_);
    printf("EHR: %e\n",EHR);
    printf("Area: %e\n",area);
    fflush(stdout);
  }

  return area;
}

double cvOneDMaterialLinear::GetPressure(double S, double z)const{
  // Again we need to get So_ from the subdomain. 
  // Then we impliment Olufsen's constitutive law...
  double So_   = GetS1(z);
  double EHR   = GetEHR(z);  // From Olufsen's paper
  // double press = (p1_ + EHR*(1.-sqrt(So_/S)));//dynes/cm^2
  double press = p1_ + EHR*(sqrt(S/So_)-1.0);//VIE for linear model dynes/cm^2
  return press;
}


double cvOneDMaterialLinear::GetDpDS(double S, double z)const{
  double EHR = GetEHR(z);
  double So_ = GetS1(z);
  double ro  = Getr1(z);
  // double dpds=0.5* EHR * sqrt(So_/S/S/S) ;
  double dpds=0.5* EHR/sqrt(So_*S) ;//VIE for linear model 
  return dpds;
}

// not used
double cvOneDMaterialLinear::GetD2pDS2( double area, double z) const{
  double EHR = GetEHR(z);
  double So_ = GetS1(z);
  return - EHR /4.0 /sqrt(So_)/sqrt(pow(area, 3));//VIE for linear model  
}

double cvOneDMaterialLinear::GetOutflowFunction(double pressure, double z)const{
  return 0.; // This is not used in our model
}

double cvOneDMaterialLinear::GetDOutflowDp(double pressure, double z)const{
  return 0.; // Nor is this.
}

// Careful!! this D2p(S,z)Dz first derivative, 2nd variable IV 08-02-04 
double cvOneDMaterialLinear::GetDpDz(double S, double z)const{
  double So_   = GetS1(z);
  double EHR   = GetEHR(z);  
  double r     = sqrt(S/M_PI);
  double ro    = Getr1(z);
  double drodz = GetDr1Dz(z);//if straight tube ->0.0 //RLS: Returns dro/dz
  // if(NON_DIM) dEHRdro*=1/density/981.;

  // double dpdz = drodz*(dEHRdro*(1.0-(ro/r))-EHR/r); //RLS: I'm not using this constitutive model.

  // RLS: The following calculation is actually dp_hat/dz as written in Jing et al's equation 18.
  double dpdz = drodz*(-EHR*r/ro/ro); //check by VIE 08-02-04
    
  // mette's causes pressure to increase in longer tubes
  // double theta= M_PI/2.*(S/So_ - 1.);
  // if(NON_DIM) dEHRdr*=1/density/981.;
  // double p=EHR/M_PI;
  // dpdz=drdz*(dEHRdr*tan(theta)- p*S/(ro*ro*ro*cos(theta)*cos(theta)));
  return dpdz;
}

