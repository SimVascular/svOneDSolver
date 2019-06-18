
//  MaterialLinear.cxx - Source for a class to maintain material properties
//  ~~~~~~~~~~~~
//
//  SYNOPSIS...This class maintains MaterialLinear Properties of the subdomain.

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
  double EHR = GetEHR(z); //for this model it is a constant
  double So_ = GetS1(z);
  double IntegralpS = EHR/3.0*So_*(area/So_*sqrt(area/So_)-1.0);
  return IntegralpS;
}

//Integral of dS(p,z,t)dz from ref P to P(t)
double cvOneDMaterialLinear::GetIntegralpD2S (double area, double z)const{
  double EHR = GetEHR(z); //for this model it is a constant
  double So_ = GetS1(z);
  double dSo_dz = GetDS1Dz(z);
  double IntegralpD2S = EHR/3.0*dSo_dz*(area/So_*sqrt(area/So_)-1.0);
  return IntegralpD2S;
}

cvOneDMaterialLinear& cvOneDMaterialLinear::operator=(const cvOneDMaterialLinear &that){
  if (this != &that) {
    cvOneDMaterial::operator=(that);
    ehr = that.ehr;
    p1_=that.PP1_; //impose P refrence here otherwise p1_ is not the value set in the input file.
    printf("this that set ehr =%f p1_=%f \n",ehr, p1_);
  }
  return *this;
}

void cvOneDMaterialLinear::SetEHR(double ehr_val, double pref_val){
  ehr = ehr_val;
  PP1_= pref_val;  // add additional commend to set P refrence otherwise p1_ is not the value set in the input file.
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
  return 0.0;
}

double cvOneDMaterialLinear::GetEHR(double z)const{
  return ehr;
}

void cvOneDMaterialLinear::SetAreas_and_length(double S_top,double S_bottom,double z){
  Stop = S_top;//inlet
  Sbot = S_bottom;//outlet
  len = z;
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
  
  return r;
}


double cvOneDMaterialLinear::GetDr1Dz(double z) const{
  double r_top = sqrt(Stop/M_PI);
  double r_bot = sqrt(Sbot/M_PI);
  double drodz = ((r_bot - r_top)/len) ; // These values are the initial radii.

  return drodz;
}


double cvOneDMaterialLinear::GetArea(double pressure, double z)const{
  // NOTE: o "So_" is the LSA under pressure p1_.
  //         This property comes from the subdomain
  //
  //       o Po is the the zero transmural pressure
  double pres = pressure;
  double So_  = GetS1(z);
  double EHR  = GetEHR(z);//*4/3
  
  // This is the area computation using the "pressure-strain" modulus, EHR.
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
  double press = p1_ + EHR*(sqrt(S/So_)-1.0);// for linear model dynes/cm^2
  return press;
}


double cvOneDMaterialLinear::GetDpDS(double S, double z)const{
  double EHR = GetEHR(z);
  double So_ = GetS1(z);
  double ro  = Getr1(z);
  double dpds=0.5* EHR/sqrt(So_*S) ;// for linear model
  return dpds;
}


// Careful!! this D2p(S,z)Dz first derivative, 2nd variable
double cvOneDMaterialLinear::GetDpDz(double S, double z)const{
  double So_   = GetS1(z);
  double EHR   = GetEHR(z);
  double r     = sqrt(S/M_PI);
  double ro    = Getr1(z);
  double drodz = GetDr1Dz(z); // if straight tube ->0.0 // Returns dro/dz
  double dpdz = drodz*(-EHR*r/ro/ro); 
  
  return dpdz;
}

