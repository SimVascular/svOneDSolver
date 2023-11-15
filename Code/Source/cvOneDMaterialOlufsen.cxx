/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//
//  MaterialOlufsen.cxx - Source for a class to maintain material properties
//
//  SYNOPSIS...This class maintains MaterialOlufsen Properties of the subdomain.
//

# include <cstring>
# include <cstdlib>
# include <cmath>
# include <iostream>

# include "cvOneDMaterialOlufsen.h"
# include "cvOneDUtility.h"

const double PI = 4.0 * atan(1.0);

using namespace std;

cvOneDMaterialOlufsen::cvOneDMaterialOlufsen(){
  cout << "cvOneDMaterialOlufsen function1: P_ambient is: " << P_ambient << endl;
  p1_ = 1333.2237 * 85.0; // initial p.. systemic side
 // p1_ = 1333.2237 * 20.0; // initial p..pulmonary side
  zstar = 70.0;
  alpha = 9.29e-3/1333.2237;
  pcrit = 33250.0;

  //default material values
  K1_ = 2.00e7;
  K2_ = -22.5267;
 // K3_ = 8.65e5;

  rigid = 0;
 // printf("call cvOneMaterialsOlufsen p1_=%f K3_=%f \n",p1_,K3_);
}

cvOneDMaterialOlufsen::~cvOneDMaterialOlufsen(){
  cout << "~cvOneDMaterialOlufsen function: P_ambient is: " << P_ambient << endl;
}

cvOneDMaterialOlufsen::cvOneDMaterialOlufsen (const cvOneDMaterialOlufsen &rhs){
  cout << "cvOneDMaterialOlufsen function2: P_ambient is: " << P_ambient << endl;
 // this seems not used
  cvOneDMaterial::operator=(rhs);
  double k1 = 0,k2 = 0,k3 = 0;
  double pref=0; double P_amb_ref=0; double L_P_ref = 0;
  rhs.GetParams( &k1, &k2, &k3, &pref, &P_amb_ref, &L_P_ref);
  K1_ = k1;
  K2_ = k2;
  K3_ = k3;
  p1_= pref;
  P_ambient = P_amb_ref;
  L_P = L_P_ref;
 // printf("call cvOneDmaterialOlufsen &rhs\n");

}

cvOneDMaterialOlufsen& cvOneDMaterialOlufsen::operator= (const cvOneDMaterialOlufsen &that){
  cout << "operator= function: P_ambient is: " << P_ambient << endl;
  if (this != &that) {
    cvOneDMaterial::operator=(that);
    double k1 = 0,k2 = 0,k3 = 0;

    double pref=0; double P_amb_ref=0; double L_P_ref = 0;
    that.GetParams( &k1, &k2, &k3, &pref, &P_amb_ref, &L_P_ref);
    K1_ = k1;
    K2_ = k2;
    K3_ = k3;
    p1_= pref;
    P_ambient = P_amb_ref;
    L_P = L_P_ref;
  //  printf("call cvOneDMaterialOlufsen that this K3_=%f p1_=%f \n",K3_,p1_ );
  }
  return *this;
}

void cvOneDMaterialOlufsen::SetMaterialType(double *mType,double Pref){
  cout << "SetMaterialType function: P_ambient is: " << P_ambient << endl;
  K1_ = mType[0];
  K2_ = mType[1];
  K3_ = mType[2];
  PP1_=Pref;
  cout<< "Setting material K's "<< K1_ <<" "<< K2_<<" "<< K3_<< " ..." << endl;
  cout<< "Setting reference Pressure "<< PP1_<<endl;
 //  printf("call SetMaterialType K3_ %f \n",K3_);
}

void cvOneDMaterialOlufsen::SetPeriod(double per){
  cout << "SetPeriod function: P_ambient is: " << P_ambient << endl;
  Period=per;
}

double cvOneDMaterialOlufsen::GetProperty(char* what)const{
  // Nothing to change here...
  if( strcmp( what, "density") == 0)
    return density;
  else if( strcmp( what, "dynamic viscosity") == 0)
    return dynamicViscosity;
  else if( strcmp( what, "kinematic viscosity") == 0)
    return kinematicViscosity;
  else if( strcmp( what, "constant n") == 0)
    return profile_exponent;
  else if( strcmp( what, "delta") == 0)
    return delta;
  else if( strcmp( what, "N") == 0)
    return N;
  else{
    abort();
    return 0.0;
  }
}

void cvOneDMaterialOlufsen::SetHydraulicConductivity(double value) {
  L_P = value;
  cout << "L_P: "<< L_P << endl;
}

void cvOneDMaterialOlufsen::SetStarlingAmbientPressure(double value) {
  cout << "value: " << value << endl;
  P_ambient = value; 
  cout << "P_ambient: "<< P_ambient << endl;
}

double cvOneDMaterialOlufsen::GetStarlingAmbientPressure() {
  cout << "GetStarlingAmbientPressure function: P_ambient is: " << P_ambient << endl;
  return P_ambient;
}

double cvOneDMaterialOlufsen::GetEHR(double z)const{
  double ro = Getr1(z);
  double ans =4./3.*( K1_*exp(K2_*ro) + K3_);//dyne/cm^2=g/cm/s^2
  return ans;
}

void cvOneDMaterialOlufsen::SetAreas_and_length(double S_top,double S_bottom,double z){
  Stop = S_top;//inlet
  Sbot = S_bottom;//outlet
  len = z;
}

double cvOneDMaterialOlufsen::GetS1(double z)const{
  double area;
  double r=Getr1(z);
  area= r*r*PI;
  return area;  // slightly increased pressure/decreased area
}

double cvOneDMaterialOlufsen::Getr1(double z)const{
  // linearly interpolated r
  double r_top=sqrt(Stop/PI);
  double r_bot=sqrt(Sbot/PI);
  double r=((z-len)/(-len))*(r_top - r_bot) + r_bot;


  return r;
}

double cvOneDMaterialOlufsen::GetDS1Dz(double z)const{
  double drdz=GetDr1Dz(z) ;
  double dsdr= 2.0*PI*Getr1(z);
  return dsdr*drdz;  // slightly increased pressure/decreased area
}


//this is in the reference state dr1dz
double cvOneDMaterialOlufsen::GetDr1Dz(double z)const{
  // linearly vary radius
  double r_top=sqrt(Stop/PI);
  double r_bot=sqrt(Sbot/PI);
  double drdz=((r_bot - r_top)/len) ;

  return drdz;
}


double cvOneDMaterialOlufsen::GetArea(double pressure, double z)const{
  // NOTE: o "So_" is the LSA under pressure p1_.
  //         This property comes from the subdomain
  //
  //       o Po is the the zero transmural pressure
  double pres = pressure;
  double So_  = GetS1(z);
  double EHR  = GetEHR(z);
  double area = So_/pow(1-(pres-p1_)/EHR,2);
  return area;
}

double cvOneDMaterialOlufsen::GetPressure(double S, double z)const{
  // Again we need to get So_ from the subdomain.
  // Then we implement Olufsen's constitutive law...
  double So_   = GetS1(z);
  double EHR   = GetEHR(z);  // From Olufsen's paper
  double press = p1_ + EHR*(1.-sqrt(So_/S)); // dynes/cm^2

  return press;
}

double cvOneDMaterialOlufsen::GetDpDS(double S, double z)const{
  double EHR = GetEHR(z);
  double So_ = GetS1(z);
  double ro=Getr1(z);
  double dpds=0.5* EHR * sqrt(So_/S)/S ;
  cout << EHR << " " << So_ << " " << ro << endl;
  return dpds;
}

double cvOneDMaterialOlufsen::GetD2pDS2(double area, double z)const{
  double EHR = GetEHR(z);
  double So_ = GetS1(z);
  return - 0.75 * EHR * sqrt(So_) / sqrt(pow(area, 5));
}

double cvOneDMaterialOlufsen::GetOutflowFunction(double pressure, double z)const{
  return L_P*(pressure - P_ambient); // JR 10/11/23: added function for outflow term
}

double cvOneDMaterialOlufsen::GetDOutflowDp(double pressure, double z)const{
  return 0.0; // Nor is this.
}

//used for viscosity term in matrix outlet flux term
double cvOneDMaterialOlufsen::GetDD2PDzDS(double area, double z)const{
  double EHR   = GetEHR(z);
  double ro    = Getr1(z);
  double drodz = GetDr1Dz(z);
  double dEHRdro = 4./3.*K2_*K1_*exp(K2_*ro);
  double derP = drodz*sqrt(PI)*0.5/area/sqrt(area)*(dEHRdro*ro+EHR);
  return derP;
}


double cvOneDMaterialOlufsen::GetIntegralpD2S(double area, double z)const{
  double EHR   = GetEHR(z);
  double DS1Dz = GetDS1Dz(z);
  double So_   = GetS1(z);
  double DroDz = GetDr1Dz(z);
  double ro    = Getr1(z);
  double DEHRinvDr = -4./3.*K1_*K2_*exp(K2_*ro); //should be /(EHR)^2 but included in IntegralpD2S
  double DEHRinvDz = DEHRinvDr*DroDz;
  double termA = sqrt(area/So_)-1.0;
  double IntegralpD2S = DS1Dz*EHR*termA+So_*termA*termA*DEHRinvDz;
  return IntegralpD2S;
}

double cvOneDMaterialOlufsen::GetIntegralpS(double area, double z)const{
  double EHR   = GetEHR(z);
  double So_   = GetS1(z);
  double IntegralpS = EHR*So_*(sqrt(area/So_)-1.);

  return IntegralpS;
}

// Careful!! this D2p(S,z)Dz first derivative, 2nd variable
double cvOneDMaterialOlufsen::GetDpDz(double S, double z)const{
  double So_   = GetS1(z);
  double EHR   = GetEHR(z);
  double r     = sqrt(S/PI);
  double ro    = Getr1(z);
  double drodz  = GetDr1Dz(z); // if straight tube ->0.0
  double dEHRdr = 4./3.*K2_*K1_*exp(K2_*ro); // dfdr
  double dpdz = drodz*(dEHRdr*(1.0-(ro/r))-EHR/r);

  return dpdz; //careful this is D2P(S,z)Dz
}

double cvOneDMaterialOlufsen::GetN(double S)const{
  double R    = sqrt(S/PI);
  double T = Period;
  double del  = sqrt(kinematicViscosity*T/2.0/PI); // from Lighthill, olufsen's thesis, pg56
  double altN = -2.0*PI*kinematicViscosity*R/del; // == -1.45838*r

  return N;

}

double cvOneDMaterialOlufsen::GetLinCompliance(double z)const{
  double EHR   = GetEHR(z);
  double S_o   = GetS1(z);
  double LinCompliance = 2.0*S_o/EHR; // =dSdP at po,S_o
  return LinCompliance;
}

double cvOneDMaterialOlufsen::GetnonLinCompliance(double area, double z)const{
  double EHR   = GetEHR(z);
  double S_o   = GetS1(z);
  double ratioA= sqrt(area/S_o);
  double nonLinCompliance = 2.0*S_o/EHR*ratioA*ratioA*ratioA;
  return nonLinCompliance;
}

