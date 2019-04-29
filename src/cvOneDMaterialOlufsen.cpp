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

//#define K1_ 2.00e7
//#define K2_ -22.53
//#define K3_ 8.65e5

//#define K3_ 22.0e5 // not original
const double PI = 4.0 * atan(1.0);

using namespace std;

cvOneDMaterialOlufsen::cvOneDMaterialOlufsen(){
 // p1_ = 0.0;//debug wave IV 081903
 // note that this overwrites p1_ reference pressure set up the input file, when GetNewInstance is called for each vessel, p1_ is overwritten,
 //K1_, K2_ and k3_ will receover to the user defined values by *olfmat = *((cvOneDMaterialOlufsen*)(materials[matID]))
 // a modification to GetParams is made to reset p1_ otherwise p1_ won't take the value set by the input file if that is not 85mmHg. wgyang Dec 2018
  p1_ = 1333.2237 * 85.0; // initial p.. systemic side
 // p1_ = 1333.2237 * 20.0; // initial p..pulmonary side
  zstar = 70.0;
  alpha = 9.29e-3/1333.2237;
  pcrit = 33250.0;
  //Period=1.1;

  //default material values
  K1_ = 2.00e7;//1.99925e07;
  K2_ = -22.5267;//-22.53
 // K3_ = 8.65e5; //865251;

  rigid = 0;
 // printf("call cvOneMaterialsOlufsen p1_=%f K3_=%f \n",p1_,K3_);
}

cvOneDMaterialOlufsen::~cvOneDMaterialOlufsen(){
}

cvOneDMaterialOlufsen::cvOneDMaterialOlufsen (const cvOneDMaterialOlufsen &rhs){
 // this seems not used
  cvOneDMaterial::operator=(rhs);
  double k1 = 0,k2 = 0,k3 = 0;
  double pref=0;
  rhs.GetParams( &k1, &k2, &k3, &pref);
  K1_ = k1;
  K2_ = k2;
  K3_ = k3;
  p1_= pref;
 // printf("call cvOneDmaterialOlufsen &rhs\n");

}

cvOneDMaterialOlufsen& cvOneDMaterialOlufsen::operator= (const cvOneDMaterialOlufsen &that){
  if (this != &that) {
    cvOneDMaterial::operator=(that);
    double k1 = 0,k2 = 0,k3 = 0;

    double pref=0;
     // that.GetParams( &k1, &k2, &k3); // modified to reset p1_ to use the value given by the input file see line 25.  wgyang Dec 2018
    that.GetParams( &k1, &k2, &k3, &pref);
    K1_ = k1;
    K2_ = k2;
    K3_ = k3;
    p1_= pref;
  //  printf("call cvOneDMaterialOlufsen that this K3_=%f p1_=%f \n",K3_,p1_ );
  }
  return *this;
}

void cvOneDMaterialOlufsen::SetMaterialType(double *mType,double Pref){
  K1_ = mType[0];
  K2_ = mType[1];
  K3_ = mType[2];
  PP1_=Pref;
  cout<< "Setting material K's "<< K1_ <<" "<< K2_<<" "<< K3_<< " ..." << endl;
  cout<< "Setting reference Pressure "<< PP1_<<endl;
 //  printf("call SetMaterialType K3_ %f \n",K3_);
}

void cvOneDMaterialOlufsen::SetPeriod(double per){
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

double cvOneDMaterialOlufsen::GetEHR(double z)const{
  // double ro = sqrt(GetS1(z)/PI);
  double ro = Getr1(z);
  double ans =4./3.*( K1_*exp(K2_*ro) + K3_);//dyne/cm^2=g/cm/s^2
  return ans;
}

void cvOneDMaterialOlufsen::SetAreas_and_length(double S_top,double S_bottom,double z){
  Stop = S_top;//inlet
  Sbot = S_bottom;//outlet
  len = z;
  //  cout<<"length "<< len<<" areas "<< Stop<<" "<<Sbot<<endl;
}

double cvOneDMaterialOlufsen::GetS1(double z)const{
  // Okay...now I'm assuming that S1 is the initial area.
  // So we'll linearly interpolate the area
  double area;

  // the area/pressure results do not appear to be correct.  Want radius to change linearly, not area

  // double areaold = ((z-len)/(-len))*(Stop - Sbot) + Sbot;
  // return areaold;

  // doesn't appear to change solution
  double r=Getr1(z);
  area= r*r*PI;
  return area;  // slightly increased pressure/decreased area
}

double cvOneDMaterialOlufsen::Getr1(double z)const{
  // linearly interpolated r
  double r_top=sqrt(Stop/PI);
  double r_bot=sqrt(Sbot/PI);
  double r=((z-len)/(-len))*(r_top - r_bot) + r_bot;

  // mette's method
  //r=r_top*exp(z/len*log(r_bot/r_top));


  return r;
}

//not used
double cvOneDMaterialOlufsen::GetDS1Dz(double z)const{
  // Since vessel geometry will always be linear in
  // space, this is just a constant
  /*
  double der;
  der = ((Sbot - Stop)/(len));
  return der;
  */
  // this shouldn't change solution, but try using dr/dz and compute S
  double drdz=GetDr1Dz(z) ;
  double dsdr= 2.0*PI*Getr1(z);
  return dsdr*drdz;  // slightly increased pressure/decreased area
}

//this is in the reference state dr1dz IV 04-23-03
double cvOneDMaterialOlufsen::GetDr1Dz(double z)const{
  // linearly vary area
  // double s = GetS1(z);
  // double dsdz=((Sbot - Stop)/(len));
  // double olddrdz= dsdz/(2*sqrt(PI*s));

  // linearly vary radius
  double r_top=sqrt(Stop/PI);
  double r_bot=sqrt(Sbot/PI);
  double drdz=((r_bot - r_top)/len) ;

  //mette's method
  //drdz=1/len*log(r_bot/r_top)*Getr1(z);

  return drdz;
}

//IV 080703 can work anywhere in the computational domain but not in the downstream domain
double cvOneDMaterialOlufsen::GetWaveSpeed(double area, double z)const{   double c = sqrt(area * GetDpDS(area,z)/density);
  return c;
}

//IV 080703 created for wave BC but can be used in gal; compute the wave speed with the reference state
//assumes  Mette's material model
double cvOneDMaterialOlufsen::GetRefWaveSpeed(double area)const{
  double ro = sqrt(area/PI);
  double EHRo = 4./3.*( K1_*exp(K2_*ro) + K3_);
  double co = sqrt(EHRo/2/density);
  return co;
}

double cvOneDMaterialOlufsen::GetArea(double pressure, double z)const{
  // NOTE: o "So_" is the LSA under pressure p1_.
  //         This property comes from the subdomain
  //
  //       o Po is the the zero transmural pressure

  double pres = pressure;
  double So_  = GetS1(z);
  double EHR  = GetEHR(z);//*4/3
  // double area1 = (16.0*So_*pow(EHR,2))/pow(-3.0*pres + 3.0*p1_ + 4.0*EHR,2);//I moved 4/3 into EHR term, so this will be off now, but did work.
 // printf("GetArea press=%f z=%f So= %f K3_=%f, p1_=%f \n",pres,z,So_,K3_,p1_);
  double area = So_/pow(1-(pres-p1_)/EHR,2);// this makes more sense to me, bns 8/22/02. gives same answer
  return area;
}

double cvOneDMaterialOlufsen::GetPressure(double S, double z)const{
  // Again we need to get So_ from the subdomain.
  // Then we implement Olufsen's constitutive law...
  double So_   = GetS1(z);
  double EHR   = GetEHR(z);  // From Olufsen's paper
  double press = p1_ + EHR*(1.-sqrt(So_/S));//dynes/cm^2
//  printf("GetPressure press=%f z=%f So= %f S=%f, K3_=%f, p1_=%f \n",press,z,So_,S,K3_,p1_);
  // mette's method - increases pressure
  // double theta = PI/2.*(S/So_ - 1.);
  // press = p1_ + EHR/PI*tan(theta);//I added p1_

  return press;
}

double cvOneDMaterialOlufsen::GetDpDS(double S, double z)const{
  double EHR = GetEHR(z);
  double So_ = GetS1(z);
  double ro=Getr1(z);
  double dpds=0.5* EHR * sqrt(So_/S)/S ;

  // mette's version
  // double theta = PI/2*(S/So_ -1);
  // double p=EHR/PI;
  // dpds = 0.5*p/(ro*cos(theta)*ro*cos(theta));

  return dpds;
}

// not used //0.75 added by IV 04-23-03, checked with matlab
double cvOneDMaterialOlufsen::GetD2pDS2(double area, double z)const{
  double EHR = GetEHR(z);
  double So_ = GetS1(z);
  return - 0.75 * EHR * sqrt(So_) / sqrt(pow(area, 5));
}

//used for viscosity term in matrix outlet flux term , checked
double cvOneDMaterialOlufsen::GetDD2PDzDS(double area, double z)const{
  double EHR   = GetEHR(z);
  double ro    = Getr1(z);
  double drodz = GetDr1Dz(z);
  double dEHRdro = 4./3.*K2_*K1_*exp(K2_*ro);
  double derP = drodz*sqrt(PI)*0.5/area/sqrt(area)*(dEHRdro*ro+EHR);
  return derP;
}

double cvOneDMaterialOlufsen::GetOutflowFunction(double pressure, double z)const{
  return 0.0; // This is not used in our model
}

double cvOneDMaterialOlufsen::GetDOutflowDp(double pressure, double z)const{
  return 0.0; // Nor is this.
}

//IV added it 01-23-03
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

//IV added it 01-24-03
double cvOneDMaterialOlufsen::GetIntegralpS(double area, double z)const{
  double EHR   = GetEHR(z);
  double So_   = GetS1(z);
  double IntegralpS = EHR*So_*(sqrt(area/So_)-1.);
  //cout<<EHR<<" ";
  //cout<<IntegralpS<<endl;
  return IntegralpS;
}

// Careful!! this D2p(S,z)Dz first derivative, 2nd variable IV 03-31-03
double cvOneDMaterialOlufsen::GetDpDz(double S, double z)const{
  double So_   = GetS1(z);
  double EHR   = GetEHR(z);
  double r     = sqrt(S/PI);
  double ro    = Getr1(z);
  double drodz  = GetDr1Dz(z);//if straight tube ->0.0
  double dEHRdr = 4./3.*K2_*K1_*exp(K2_*ro); // dfdr
  double dpdz = drodz*(dEHRdr*(1.0-(ro/r))-EHR/r);

  // mette's causes pressure to increase in longer tubes
  // double theta= PI/2.*(S/So_ - 1.);
  // dEHRdr=4./3.*K2_*K1_*exp(K2_*ro)/PI;
  // double p=EHR/PI;
  // dpdz=drdz*(dEHRdr*tan(theta)- p*S/(ro*ro*ro*cos(theta)*cos(theta)));

  return dpdz; //careful this is D2P(S,z)Dz
}

// Want to see what happens if we use boundary layer flow profile instead of parabolic flow profile->not much
// From olufsen
double cvOneDMaterialOlufsen::GetN(double S)const{
  double R    = sqrt(S/PI);
  // mette uses period=T/10
  double T = Period;//*10;// from non-dim,
  double del  = sqrt(kinematicViscosity*T/2.0/PI);// from Lighthill, olufsen's thesis, pg56
  double altN = -2.0*PI*kinematicViscosity*R/del; // == -1.45838*r
  // N = -2.0*PI*kinematicViscosity*(profile_exponent+2.0) ; // == -1.1673

  return N;
  //  return altN;

}

//IV added it 02-03-03
double cvOneDMaterialOlufsen::GetLinCompliance(double z)const{
  double EHR   = GetEHR(z);
  double S_o   = GetS1(z);
  double LinCompliance = 2.0*S_o/EHR;//=dSdP at po,S_o checked IV 02/13/03
  return LinCompliance;
}

//IV added it 02-13-03
double cvOneDMaterialOlufsen::GetnonLinCompliance(double area, double z)const{
  double EHR   = GetEHR(z);
  double S_o   = GetS1(z);
  double ratioA= sqrt(area/S_o);
  double nonLinCompliance = 2.0*S_o/EHR*ratioA*ratioA*ratioA;
  return nonLinCompliance;
}

