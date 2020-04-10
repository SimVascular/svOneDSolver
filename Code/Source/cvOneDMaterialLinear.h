#ifndef CVONEDMATERIALLINEAR_H
#define CVONEDMATERIALLINEAR_H

//
//  cvOneDMaterialLinear.h - Header for a class to maintain material properties
//
//  This class maintains MaterialLinear Properties of the subdomain.
//

# include "cvOneDMaterial.h"
# include "cvOneDEnums.h"

class cvOneDMaterialLinear:public cvOneDMaterial{

  public:

    cvOneDMaterialLinear();
    ~cvOneDMaterialLinear();
    cvOneDMaterialLinear (const cvOneDMaterialLinear &rhs);
    cvOneDMaterialLinear& operator= ( const cvOneDMaterialLinear &that);

    void   SetAreas_and_length(double S_top, double S_bottom, double z);
    void   SetStop(double S){Stop = S;}
    void   SetSbottom(double S){Sbot = S;}
    void   SetLength(double length){len = length;}
    double GetProperty( char* what) const;
    double GetArea( double pressure, double z) const;
    double GetPressure( double S, double z) const;
    double GetDpDS( double area, double z) const;
    double GetD2pDS2( double area, double z) const;
    double GetDD2PDzDS( double area, double z) const; // for viscosity term
    double GetOutflowFunction( double pressure, double z) const;
    double GetDOutflowDp( double pressure, double z) const;
    double GetIntegralpD2S ( double area, double z) const;
    double GetIntegralpS ( double area, double z) const;
    double GetDpDz( double area, double z) const;
    double GetTopArea() const {return Stop;}
    double GetBotArea() const {return Sbot;}
    void   SetEHR(double ehr_val, double pref_val);
    double GetEHR(double z) const;
    double GetMette2(double area,double z) const;
    double GetLinCompliance(double z) const;
    double GetnonLinCompliance(double area,double z) const;
    double GetN(double S) const{return 0.0;};
    void   SetPeriod(double period){};

  private:

    double Stop;
    double Sbot;
    double len;

    double ehr;
    double PP1_;

    double GetS1( double z) const;
    double Getr1( double z) const;
    double GetDS1Dz( double z) const;
    double GetDr1Dz(double z) const;
};

#endif // CVONEDMATERIALLINEAR_H

