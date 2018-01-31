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
    double GetPsi( double pressure, double z) const;//not used
    double GetPhi( double pressure, double z) const; //not used
    double GetDSDp( double pressure, double z) const; //not used
    double GetDpDS( double area, double z) const;
    double GetD2pDS2( double area, double z) const;
    double GetDD2PDzDS( double area, double z) const; //IV added 03-23-03 for viscosity term 
    double GetOutflowFunction( double pressure, double z) const;
    double GetDOutflowDp( double pressure, double z) const;
    double GetD2S( double pressure, double z) const; //not used
    double GetIntegralpD2S ( double area, double z) const;//IV added 01-23-03
    double GetIntegralpS ( double area, double z) const;//IV added 01-24-03
    double GetDpDz( double area, double z) const;
    double GetTopArea() const {return Stop;}
    double GetBotArea() const {return Sbot;}
    void   SetEHR(double ehr_val);
    double GetEHR(double z) const;
    double GetMette2(double area,double z) const;
    double GetLinCompliance(double z) const;//IV added 02-03-03
    double GetnonLinCompliance(double area,double z) const;//IV added 02-13-03
    double GetWaveSpeed( double area, double z) const;//added by IV 080703
    double GetRefWaveSpeed( double area) const;//IV 080703 created for wave BC but can be used in gal; 
    double GetN(double S) const{return 0.0;};
    void   SetPeriod(double period){};
    
  private:

    double Stop;
    double Sbot;
    double len;

    double ehr;

    double GetS1( double z) const;
    double Getr1( double z) const;
    double GetDS1Dz( double z) const;
    double GetDr1Dz(double z) const;//was GetDrDz, but IV changed it 04-23-03 because it was misleading 
};

#endif // CVONEDMATERIALLINEAR_H

