#ifndef CVONEDMATERIAL_H
#define CVONEDMATERIAL_H

//
//  cvOneDMaterial.h - Header for a class to maintain material properties
//
//  This class maintains Material Properties of the subdomain.
//

# include "cvOneDEnums.h"
# include <math.h>

class cvOneDMaterial{

  public:

    cvOneDMaterial(){
      // using CGS units - default values
      density = 1.06;
      dynamicViscosity = 0.04;
      this->UpdateKinematicViscosity();
      this->SetProfileExponent(2);
    }

    // Copy Constructor
    cvOneDMaterial(const cvOneDMaterial &rhs ){
      density = rhs.GetDensity();
      dynamicViscosity = rhs.GetDynamicViscosity();
      this->UpdateKinematicViscosity();
      this->SetProfileExponent(rhs.GetProfileExponent());
    }

    // Destructor
    virtual ~cvOneDMaterial(){}

    cvOneDMaterial& operator= (const cvOneDMaterial& that){
      if(this != &that){
        density = that.GetDensity();
        dynamicViscosity = that.GetDynamicViscosity();
        this->UpdateKinematicViscosity();
        this->SetProfileExponent(that.GetProfileExponent());
      }
      return *this;
    }

    virtual double GetArea(double pressure, double z) const = 0;
    virtual double GetPressure(double S, double z) const = 0;
    virtual double GetDpDS(double area, double z) const = 0;
    virtual double GetD2pDS2(double area, double z) const = 0;

    virtual double GetProperty(char* what) const = 0;
    virtual double GetIntegralpS (double area, double z) const = 0;

    virtual double GetN(double S) const = 0; //not really dependent on S actually 
    virtual double GetOutflowFunction(double pressure, double z) const = 0;
    virtual double GetDpDz(double area, double z) const = 0;
    virtual double GetDOutflowDp(double pressure, double z) const = 0;
    virtual double GetIntegralpD2S (double area, double z) const = 0;
    virtual double GetRefWaveSpeed(double area) const = 0;
    virtual void   SetPeriod(double period) = 0;
    virtual void   SetAreas_and_length(double S_top, double S_bottom, double z) = 0;

    double GetProfileExponent() const {return profile_exponent;}
    double GetDensity() const {return density;}
    double GetDynamicViscosity() const {return dynamicViscosity;}
    void   SetProfileExponent(double value) {profile_exponent = value;
                                             delta = 1.0/(1.0+profile_exponent);
                                             double myPI = 4.0 * atan(1.0);
                                             N = -2.0*myPI*kinematicViscosity*(profile_exponent+2.0);
                                             return;}
    void   SetDensity(double value) {density = value;this->UpdateKinematicViscosity();return;}
    void   SetDynamicViscosity(double value) {dynamicViscosity = value;this->UpdateKinematicViscosity();return;}

    void   UpdateKinematicViscosity() {kinematicViscosity = dynamicViscosity/density;return;}

    void   SetReferencePressure(double refP) {p1_ = refP;}
    double GetReferencePressure() {return p1_; }
    double GetReferencedPressure_dt() {return 0; }
    
  protected:

    double density;
    double dynamicViscosity;
    double kinematicViscosity;

    double profile_exponent;
    double delta;
    double N;

    double p1_;

};

#endif // CVONEDMATERIAL_H

