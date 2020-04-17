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


    virtual double GetN(double S) const = 0;//not really dependent on S actually IV 080703
    virtual double GetEHR(double z) const = 0;

    virtual double GetOutflowFunction(double pressure, double z) const = 0;
    virtual double GetDpDz(double area, double z) const = 0;
    virtual double GetDOutflowDp(double pressure, double z) const = 0;
    virtual double GetIntegralpD2S (double area, double z) const = 0;
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

