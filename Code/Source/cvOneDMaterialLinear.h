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
    void   SetHydraulicConductivity(double value) {L_P = value;};
    void   SetStarlingAmbientPressure(double value) {P_ambient = value;};
    double GetStarlingAmbientPressure() {return P_ambient;}
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
    double GetN(double S) const{return N;}; // JR 15/11/23: not sure why this returned 0.0 instead of just returning N. This has now bee fixed.
    void   SetPeriod(double period){};

  private:

    double Stop;
    double Sbot;
    double len;

    double L_P;
    double P_ambient;

    double ehr;
    double PP1_;

    double GetS1( double z) const;
    double Getr1( double z) const;
    double GetDS1Dz( double z) const;
    double GetDr1Dz(double z) const;
};

#endif // CVONEDMATERIALLINEAR_H

