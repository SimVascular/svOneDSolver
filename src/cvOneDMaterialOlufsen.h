#ifndef CVONEDMATERIALOLUFSEN_H
#define CVONEDMATERIALOLUFSEN_H

//
//  cvOneDMaterialOlufsen.h - Header for a class to maintain material properties
//  
//  This class maintains MaterialOlufsen Properties of the subdomain.
//
//  History:
//  Feb. 13, 2003, I. Vignon
//      linear and nonlinear compliance.
//  Jul. 2000, J.Wan
//      Stenosis model properties
//  May 1999, J.Wan, G.R.Feijoo, S.A.Spicer and S.Strohband
//      Creation of file, class project of ME234C of T.J.R. Hughes and C.Taylor

# include "cvOneDMaterial.h"
# include "cvOneDEnums.h"

class cvOneDMaterialOlufsen:public cvOneDMaterial{
 
  public:

    cvOneDMaterialOlufsen();
    ~cvOneDMaterialOlufsen();
    cvOneDMaterialOlufsen (const cvOneDMaterialOlufsen &rhs );
    cvOneDMaterialOlufsen& operator= ( const cvOneDMaterialOlufsen &that);

    void   SetAreas_and_length(double S_top, double S_bottom, double z);
    void   SetStop(double S){Stop = S;}
    void   SetSbottom(double S){Sbot = S;}
    void   SetLength(double length){len = length;}
    void   SetMaterialType(double*);
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
    double GetEHR(double z) const;
    double GetMette2(double area,double z) const;
 	double GetN(double S) const;//not really dependent on S actually IV 080703
    

	int IsRigid() {return rigid;}
    
	double GetLinCompliance(double z) const;//IV added 02-03-03
	double GetnonLinCompliance(double area,double z) const;//IV added 02-13-03
	void SetPeriod(double period);//static
	
    double GetWaveSpeed( double area, double z) const;//added by IV 080703
    double GetRefWaveSpeed( double area) const;//IV 080703 created for wave BC but can be used in gal; 
    //compute wave speed at the ref. state, assumes  Mette's material model

    void GetParams(double *K1, double *K2, double *K3) const {*K1 = K1_; *K2 = K2_; *K3 = K3_;}
    
 private:

    double Stop;
    double Sbot;
    double len;

    double Period;

    double zstar;
    double alpha;
    double pcrit;

    double K1_;
    double K2_;
    double K3_;
	
	int rigid;

    double GetS1( double z) const;
    double Getr1( double z) const;
    double GetDS1Dz( double z) const;
    double GetDr1Dz(double z) const;//was GetDrDz, but IV changed it 04-23-03 because it was misleading
    
};

#endif // CVONEDMATERIALOLUFSEN_H

