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
//  Subdomain.cxx - Source for a class to contain the discretization of
//  ~~~~~~~~~~~~~   of the Geometry.
//
//  This class will contain the discretization of the geometry.
//

# include "cvOneDGlobal.h"
# include "cvOneDSubdomain.h"
# include "cvOneDMaterialManager.h"

const double PI = 4.0 * atan(1.0);//IV 080703

using namespace std;

cvOneDSubdomain::cvOneDSubdomain(){
  mat = NULL;
  pressureWave = NULL;
  pressureTime = NULL;
  resistanceWave = NULL;
  resistanceTime = NULL;
  numPressurePts = 0;
  // Coronary boundarycondition
  PressLVWave=NULL;
  PressLVTime=NULL;
  numPressLVPts=0;
  branchAngle = 90.0;
}

cvOneDSubdomain::~cvOneDSubdomain(){
  if(connectivities != NULL)  delete [] connectivities;
  if(nodes != NULL)      delete [] nodes;
  if(pressureWave != NULL)  delete [] pressureWave;
  if(pressureTime != NULL)  delete [] pressureTime;
  if(resistanceWave != NULL)  delete [] resistanceWave;
  if(resistanceTime != NULL)  delete [] resistanceTime;
  if(presslv != NULL) delete [] presslv;// Added by Jongmin Seo 04062020 & Hyunjin Kim 09022005
  if(PressLVWave != NULL)  delete [] PressLVWave;// Added by Jongmin Seo 04062020 & Hyunjin Kim 09022005
  if(PressLVTime != NULL)  delete [] PressLVTime;// Added by Jongmin Seo 04062020 & Hyunjin Kim 09022005
}

void cvOneDSubdomain::SetInitInletS(double So){
  S_initial = So;
}

void cvOneDSubdomain::SetInitOutletS(double Sn){
  S_final = Sn;
}

void cvOneDSubdomain::SetInitialFlow(double Qo){
  Q_initial = Qo;
}

void cvOneDSubdomain::SetInitialPressure(double Po){
  P_initial = Po;
}

void cvOneDSubdomain::SetInitialdFlowdT(double dQ0dT){
  //dQ_dT_initial = dQ0dT;
  dQ_dT_initial = 0;
}

void cvOneDSubdomain::SetupMaterial(int matID){
  mat = cvOneDGlobal::gMaterialManager->GetNewInstance(matID);
 // printf("subdomain cpp setupMaterial matID=%i  \n", matID);
  mat->SetAreas_and_length(S_initial, S_final, fabs(z_out - z_in));
  cout << "line 98 of cvOneDSubdomain.cxx" << endl;
  cout << mat->GetEHR((double) 0.3) << endl;
  // mat->GetStarlingAmbientPressure();
}

void cvOneDSubdomain::SetBoundValue(double boundV){
  switch(boundType){
    case BoundCondTypeScope::PRESSURE:
      boundValue = mat->GetArea(boundV, fabs(z_out-z_in));
      break;
    default:
      boundValue = boundV;
      break;
  }
}

double cvOneDSubdomain::GetInitInletS(void) {return S_initial;}
double cvOneDSubdomain::GetInitOutletS(void) {return S_final;}
double cvOneDSubdomain::GetInitialFlow(void) {return Q_initial;}
double cvOneDSubdomain::GetInitialPressure(void) {return P_initial;}

//  Added by Jongmin Seo 04062020 & Hyunjin Kim 09022005
double cvOneDSubdomain::GetInitialdFlowdT(void){
  return dQ_dT_initial;
  // return 0;
}

void cvOneDSubdomain::SetNumberOfNodes(long nNodes){
  numberOfNodes = nNodes;
}

void cvOneDSubdomain::SetNumberOfElements(long nElements){
  numberOfElements = nElements;
}

void cvOneDSubdomain::SetMeshType(MeshType mType){
  meshType = mType;
}

void cvOneDSubdomain::Init(double x0, double xL){
  z_in = x0;
  z_out = xL;
  char errStr[256];
  nodes = new double[numberOfNodes];
  connectivities = new long[ 2 * numberOfElements];
  assert( nodes != 0 && connectivities != 0);

  // Only handle uniform mesh for now...
  if(meshType != MeshTypeScope::UNIFORM){
    cvOneDError::setErrorNumber(ErrorTypeScope::UNSUPPORTED);
    strcpy(errStr,"In Subdomain::Init(...), Unsuported mesh type requested, continuing with a uniform mesh.");
    cvOneDError::setErrorString(errStr);
    cvOneDError::CallErrorHandler();
  }

  // Set the mesh to be uniform
  double h = (xL - x0) / static_cast<double>(numberOfElements);

  for( long i = 0; i < numberOfNodes; i++){
    nodes[i] = x0 + i * h;
  }

  long nd = 0;
  for( long element = 0; element < numberOfElements; element++, nd++){
    connectivities[ 2 * element    ] = nd;
    connectivities[ 2 * element + 1] = nd + 1;

  }

  finiteElement = new cvOneDFiniteElement();
  assert( finiteElement != 0);
}

long cvOneDSubdomain::GetNumberOfNodes()const{
  return numberOfNodes;
}

long cvOneDSubdomain::GetNumberOfElements()const{
  return numberOfElements;
}

void cvOneDSubdomain::GetConnectivity(long element, long* conn)const{
   conn[0] = connectivities[ 2 * element    ] + global1stNodeID;
   conn[1] = connectivities[ 2 * element + 1] + global1stNodeID;
}

void cvOneDSubdomain::GetNodes( long element, double* nd)const{
  nd[0] = nodes[ connectivities[2*element]];
  nd[1] = nodes[ connectivities[2*element+1]];
}

cvOneDFiniteElement* cvOneDSubdomain::GetElement(long element)const{
  long conn[2];
  double nd[2];
  GetConnectivity( element, conn);
  GetNodes( element, nd);
  finiteElement->Set( nd, conn);
  return finiteElement;
}

double cvOneDSubdomain::GetNodalCoordinate( long node)const{
  return nodes[node];
}

void cvOneDSubdomain::SetBoundPresWave(double *time, double *pres, int num){
  int i ;
  pressureWave = new double[num];
  pressureTime = new double[num];
  for(i = 0; i < num; i++){
      pressureTime[i] = time[i];
      pressureWave[i] = pres[i];
  }
  numPressurePts = num;
}

void cvOneDSubdomain::SetBoundResistanceWave(double *time, double *resist, int num){
  int i ;
  resistanceWave = new double[num];
  resistanceTime = new double[num];
  for(i = 0; i < num; i++){
      resistanceTime[i] = time[i];
      resistanceWave[i] = resist[i];
  }
  numPressurePts = num;
}

// Added by Jongmin Seo 04062020 & Hyunjin Kim 09022005
double cvOneDSubdomain::getBoundCoronaryValues(double currentTime){
  if(PressLVTime == NULL || PressLVWave == NULL){
      cout << "ERROR: LV pressure information is not prescribed"<< endl;
      exit(1);
  }
  // flow rate is assumed to be periodic
  double cycleTime = PressLVTime[numPressLVPts-1];
  double correctedTime = currentTime - static_cast<long>(currentTime / cycleTime) * cycleTime;
  double presslv=0.;
  if(correctedTime >= 0 && correctedTime <= PressLVTime[0]){
    double xi = (correctedTime - PressLVTime[numPressLVPts-1]) / (PressLVTime[0] - PressLVTime[numPressLVPts-1]);
    presslv = (PressLVWave[numPressLVPts-1] + xi * (PressLVWave[0] - PressLVWave[numPressLVPts-1]));
  }else{
    int ptr = 0;
    bool wasFound = false;
    while( !wasFound){
      if(correctedTime >= PressLVTime[ptr] && correctedTime <= PressLVTime[ptr+1]){
        wasFound = true;
      }else{
        ptr++;
      }
    }

    // linear interpolation between values
    double xi = (correctedTime - PressLVTime[ptr]) / (PressLVTime[ptr+1] - PressLVTime[ptr]);

    presslv = (PressLVWave[ptr] + xi * (PressLVWave[ptr+1] - PressLVWave[ptr]));
  }
  return presslv;
}

// Added by Jongmin Seo 04062020 & Hyunjin Kim 09022005
void cvOneDSubdomain::SetBoundCoronaryValues(double *time, double *p_lv, int num){

  int numpts=num;
  corTime=0.0;
  //a=Ra1*Ra2*Ca*Cc;
  //b=Ca*Ra1+Cc*(Ra1+Ra2);
  //expo1=(-b+sqrt(b*b-4*a))/2/a;
  //expo2=(-b-sqrt(b*b-4*a))/2/a;
  if (time[0]< 0.0) {Ra1=p_lv[0];}
  if (time[1]< 0.0) {Ra2=p_lv[1];}
  if (time[2]< 0.0) {Ca=p_lv[2];}
  if (time[3]< 0.0) {Cc=p_lv[3];}
  if (time[4]< 0.0) {Rv1=p_lv[4];}
  if (time[5]< 0.0) {P_v=p_lv[5];}
  if (time[6]< 0.0) {
    fprintf(stdout, "Wrong file format");
    exit(1);
  }
  numpts=num-6;
  numPressLVPts=numpts;
  PressLVWave=new double[numpts];
  PressLVTime=new double[numpts];
  for (int i=0;i<numpts;i++){
    PressLVTime[i]=time[i+6];
    PressLVWave[i]=p_lv[i+6];
//    fprintf(stdout,"Time[%i]: %le Pressure_LV[%i]: %le\n",i, PressLVTime[i], i,PressLVWave[i]);
  }
  p0COR=1;
  p1COR=Ra2*Ca+(Rv1)*(Ca+Cc);
  p2COR=Ca*Cc*Ra2*(Rv1);
  q0COR=Ra1+Ra2+Rv1;
  q1COR=Ra1*Ca*(Ra2+Rv1)+Cc*(Rv1)*(Ra1+Ra2);
  q2COR=Ca*Cc*Ra1*Ra2*(Rv1);
  b0COR=0;
  b1COR=Cc*(Rv1);
  detCOR=sqrt(q1COR*q1COR-4*q0COR*q2COR);
  expo2COR=-(q1COR+detCOR)/2/q2COR;
  expo1COR=q0COR/q2COR/expo2COR;
  CoefZ1=(p0COR+p1COR*expo1COR+p2COR*expo1COR*expo1COR)/detCOR;
  CoefY1=-(b0COR+b1COR*expo1COR)/detCOR;
  CoefZ2=(p0COR+p1COR*expo2COR+p2COR*expo2COR*expo2COR)/detCOR;
  CoefY2=-(b0COR+b1COR*expo2COR)/detCOR;
  CoefR=p2COR/q2COR;
}


void cvOneDSubdomain::SetBoundRCRValues(double *rcr, int num){//added IV 050803

  rcrTime = 0.0;
  rcrTime2=0.0;
  assert(num>=3);
  proximalResistance = rcr[0];
  capacitance = rcr[1] ;
  distalResistance= rcr[2];
  alphaRCR = (proximalResistance+distalResistance)/(proximalResistance*distalResistance*capacitance);
  if (num>3){
  Pd=rcr[3];
  }else {
  Pd=0.0;
  }
}


void cvOneDSubdomain::SetBoundResistPdValues(double *value, int num){//resistance with Pd wgyang 2019/4
  assert(num<=2);
  resistancevalue=value[0];
  if (num==2){
  Pd=value[1];
  }else{
  Pd=0.0;//if no Pd is provided
  }


}


// Returns interpolated pressure from specific boundary pressure wave
double cvOneDSubdomain::GetPressure(double currentTime){
  if (currentTime == 0){
    double pressure= P_initial;
    return pressure;
  }

  if(pressureTime == NULL || pressureWave == NULL){
    cout << "ERROR: pressure information is not prescribed"<< endl;
    exit(1);
  }

  // flow rate is assumed to be periodic
  double cycleTime = pressureTime[numPressurePts-1];
  double correctedTime = currentTime - static_cast<long>(currentTime / cycleTime) * cycleTime;

  double pressure=0.;
  if(correctedTime >= 0 && correctedTime <= pressureTime[0]){
    double xi = (correctedTime - pressureTime[numPressurePts-1]) / (pressureTime[0] - pressureTime[numPressurePts-1]);
    pressure = (pressureWave[numPressurePts-1] + xi * (pressureWave[0] - pressureWave[numPressurePts-1]));
  }else{
    int ptr = 0;
    bool wasFound = false;
    while( !wasFound){
      if( correctedTime >= pressureTime[ptr] && correctedTime <= pressureTime[ptr+1]){
        wasFound = true;
      }else{
        ptr++;
      }
    }

    // linear interpolation between values
    double xi = (correctedTime - pressureTime[ptr]) / (pressureTime[ptr+1] - pressureTime[ptr]);
    pressure = (pressureWave[ptr] + xi * (pressureWave[ptr+1] - pressureWave[ptr]));
  }
  return pressure;
}

// Returns interpolated resistance from specific boundary resistance wave
double cvOneDSubdomain::GetBoundResistance(double currentTime){

  if(resistanceTime == NULL || resistanceWave == NULL){
    cout << "ERROR: resistance information is not prescribed"<< endl;
    exit(1);
  }

  // Flow rate is assumed to be periodic
  double cycleTime = resistanceTime[numPressurePts-1];
  double correctedTime = currentTime - static_cast<long>(currentTime / cycleTime) * cycleTime;

  double resistance=0.;
  if(correctedTime >= 0 && correctedTime <= resistanceTime[0]){
    double xi = (correctedTime - resistanceTime[numPressurePts-1]) / (resistanceTime[0] - resistanceTime[numPressurePts-1]);
    resistance = (resistanceWave[numPressurePts-1] + xi * (resistanceWave[0] - resistanceWave[numPressurePts-1]));
  }else{
    int ptr = 0;
    bool wasFound = false;
    while( !wasFound){
      if( correctedTime >= resistanceTime[ptr] && correctedTime <= resistanceTime[ptr+1]){
        wasFound = true;
      }else{
        ptr++;
      }
    }
    // linear interpolation between values
    double xi = (correctedTime - resistanceTime[ptr]) / (resistanceTime[ptr+1] - resistanceTime[ptr]);
    resistance = (resistanceWave[ptr] + xi * (resistanceWave[ptr+1] - resistanceWave[ptr]));
  }
  return resistance;
}

double cvOneDSubdomain::GetBoundAreabyPresWave(double currentTime){
  double pressure =  GetPressure(currentTime);
  return mat->GetArea(pressure, GetLength());
}





/*
* Use the following functions MemIntRCR, MemAdvRCR, ConvPressRCR for RCR BC conditions
* Mem is for "memory" - time integrals
* all added IV 050803
*/

double cvOneDSubdomain::MemIntRCR(double currP, double previousP, double deltaTime, double currentTime){
  double MemIrcr;
  MemIrcr = ConvPressRCR(previousP, deltaTime, currentTime)*expmDtOne(deltaTime)/(alphaRCR*alphaRCR)
    + (currP)/alphaRCR*(deltaTime-expmDtOne(deltaTime)/alphaRCR);
  return MemIrcr;
}

double cvOneDSubdomain::MemAdvRCR(double currP, double previousP, double deltaTime, double currentTime){
  double MemK;
  double Coeff = MemC(currP, previousP, deltaTime, currentTime);
  MemK = Coeff*Coeff/2/alphaRCR*(1-exp(-2*alphaRCR*deltaTime))
    + 2*Coeff*currP/alphaRCR/(proximalResistance + distalResistance)*expmDtOne(deltaTime)
    + deltaTime*pow(currP/(proximalResistance+distalResistance),2);
  return MemK;
}

double cvOneDSubdomain::MemC(double currP, double previousP, double deltaTime, double currentTime){
  double prevTime = currentTime-deltaTime;
  double Cadv;
  Cadv = - ConvPressRCR(previousP, deltaTime, currentTime)/(alphaRCR*proximalResistance*proximalResistance*capacitance)+
    (Q_initial - mat->GetReferencePressure()/proximalResistance)*exp(-alphaRCR*prevTime)
    + currP/(alphaRCR*proximalResistance*proximalResistance*capacitance);
  return Cadv;
}

double cvOneDSubdomain::dMemCdP(void){
  double diffMemCdP;
  diffMemCdP = 1/(alphaRCR*capacitance)/pow(proximalResistance,2);
  return diffMemCdP;
}

double cvOneDSubdomain::ConvPressRCR(double previousP, double deltaTime, double currentTime){
  // Initialization
  if(currentTime <= deltaTime ){MemD = 0.0;}
  // Convoluted pressure
  if(currentTime > deltaTime && currentTime != rcrTime){
    MemD = MemD*exp(-alphaRCR*deltaTime)+previousP*expmDtOne(deltaTime);//all the deltaTime are from previous time step
    rcrTime = currentTime;
  }
  return MemD;
}
// wgyang convolution int(p(t')exp(-alphaRCR(t-t')dt'
double cvOneDSubdomain::MemIntPexp( double previousP, double deltaTime, double currentTime){
  double MemIrcr;
  MemIrcr = ConvPressexp(previousP, deltaTime, currentTime);
  return MemIrcr;
}

/*  Added by Jongmin Seo 04062020 & Hyunjin Kim 09022005
* Adding a new boundary condition for coronary arteries
*/
double cvOneDSubdomain::ConvPressCoronary(double previousP, double deltaTime, double currentTime, double exponent){

  double prevTime=currentTime-deltaTime;
  if (currentTime <= deltaTime) {
    MemD1=0.0;
    MemD2=0.0;
  }

  if (currentTime > deltaTime && currentTime!=corTime) {
    MemD1=MemD1*exp(expo1COR*deltaTime)+expmDtOneCoronary(deltaTime, expo1COR)*(CoefZ1*previousP+CoefY1*getBoundCoronaryValues(prevTime));
    MemD2=MemD2*exp(expo2COR*deltaTime)+expmDtOneCoronary(deltaTime, expo2COR)*(CoefZ2*previousP+CoefY2*getBoundCoronaryValues(prevTime));
    corTime=currentTime;
  }
  if (exponent==expo1COR){ return MemD1; }
  else if (exponent==expo2COR){ return MemD2; }
  cout << "ConvPressCoronary return value undefined!" << std::endl;
  assert(0);
  return 0.0;
}

double cvOneDSubdomain::expmDtOneCoronary(double deltaTime, double exponent){
  double lambda=exponent;
  double expmDt=(exp(lambda*deltaTime)-1)/lambda;
  return expmDt;
}

double cvOneDSubdomain::MemIntCoronary(double currP, double previousP, double deltaTime, double currentTime, double exponent){
  double MemIcor = 0.0;
  double lambda = exponent;
  if (lambda==expo1COR) {
    MemIcor = ConvPressCoronary(previousP, deltaTime, currentTime, lambda)*expmDtOneCoronary(deltaTime,lambda)
    + (CoefZ1*currP+CoefY1*getBoundCoronaryValues(currentTime))*(-deltaTime+expmDtOneCoronary(deltaTime,lambda))/lambda;
  }
  else if (lambda==expo2COR) {
    MemIcor = ConvPressCoronary(previousP, deltaTime, currentTime, lambda)*expmDtOneCoronary(deltaTime,lambda)
    + (CoefZ2*currP+CoefY2*getBoundCoronaryValues(currentTime))*(-deltaTime+expmDtOneCoronary(deltaTime,lambda))/lambda;
  }
  return MemIcor;
}

double cvOneDSubdomain::CORic1(void){
  double CORic1=1/detCOR*(q2COR*(dQ_dT_initial-expo2COR*Q_initial)
       -(p1COR+p2COR*expo1COR)*GetPressure(0)-p2COR*GetdPressuredt(0)
       +b1COR*getBoundCoronaryValues(0));
  return CORic1;
}

double cvOneDSubdomain::CORic2(void){
  double CORic2=1/detCOR*(q2COR*(dQ_dT_initial-expo1COR*Q_initial)
       -(p1COR+p2COR*expo2COR)*GetPressure(0)-p2COR*GetdPressuredt(0)
       +b1COR*getBoundCoronaryValues(0));
  return CORic2;
}

double cvOneDSubdomain::MemCoronary1(double currP, double previousP, double deltaTime, double currentTime){
  double prevTime = currentTime-deltaTime;
  double Cadv;
  Cadv = ConvPressCoronary(previousP, deltaTime, currentTime, expo1COR)+CORic1()*exp(expo1COR*prevTime)
    + (currP*CoefZ1+CoefY1*getBoundCoronaryValues(currentTime))/expo1COR;
  return Cadv;
}

double cvOneDSubdomain::MemCoronary2(double currP, double previousP, double deltaTime, double currentTime){
  double prevTime = currentTime-deltaTime;
  double Cadv;
  Cadv = ConvPressCoronary(previousP, deltaTime, currentTime, expo2COR)+CORic2()*exp(expo2COR*prevTime)
    + (currP*CoefZ2+CoefY2*getBoundCoronaryValues(currentTime))/expo2COR;
  return Cadv;
}

double cvOneDSubdomain::dMemCoronary1dP(void){
  double diffMemCdP;
  diffMemCdP = CoefZ1/expo1COR;
  return diffMemCdP;
}

double cvOneDSubdomain::dMemCoronary2dP(void){
  double diffMemCdP;
  diffMemCdP = CoefZ2/expo2COR;
  return diffMemCdP;
}

double cvOneDSubdomain::dMemIntCoronarydP(double deltaTime, double exponent){
  double dMemIdP = 0.0;
  double lambda=exponent;
  if (lambda==expo1COR) { dMemIdP = CoefZ1*(-deltaTime+expmDtOneCoronary(deltaTime,lambda))/lambda; }
  else if (lambda==expo2COR) { dMemIdP = CoefZ2*(-deltaTime+expmDtOneCoronary(deltaTime,lambda))/lambda; }
  return dMemIdP;
}

double cvOneDSubdomain::MemAdvCoronary(double currP, double previousP, double deltaTime, double currentTime){
  double MemK;
  double Coeff1 = MemCoronary1(currP, previousP, deltaTime, currentTime);
  double Coeff2 = MemCoronary2(currP, previousP, deltaTime, currentTime);
  double currPlv = getBoundCoronaryValues(currentTime);
  double factor = (p0COR*currP-b0COR*currPlv)/q0COR;
  MemK = Coeff1*Coeff1/2/expo1COR*(exp(2*expo1COR*deltaTime)-1)
    - 2*Coeff1*Coeff2*expmDtOneCoronary(deltaTime, expo1COR+expo2COR)
    + Coeff2*Coeff2/2/expo2COR*(exp(2*expo2COR*deltaTime)-1)
    +pow(factor,2)*deltaTime
    +2*Coeff1*factor*expmDtOneCoronary(deltaTime, expo1COR)
    -2*Coeff2*factor*expmDtOneCoronary(deltaTime, expo2COR);
  return MemK;
}

double cvOneDSubdomain::dMemAdvCoronarydP(double currP, double prevP, double deltaTime, double currentTime){
  double Coeff1 = MemCoronary1(currP, prevP, deltaTime, currentTime);
  double Coeff2 = MemCoronary2(currP, prevP, deltaTime, currentTime);
  double currPlv = getBoundCoronaryValues(currentTime);
  double factor = (p0COR*currP-b0COR*currPlv)/q0COR;
  double dMemKdP;
  dMemKdP = (Coeff1/expo1COR*(exp(2*expo1COR*deltaTime)-1)
    -2*Coeff2*expmDtOneCoronary(deltaTime, expo1COR+expo2COR)
    +2*factor*expmDtOneCoronary(deltaTime, expo1COR))*dMemCoronary1dP()
    +(-2*Coeff1*expmDtOneCoronary(deltaTime, expo1COR+expo2COR)
    +Coeff2/expo2COR*(exp(2*expo2COR*deltaTime)-1)
    -2*factor*expmDtOneCoronary(deltaTime, expo2COR))*dMemCoronary2dP()
    +(2*factor*deltaTime+2*Coeff1*expmDtOneCoronary(deltaTime, expo1COR)
    -2*Coeff2*expmDtOneCoronary(deltaTime, expo2COR))*p0COR/q0COR;
  return dMemKdP;
}

double cvOneDSubdomain::ConvPressexp(double previousP, double deltaTime, double currentTime){
  // Initialization
  if(currentTime <= deltaTime ){MemConvP = 0.0;}
  // Convoluted pressure
  if(currentTime > deltaTime && currentTime != rcrTime2){
    MemConvP = MemConvP*exp(-alphaRCR*deltaTime)+previousP/alphaRCR*expmDtOne(deltaTime);//all the deltaTime are from previous time step
    rcrTime2 = currentTime;
  }
  return MemConvP;
}
// wgyang convolution int(S(t')^(-3/2)exp(-alphaRCR(t-t')dt'
double cvOneDSubdomain::MemIntSexp( double previousS, double deltaTime, double currentTime){
  double MemIrcr;
  MemIrcr = ConvSexp(previousS, deltaTime, currentTime);
  return MemIrcr;
}

double cvOneDSubdomain::ConvSexp(double previousS, double deltaTime, double currentTime){
  // Initialization
  if(currentTime <= deltaTime ){MemConvS = 0.0;}
  // Convoluted pressure
  if(currentTime > deltaTime && currentTime != rcrTime3){
    //MemConvS = MemConvS*exp(-alphaRCR*deltaTime)+pow(previousS,-1.5)/alphaRCR*expmDtOne(deltaTime);//all the deltaTime are from previous time step  convolution  int S(t)^(-3/2)exp(-alpha*(t-t')dt
    MemConvS=MemConvS*exp(-alphaRCR*deltaTime)+previousS/alphaRCR*expmDtOne(deltaTime);
    rcrTime3 = currentTime;
  }

  return MemConvS;
}

double cvOneDSubdomain::dMemIntRCRdP(double deltaTime){
  double dMemIdP;
  dMemIdP = (deltaTime - expmDtOne(deltaTime)/alphaRCR)/alphaRCR;
  // Cout << dMemIdP ;
  return dMemIdP;
}


double cvOneDSubdomain::dMemAdvRCRdP(double currP, double prevP, double deltaTime, double currentTime){
  double MemCoeff = MemC(currP, prevP, deltaTime, currentTime);
  double dMemKdP;
  dMemKdP = (1-exp(-2*alphaRCR*deltaTime))/alphaRCR*MemCoeff*dMemCdP()
    + expmDtOne(deltaTime)/alphaRCR/(proximalResistance + distalResistance)*(MemCoeff + currP*dMemCdP())
    + 2*deltaTime*currP/pow(proximalResistance + distalResistance,2);
  //cout<<" dMemKdP"<<" "<<dMemKdP<<endl;
  return dMemKdP;
}
