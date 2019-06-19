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
  K = new double[1000];
  for(int i=0;i<1000;i++){
      K[i]=0.0;
  }
  branchAngle = 90.0;
}

cvOneDSubdomain::~cvOneDSubdomain(){
  if(connectivities != NULL)  delete [] connectivities;
  if(nodes != NULL)      delete [] nodes;
  if(pressureWave != NULL)  delete [] pressureWave;
  if(pressureTime != NULL)  delete [] pressureTime;
  if(resistanceWave != NULL)  delete [] resistanceWave;
  if(resistanceTime != NULL)  delete [] resistanceTime;
  if(presslv != NULL) delete [] presslv;//kimhj 09032005
  if(PressLVWave != NULL)  delete [] PressLVWave;//kimhj 09032005
  if(PressLVTime != NULL)  delete [] PressLVTime;//kimhj 09032005
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
  printf("subdomain cpp setupMaterial matID=%i  \n", matID);
  mat->SetAreas_and_length(S_initial, S_final, fabs(z_out - z_in));
}

void cvOneDSubdomain::SetBoundValue(double boundV){
  switch(boundType){
    case BoundCondTypeScope::PRESSURE:
      boundValue = mat->GetArea(boundV*1333.2237, fabs(z_out-z_in));
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

// kimhj added for coronary boundary conditions kimhj 09022005
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
      pressureWave[i] = pres[i]*1333.2237;
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

void cvOneDSubdomain::SetBoundRCRValues(double *rcr, int num){//added IV 050803

  rcrTime = 0.0;

  assert(num==3);
  proximalResistance = rcr[0];
  capacitance = rcr[1] ;
  distalResistance= rcr[2];
  alphaRCR = (proximalResistance+distalResistance)/(proximalResistance*distalResistance*capacitance);
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


void cvOneDSubdomain::SaveK(double k, int i){
  K[i]=k;
}



/*
* Use the following functions MemIntRCR, MemAdvRCR, ConvPressRCR for RCR BC conditions
* Mem is for "memory" - time integrals
* all added IV 050803
*/

double cvOneDSubdomain::MemIntRCR(double currP, double previousP, double deltaTime, double currentTime){
  double MemIrcr;
  MemIrcr = ConvPressRCR(previousP, deltaTime, currentTime)*expmDtOne(deltaTime)/(alphaRCR*alphaRCR)
    + currP/alphaRCR*(deltaTime-expmDtOne(deltaTime)/alphaRCR);
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