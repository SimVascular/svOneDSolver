//
//  MthSegmentModel.cxx
//

# include "cvOneDMthSegmentModel.h"
# include "cvOneDGlobal.h"

// stabilization parameter blows up if you have a very small or zero
// kinematic viscosity
#define SMALL_KINEMATIC_VISCOSITY 0.0001

int cvOneDMthSegmentModel::STABILIZATION;

cvOneDMthSegmentModel::cvOneDMthSegmentModel(const vector<cvOneDSubdomain*>& subdList,
                                             const vector<cvOneDFEAJoint*>& jtList,
                                             const vector<int>& outletList, long quadPoints_):
                       cvOneDMthModelBase(subdList, jtList, outletList), quadrature_(quadPoints_){
  quadPoints = quadPoints_;
  weight = new double[quadPoints];
  xi = new double[quadPoints];
}

cvOneDMthSegmentModel::~cvOneDMthSegmentModel(){
  delete [] xi;
  delete [] weight;
}

void cvOneDMthSegmentModel::SetEquationNumbers(long element, cvOneDDenseMatrix* elementMatrix, int ithSubdomain){
  long eqNumbers[4];
  GetEquationNumbers(element, eqNumbers, ithSubdomain);
  elementMatrix->SetEquationNumbers(eqNumbers);
}

void cvOneDMthSegmentModel::FormNewtonLHS(cvOneDFEAMatrix* lhsMatrix){
  int i;
  lhsMatrix->Clear();

  cvOneDDenseMatrix elementMatrix(4, "eLhsMatrix");
  // no boundary related terms
  for(i = 0; i < subdomainList.size(); i++){
    for(long element = 0; element < subdomainList[i]->GetNumberOfElements();element++){
      FormElementLHS(element, &elementMatrix, i);
      lhsMatrix->Add(elementMatrix);
      // cout<< elementMatrix << endl;
    }
  }
}

void cvOneDMthSegmentModel::FormNewtonRHS(cvOneDFEAVector* rhsVector){
  rhsVector->Clear();
  cvOneDFEAVector elementVector(4, "eRhsVector");
  for(int i = 0; i < subdomainList.size(); i++){
    for(long element = 0; element<subdomainList[i]->GetNumberOfElements(); element++){
      FormElementRHS( element, &elementVector, i);
      rhsVector->Add( elementVector);
    }
  }
}

double cot(double x){
  return cos(x)/sin(x);
}

double cvOneDMthSegmentModel::N_MinorLoss(long ith){

  char propName[256];

  cvOneDSubdomain *sub = subdomainList[ith];

  MinorLoss minorLoss=sub->GetMinorLossType();

  // localize the values of the current approximation for U on this element
  long eqNumbers[4];
  double S[3];    // area nodal values    (S = U[0])
  double Q[3];    // flow rate nodal values  (Q = U[1])

   // first node of this segment(1), the stenosis (minor loss) segment,not used for converging flow
  GetEquationNumbers(0, eqNumbers, ith);
  S[1] = currSolution->Get( eqNumbers[0]);
  Q[1] = currSolution->Get( eqNumbers[1]);

  if(minorLoss == MinorLossScope::NONE ){//|| minorLoss != MinorLossScope::STENOSIS ){
    return sub->GetMaterial()->GetN(S[1]);
  }

  //   In case of branch, previous seg might not be adjacent seg
  // Subdomain *adjacent =subdomainList[ith-1];
  cvOneDSubdomain *upstream = subdomainList[sub->GetUpstreamSeg()];
  cvOneDSubdomain *branch = subdomainList[sub->GetBranchSeg()];

  double L = sub->GetLength();
  double Kv, Kt, Re0, La, D[2], func, N, theta, min;
  strcpy(propName,"kinematic viscosity");
  double kinViscosity0 = upstream->GetMaterial()->GetProperty(propName);

  min=0.6; // if a for branches is too small, troubles abound

  // The last node of segment (3) at the upstream of the stenosis (minor loss) segment // not true for converging flow
  GetEquationNumbers(upstream->GetNumberOfElements()-1, eqNumbers, sub->GetUpstreamSeg());
  S[0] = currSolution->Get( eqNumbers[2]);
  Q[0] = currSolution->Get( eqNumbers[3]);

  // first node of the branch (minor loss) segment, only used for through(2) calcs
   GetEquationNumbers(0, eqNumbers, sub->GetBranchSeg());
  S[2] = currSolution->Get( eqNumbers[0]);
  Q[2] = currSolution->Get( eqNumbers[1]);

  // ratio of flow rate between branch and combined flow
  double q = Q[1]/Q[0];
  // ratio of area between branch and combined flow
  double a = S[1]/S[0];

  // why are we doing this?  Sometimes flow is negative.
  if(Q[0] < 0.000001 || Q[1] < 0.000001){
    return sub->GetMaterial()->GetN(S[1]);
  }

  switch(minorLoss){
    case MinorLossScope::STENOSIS:
      Kt = 1.52;
      // modify to lessen stenosis
      D[0] = sqrt(4*S[0]/3.14159265358979323846);
      D[1] = sqrt(4*S[1]/3.14159265358979323846);
      La = 0.83*L + 1.64*D[1];
      // La = L;
      Kv = 32. * La / D[0] * 1.0/a*1.0/a;
      Re0 = D[0] * Q[0] / (S[0] * kinViscosity0);
      func= Kv/Re0 + Kt / 2. * (1.0/a - 1.) * (1.0/a - 1.);
      // func= Kv/Re0 + Kt / 2. * (1/a - 1.) * (1/a - 1.)*abs(Q[0])/Q[0]+Ku*1.06*L*sub->getDQDt()/Q;  // for non-steady flow
      func = 2.0 * func*a*a/q/q;   // a^2/q^2 switch from upstream segment to stenosed segment
      break;

    case MinorLossScope::BRANCH_THROUGH_DIVIDING:
      theta = sub -> GetBranchAngle();

      q = Q[2]/Q[0]; // branch(1)/combined(0)
      a = S[2]/S[0];

      // having problems with a too small, try hokey fix for now 02-08-02
      if(a<min){
        // a=min;
        // return sub->GetMaterial()->GetN(S[1]);
      }
      // Wood K32,  modified, times combined flow div by leg flow
      func = (-0.03*(1-q)*(1-q) - 0.35 * q*q + 0.2*q*(1-q));  // no difference
      // modified coeff -- too strong
      // func=func * Q[0]*Q[0]/(Q[1]*Q[1]);
      // cout<<"through dividing   q "<<q<< " func "<<func<<endl;
      break;

    // problem with this case.
    case MinorLossScope::BRANCH_SIDE_DIVIDING:
      theta = sub -> GetBranchAngle();

      // having problems with a too small, try hokey fix for now 02-08-02
      if(a<min){
      // a=min;
      // return sub->GetMaterial()->GetN(S[1]);
      }
      // q= Q[1]/Q[0]; // branch(this=1)/combined(0)

      // Wood K31
      func = (-0.95*(1-q)*(1-q) + q*q*(1.3*cot((3.14159265358979323846-theta)/2) - 0.3 + (.4 - .1*a)/(a*a)) +
             -0.4*q*(1-q)*(1+(1/a))*cot((3.14159265358979323846 - theta)/2));// one sign change
      // func = 0.95*(1-q)*(1-q) + q*q*(1.3*cot((3.14159265358979323846-theta)/2) - 0.3 + (.4 - .1*a)/(a*a)) +
      //        0.4*q*(1-q)*(1+(1/a))*cot((3.14159265358979323846 - theta)/2);

      // modified coefficient -- too strong
      // func=func/(q*q);
      // cout<<"side dividing   q "<<q <<" func "<<func<< endl;
      break;

    case MinorLossScope::BRANCH_THROUGH_CONVERGING:

      // The segment at the upstream of the stenosis (minor loss) segment // not true for converging flow
      GetEquationNumbers(0, eqNumbers, sub->GetUpstreamSeg());
      S[0] = currSolution->Get( eqNumbers[0]);
      Q[0] = currSolution->Get( eqNumbers[1]);

      // the stenosis (minor loss) segment // not true for converging flow
      GetEquationNumbers(subdomainList[ith]->GetNumberOfElements()-1, eqNumbers, ith);
      S[1] = currSolution->Get( eqNumbers[2]);
      Q[1] = currSolution->Get( eqNumbers[3]);

      theta = sub -> GetBranchAngle();

      q = Q[2]/Q[0]; // branch/combined
      a = S[2]/S[0];
      // having problems with a too small, try hokey fix for now 02-08-02
      if(a<min){
        a = min;
        return sub->GetMaterial()->GetN(S[1]);
      }

      // this method from Marc Serre  Using upstream area instead of downstreram area in divisor a
      // it also wants branch flow instead of through flow, so modify Q for 90 degree
      // func = 2*(1-0.2)*q-q*q+.04*q*q*(1/sqrt(a*a*a));

      // Serre for arbitrary angle
      //func = (1.61*q+(-1+ .04*(1/sqrt(a*a*a)) - 1.74*(1/a)*cos(theta))*q*q);
      //cout<<"Serre method conv-through "<<func<<endl;

      // Wood K23
      func = (-0.03*(1+q)*(1+q) + q*q*(1 + 1.62*(cos(theta)/a - 1) - 0.38*(1-a))+(2-a)*q*(1+q)); // sign change in front of last term and in last term
      // func =(0.03*(1-q)*(1-q) - q*q*(1 + 1.62*(cos(theta)/a - 1) - 0.38*(1-a))+(2-a)*q*(1-q));
      // modified value, divide by flow ratio squared. I found that this was too strong.
      // func= func * Q[0]*Q[0]/(Q[1]*Q[1]);
      // cout <<" branch converging q:"<<q<<"func:"<< func ;
      break;

    case MinorLossScope::BRANCH_SIDE_CONVERGING:
      // The segment of combined flow (downstream), want first segment
      GetEquationNumbers(0, eqNumbers, sub->GetUpstreamSeg());
      S[0] = currSolution->Get(eqNumbers[0]);
      Q[0] = currSolution->Get(eqNumbers[1]);

      // the branch segment. want last element
      GetEquationNumbers(subdomainList[ith]->GetNumberOfElements()-1, eqNumbers, ith);
      S[1] = currSolution->Get(eqNumbers[2]);
      Q[1] = currSolution->Get(eqNumbers[3]);

      q = Q[1]/Q[0]; // branch(1)/combined(3)
      a = S[1]/S[0];
      // having problems with a too small, try hokey fix for now 02-08-02
      if(a < min){
        a = min;
        return sub->GetMaterial()->GetN(S[1]);
      }
      theta = sub -> GetBranchAngle();
      // Marc Serre eq29 q=branch/upstream
      // func=(1-1.8*a)*((q*q)/(a*a) - 1);
      // cout<< "Serre method - conv-branch "<< func <<endl;
      // Wood K13 // try changing first and last terms to (1+q)
      func =(0.92*(1+q)*(1+q)+(q*q)*(1.2*(cos(theta)/a-1)+0.8*(1-1/(a*a))-(1-a)*cos(theta)/a) +(2-a)*q*(1+q) );
      //  func =-0.92*(1-q)*(1-q)-(q*q)*(1.2*(cos(theta)/a-1)+0.8*(1-1/(a*a))-(1-a)*cos(theta)/a) +(2-a)*q*(1-q) ;
      // use modified value..
      //func=func/(q*q);
      //  cout<<"side converging  " <<func<< " q "<<q ;
      break;
    case MinorLossScope::BIFURCATION_BRANCH:

      theta = sub ->GetBranchAngle();

      /* JNIEnv* env=BFSolver::jenv;
      jobject obj=BFSolver::jobj;// OneDModel object
      // cout<< "in bif-branch " << endl;
      jclass cls = env->GetObjectClass(obj);
      // cout<< "jclass " << cls << endl;

      // jmethodID mid = env->GetMethodID(cls, "callback", "(DDD)D");
      jmethodID mid = env->GetMethodID(cls, "findBifurcationK", "(DDD)D");
      if (mid == 0) {
        cout<< "method id is zero " << endl;
        return sub->GetMaterial()->GetProperty( "N");
      }
      jdouble ans = env->CallDoubleMethod(obj, mid, theta,Q[0],S[0]);
      // cout<< ">>found method id ans = "<<ans  << endl;
      func = ans;
      */
      //this is turned off because we don't use java in the debug version
      strcpy(propName,"N");
      return sub->GetMaterial()->GetProperty(propName);
      break;

    /*
    case MinorLossScope::BEND_90:
      func = .45;
      break;
    case MinorLossScope::BEND_45:
      func = .2;
      break;
    case MinorLossScope::BEND_180: // don't think I can really model this, case
      func = 0.4;
      break;
    default:
      return sub->GetMaterial()->GetProperty( "N");
    */
  }// end of switch

  // integrate (1) over z to get formula for dP
  // jing's
  // N = func * Q[0] * Q[0] * S[1] * S[1] / (S[0] * S[0] * Q[1] * L);
  // mine.. simple

  N = func *  Q[1]  / (2 * L);

  sub->SaveK(N,(int)(cvOneDBFSolver::currentTime / cvOneDBFSolver::deltaTime) );

  // don't want to have less than the default Puoseille
  strcpy(propName,"N");
  double std= sub->GetMaterial()->GetProperty(propName);
  if( -N > std){
    return std;
  }

  // cout << " -N " << -N << " Q(1) :"<< Q[1] << endl;
  return -N;

}

void cvOneDMthSegmentModel::FormElementLHS(long element, cvOneDDenseMatrix* elementMatrix, long ith){

  char propName[256];
  //Framework
  cvOneDSubdomain *sub = subdomainList[ith];
  cvOneDMaterial* material = sub->GetMaterial();
  strcpy(propName,"density");
  double density = material->GetProperty(propName);
  strcpy(propName,"delta");
  double delta = material->GetProperty(propName);
  strcpy(propName,"kinematic viscosity");
  double kinViscosity = material->GetProperty(propName);
  //  material->GetProperty( "N");
  double N = N_MinorLoss(ith);

  double k11,k12,k21,k22;
  BoundCondType bound = sub->GetBoundCondition();

  // now get the element information from the domain
  const cvOneDFiniteElement* finiteElement = sub->GetElement(element);

  // element nodes
  int numberOfNodes = 2;
  double nodes[2];
  sub->GetNodes( element, nodes);

  double shape[2];
  double DxShape[2];
  double jacobian;

  // get the weight and quadrature point
  quadrature_.Get( weight, xi);

  // localize the values of the current approximation for U on this element
  long eqNumbers[4];
  GetEquationNumbers(element, eqNumbers, ith);
  double S[2];    // area nodal values    (S = U[0])
  double Q[2];    // flow rate nodal values  (Q = U[1])
  S[0] = currSolution->Get(eqNumbers[0]);
  Q[0] = currSolution->Get(eqNumbers[1]);
  S[1] = currSolution->Get(eqNumbers[2]);
  Q[1] = currSolution->Get(eqNumbers[3]);

  if(cvOneDGlobal::debugMode){
    printf("(Debug) Assembling Element %ld\n",element);  
    printf("--- End node solutions\n");  
    printf("S[0]: %e\n",S[0]);
    printf("Q[0]: %e\n",Q[0]);
    printf("S[1]: %e\n",S[1]);
    printf("Q[1]: %e\n",Q[1]);
    //getchar();
  }

  // set the equation numbers for the element matrix
  elementMatrix->SetEquationNumbers( eqNumbers);
  elementMatrix->Clear();

  for( int l = 0; l < quadPoints; l++){

    finiteElement->Evaluate( xi[l], shape, DxShape, &jacobian);

    // evaluate U at the current integration point
    // U[0] = S, U[1] = Q
    double U[2];
    U[0] = finiteElement->Interpolate( xi[l], S);
    U[1] = finiteElement->Interpolate( xi[l], Q);

    // get the coordinate corresponding to the quadrature point
    // the elements are assumed to be isoparametric
    double z = finiteElement->Interpolate( xi[l], nodes);

    // more values coming from constitutive equations
    double pressure = material->GetPressure( U[0], z);
    double Outflow = material->GetOutflowFunction( pressure, z);//0.0
    double DpDS = material->GetDpDS( U[0], z);
    double DpDz = material->GetDpDz( U[0], z);
    double DOutflowDp = material->GetDOutflowDp( pressure, z);// 0.0
    double IntegralpD2S = 0.0;

    if(cvOneDGlobal::debugMode){
      printf("(Debug) Assembling Element %ld\n",element);  
      printf("--- Element Values\n");  
      printf("z: %e\n",z);
      printf("pressure: %e\n",pressure);
      printf("Outflow: %e\n",Outflow);
      printf("DpDS: %e\n",DpDS);
      printf("DpDz: %e\n",DpDz);
      printf("DOutflowDp: %e\n",DOutflowDp);
      printf("IntegralpD2S: %e\n",IntegralpD2S);
      //getchar();
    }


    if(cvOneDGlobal::CONSERVATION_FORM == 1) {
      IntegralpD2S = material->GetIntegralpD2S( U[0], z);//IV added IntegralpD2S 01-18-03
    }

    // evaluate the matrices A, C and K at U  CU=G
    // recall that:
    //    A11 = 0.0  A12 = 1.0  /partdF/partdU
    //    C12 = 0.0
    //    K11 = 0.0  K12 = 0.0  K21 = 0.0
    double aux = U[1]/U[0];        // aux = Q/S
    double A12 = 1.0;
    double A21 = -(1.0+delta)*aux*aux+U[0]/density*DpDS;
    double A22 = 2.0*(1.0+delta)*aux;
    double C11 = -Outflow / U[0];
    double C21 = -1.0/density*DpDz;
    double C22 = N / U[0];
    double K22 = kinViscosity;

    if(cvOneDGlobal::debugMode){
      printf("(Debug) Assembling Element %ld\n",element);  
      printf("--- Element Values 2\n");  
      printf("aux: %e\n",aux);
      printf("A12: %e\n",A12);
      printf("A21: %e\n",A21);
      printf("A22: %e\n",A22);
      printf("C11: %e\n",C11);
      printf("C21: %e\n",C21);
      printf("C22: %e\n",C22);
      printf("K22: %e\n",K22);
      //getchar();
    }

    // K22=0.0;

    // 01-18-03
    // in IV's formulation we use F and the corresponding C
    // which is different than the previous C, so we will call it CF
    // the contribution of the F-term into the tangent matrix is different
    // CF12=0.0
    double CF11 = -Outflow / U[0];
    // double CF21 = -1.0/density/U[0]*IntegralpD2S;
    double CF21 = 1.0/density/U[0]*IntegralpD2S;//sign error corrected by IV 06-15-04
    double CF22 = N / U[0];

    if(cvOneDGlobal::debugMode){
      printf("(Debug) Assembling Element %ld\n",element);  
      printf("--- Element Values 2\n");  
      printf("CF11: %e\n",CF11);
      printf("CF21: %e\n",CF21);
      printf("CF22: %e\n",CF22);
      //getchar();
    }

    // evaluate the tau matrix
    double tau[4];
    double h = nodes[1] - nodes[0];
    double A[4] = { 0.0, 1.0, A21, A22};
    double C[4] = { C11, 0.0, C21, C22};
    double modA[4];
    double modC[4];

    if(cvOneDGlobal::debugMode){
      printf("(Debug) Assembling Element %ld\n",element);  
      printf("C[0]: %e\n",C[0]);
      printf("C[1]: %e\n",C[1]);
      printf("C[2]: %e\n",C[2]);
      printf("C[3]: %e\n",C[3]);
      getchar();
    }

    if(STABILIZATION == 1){

      GetModulus(A, modA);

      // stabilization parameter blows up if you have a very small or zero
      // kinematic viscosity
      if(kinViscosity < SMALL_KINEMATIC_VISCOSITY){
        modC[0] = 0;
        modC[1] = 0;
        modC[2] = 0;
        modC[3] = 0;
      }else{
        GetModulus(C,modC);
      }

      // this is actually the inverse of tau...
      tau[0] = (2.0/deltaTime) + 2.0/h*modA[0] + modC[0];
      tau[1] = 2.0/h*modA[1] + modC[1];
      tau[2] = 2.0/h*modA[2] + modC[2];
      tau[3] = (2.0/deltaTime) + 2.0/h*modA[3] + 12.0/(h*h)*K22 + modC[3];

      if(cvOneDGlobal::debugMode){
        printf("(Debug) Assembling Element %ld\n",element);  
        printf("kinViscosity: %e\n",kinViscosity);
        printf("SMALL_KINEMATIC_VISCOSITY: %e\n",SMALL_KINEMATIC_VISCOSITY);
        printf("h: %e\n",h);
        printf("modA[0]: %e\n",modA[0]);
        printf("modC[0]: %e\n",modC[0]);
        printf("modA[1]: %e\n",modA[1]);
        printf("modC[1]: %e\n",modC[2]);
        printf("modA[2]: %e\n",modA[2]);
        printf("modA[3]: %e\n",modA[3]);
        printf("modC[3]: %e\n",modC[3]);
        getchar();
      }

      double det = tau[0]*tau[3]-tau[1]*tau[2];

      // correcting tau...
      double temp = tau[0];
      tau[0] =  tau[3]/det;
      tau[1] = -tau[1]/det;
      tau[2] = -tau[2]/det;
      tau[3] =  temp/det;

      if(cvOneDGlobal::debugMode){
        printf("(Debug) Assembling Element %ld\n",element);  
        printf("tau[0]: %e\n",tau[0]);
        printf("tau[1]: %e\n",tau[1]);
        printf("tau[2]: %e\n",tau[2]);
        printf("tau[3]: %e\n",tau[3]);
        getchar();
      }

    } // end STABILIZATION

    for( int a = 0; a < numberOfNodes; a++){

      for( int b = 0; b < numberOfNodes; b++){

        double jw = jacobian * weight[l];

        // DG terms
        if(cvOneDGlobal::CONSERVATION_FORM == 1){
          // IV's formulation 01-18-03
          k11 = deltaTime*(shape[a]*CF11*shape[b])-shape[a]*shape[b];
          k12 = deltaTime*(A12*shape[b]*DxShape[a]);
          k21 = deltaTime*(DxShape[a]*A21*shape[b]+shape[a]*CF21*shape[b]);
          k22 = deltaTime*(DxShape[a]*A22*shape[b]-DxShape[a]*K22*DxShape[b]+shape[a]*CF22*shape[b]) - shape[a]*shape[b];
        }else{
          // Here is Brooke's version that I am not using IV 01-30-03
          k11 = deltaTime*(-shape[a]*C11*shape[b])+shape[a]*shape[b];
          k12 = deltaTime*(shape[a]*DxShape[b]);
          k21 = deltaTime*(shape[a]*A21*DxShape[b]-shape[a]*C21*shape[b]);
          k22 = deltaTime*(shape[a]*A22*DxShape[b]+DxShape[a]*K22*DxShape[b]-shape[a]*C22*shape[b]) + shape[a]*shape[b];
        }

        if(STABILIZATION == 1){
          // GLS terms
          // Create an auxiliary matrix to handle some of the terms
          double auxa[4];
          auxa[0] = -shape[a]*C11; // A11 = 0.0
          auxa[1] = DxShape[a]; // A12 = 1.0, C12 = 0.0
          auxa[2] = DxShape[a]*A21-shape[a]*C21;
          auxa[3] = DxShape[a]*A22-shape[a]*C22;

          double auxb[4];
          auxb[0] = -shape[b]*C11; // A11 = 0.0
          auxb[1] = DxShape[b]; // A12 = 1.0, C12 = 0.0
          auxb[2] = DxShape[b]*A21-shape[b]*C21;
          auxb[3] = DxShape[b]*A22-shape[b]*C22;

          if(cvOneDGlobal::debugMode){
            printf("(Debug) Assembling Element %ld\n",element);  
            printf("auxa[0]: %e\n",auxa[0]);
            printf("auxa[1]: %e\n",auxa[1]);
            printf("auxa[2]: %e\n",auxa[2]);
            printf("auxa[3]: %e\n",auxa[3]);
            printf("auxb[0]: %e\n",auxb[0]);
            printf("auxb[1]: %e\n",auxb[1]);
            printf("auxb[2]: %e\n",auxb[2]);
            printf("auxb[3]: %e\n",auxb[3]);
            getchar();
          }

          //
          // Multiply the matrices to obtain the GLS contribution into auxa
          //
          // Contains the product: auxa * tau
          double auxc[4];
          auxc[0] = auxa[0]*tau[0]+auxa[1]*tau[2];
          auxc[1] = auxa[0]*tau[1]+auxa[1]*tau[3];
          auxc[2] = auxa[2]*tau[0]+auxa[3]*tau[2];
          auxc[3] = auxa[2]*tau[1]+auxa[3]*tau[3];

          // auxa will contain the product: auxa * tau * auxb
          auxa[0] = auxc[0]*auxb[0]+auxc[1]*auxb[2];
          auxa[1] = auxc[0]*auxb[1]+auxc[1]*auxb[3];
          auxa[2] = auxc[2]*auxb[0]+auxc[3]*auxb[2];
          auxa[3] = auxc[2]*auxb[1]+auxc[3]*auxb[3];

          if(cvOneDGlobal::debugMode){
            printf("(Debug) Assembling Element %ld\n",element);  
            printf("auxa[0]: %e\n",auxa[0]);
            printf("auxa[1]: %e\n",auxa[1]);
            printf("auxa[2]: %e\n",auxa[2]);
            printf("auxa[3]: %e\n",auxa[3]);
            printf("auxc[0]: %e\n",auxc[0]);
            printf("auxc[1]: %e\n",auxc[1]);
            printf("auxc[2]: %e\n",auxc[2]);
            printf("auxc[3]: %e\n",auxc[3]);
            printf("deltaTime: %e\n",deltaTime);
            getchar();
          }


          //
          // now sum the GLS terms to the DG terms
          //
          k11 += deltaTime*auxa[0];
          k12 += deltaTime*auxa[1];
          k21 += deltaTime*auxa[2];
          k22 += deltaTime*auxa[3];

        } // end stabilization

        if(cvOneDGlobal::debugMode){
          printf("(Debug) Assembling Element %ld\n",element);  
          printf("--- Element Values 2\n");  
          printf("k11: %e\n",k11);
          printf("k12: %e\n",k12);
          printf("k21: %e\n",k21);
          printf("k22: %e\n",k22);
          getchar();
        }

        elementMatrix->Add( 2*a  , 2*b  , k11*jw);
        elementMatrix->Add( 2*a  , 2*b+1, k12*jw);
        elementMatrix->Add( 2*a+1, 2*b  , k21*jw);
        elementMatrix->Add( 2*a+1, 2*b+1, k22*jw);

      }
    }
  }

  if(cvOneDGlobal::CONSERVATION_FORM){

    //Inlet flux term (at z=z_inlet) which is the linearized F-KU IV 01-28-03
    if (element == 0){
      long node = 0;
      double z = sub->GetNodalCoordinate( node);
      double aux = Q[0]/S[0];
      double DpDS = material->GetDpDS( S[0], z);
      finiteElement->Evaluate( z, shape, DxShape, &jacobian);

      int a = 0;
      for( int b = 0; b < numberOfNodes; b++){
        //(1-b)= trick because shape is note defined in z coord;
        //has to be changed if other shape functions are used IV 02-07-03
        //double x=(1.0-(double)b);
        double Inlet11 = 0;
        double Inlet12 = (1.0-(double)b);
        double Inlet21 = (1.0-(double)b)*(-(1.0+ delta)*aux*aux + S[0]/density*DpDS);
        //double Inlet22 = 2*(1.0+ delta)*aux*(1.0-(double)b)- kinViscosity*DxShape[b];
        double Inlet22 = 2*(1.0+ delta)*aux*(1.0-(double)b);//without viscosity in flux term

        elementMatrix->Add( 2*a  , 2*b  , Inlet11*deltaTime);
        elementMatrix->Add( 2*a  , 2*b+1, Inlet12*deltaTime);
        elementMatrix->Add( 2*a+1, 2*b  , Inlet21*deltaTime);
        elementMatrix->Add( 2*a+1, 2*b+1, Inlet22*deltaTime);
      }
    }

    if(bound==BoundCondTypeScope::NOBOUND||bound==BoundCondTypeScope::PRESSURE
      ||bound==BoundCondTypeScope::PRESSURE_WAVE||bound==BoundCondTypeScope::AREA||bound==BoundCondTypeScope::FLOW){
      //Outlet flux term (at z=z_outlet) which is the linearized F-KU IV 01-28-03
      if (element == (sub->GetNumberOfElements())-1){
        //cout<<"elementLHS "<<element<<endl;
        double z = sub->GetOutletZ();//checked IV 02-03-03
        finiteElement->Evaluate( z, shape, DxShape, &jacobian);
        double Pressure= material->GetPressure( S[1], z);
        double aux= Q[1]/S[1];
        double DpDS = material->GetDpDS( S[1], z);
        int a = 1;

        for( int b = 0; b < numberOfNodes; b++){
          //b= trick because shape is note defined in z coord;
          //has to be changed if other shape functions are used IV 02-07-03
          //double x=double b;
          double Outlet11 = 0.0;
          double Outlet12 = (double)b;
          double Outlet21 = (double)b*(-(1.0+ delta)*aux*aux + S[1]/density*DpDS);
          //double Outlet22 = 2*(1.0+ delta)*aux*(double)b - kinViscosity*DxShape[b];
          //cout<<Outlet21<<"-"<<Outlet22<<" ";
          double Outlet22 = 2*(1.0+ delta)*aux*(double)b ;//without viscosity in flux
          //*/
          //try no adv term
          /*
          double Outlet21 = (double)b*(S[1]/density*DpDS);
          double Outlet22 = 0.0;
          //*/
          //try no M2h2
          /*
          double Outlet21 = 0.0;
          double Outlet22 = 0.0;
          //*/
          //try linear downstream domain-Hughes
          /* double Cp = material->GetLinCompliance(z);
             double Outlet21 = (double)b*S[1]/density/Cp;//linear downstream domain-Hughes
             double Outlet22 = 0.0;//linear downstream domain hughes
          */
          elementMatrix->Add( 2*a  , 2*b  , -Outlet11*deltaTime);
          elementMatrix->Add( 2*a  , 2*b+1, -Outlet12*deltaTime);
          elementMatrix->Add( 2*a+1, 2*b  , -Outlet21*deltaTime);
          elementMatrix->Add( 2*a+1, 2*b+1, -Outlet22*deltaTime);
        }//end for
      }//end outlet term
    }//end for compute full flux if no outlet BC or Dirichlet BC
  }//end if(CONSERVATION_FORM)
}//end function elementLHS

void cvOneDMthSegmentModel::FormElementRHS(long element, cvOneDFEAVector* elementVector, long ith){

  char propName[256];

  cvOneDSubdomain *sub = subdomainList[ith];
  // get the material properties
  cvOneDMaterial* material = sub->GetMaterial();
  strcpy(propName,"density");
  double density = material->GetProperty(propName);
  strcpy(propName,"delta");
  double delta = material->GetProperty(propName);
  strcpy(propName,"kinematic viscosity");
  double kinViscosity = material->GetProperty(propName);
  //  material->GetProperty( "N");
  double N = N_MinorLoss(ith);

  BoundCondType bound = sub->GetBoundCondition();
  // now get the element information from the domain
  const cvOneDFiniteElement* finiteElement=sub->GetElement(element);

  int numberOfNodes = 2;
  double nodes[2];
  sub->GetNodes( element, nodes);

  // shape functions and derivatives
  double shape[2];
  double DxShape[2];
  double jacobian;

  // get the weight and quadrature point
  quadrature_.Get(weight, xi);

  // localize the components of the solution at t_{n} on this element
  long eqNumbers[4];
  GetEquationNumbers(element, eqNumbers, ith);
  double Sn[2];
  double Qn[2];
  Sn[0] = prevSolution->Get(eqNumbers[0]);
  Qn[0] = prevSolution->Get(eqNumbers[1]);
  Sn[1] = prevSolution->Get(eqNumbers[2]);
  Qn[1] = prevSolution->Get(eqNumbers[3]);

  // localize the components of the current approximation for U on this element
  double S[2];
  double Q[2];
  S[0] = currSolution->Get( eqNumbers[0]);
  Q[0] = currSolution->Get( eqNumbers[1]);
  S[1] = currSolution->Get( eqNumbers[2]);
  Q[1] = currSolution->Get( eqNumbers[3]);
  //cout<<"element"<<" "<<element<<" "<<S[0]<<" "<<endl;

  // set the equation numbers on the element vector
  elementVector->SetEquationNumbers(eqNumbers);
  elementVector->Clear();

  for( int l = 0; l < quadPoints; l++){

    finiteElement->Evaluate( xi[l], shape, DxShape, &jacobian);

    // evaluate Un at integration point
    double Un[2];
    Un[0] = finiteElement->Interpolate( xi[l], Sn);
    Un[1] = finiteElement->Interpolate( xi[l], Qn);

    // evaluate U at integration point
    double U[2];
    U[0] = finiteElement->Interpolate( xi[l], S);
    U[1] = finiteElement->Interpolate( xi[l], Q);

    // evaluate U,z at integration point
    double DxU[2];
    DxU[0] = DxShape[0]*S[0]+DxShape[1]*S[1];
    DxU[1] = DxShape[0]*Q[0]+DxShape[1]*Q[1];

    // get position corresponding to quadrature point
    double z = finiteElement->Interpolate( xi[l], nodes);

    // more values coming from constitutive equations
    // IV added IntegralpS and IntegralpD2S 01-24-03
    double pressure = material->GetPressure( U[0], z);
    double Outflow = material->GetOutflowFunction( pressure, z);
    double DpDS = material->GetDpDS( U[0], z);
    double DpDz = material->GetDpDz( U[0], z);
    double IntegralpS = material->GetIntegralpS( U[0], z); //0.0;
    double IntegralpD2S = 0;

    if(cvOneDGlobal::CONSERVATION_FORM==1) {
      IntegralpD2S = material->GetIntegralpD2S( U[0], z); //0.0;
    }

    // evaluate the matrices A, C and K at U where CU=G
    // recall that:
    //    A11 = 0.0  A12 = 1.0
    //    C12 = 0.0
    //    K11 = 0.0  K12 = 0.0  K21 = 0.0
    // IV 01-24-03, used only for the stabilization term
    double aux = U[1]/U[0];        // aux = Q/S
    double A21 = -(1.0+delta)*aux*aux+U[0]/density*DpDS;
    double A22 = 2.0*(1.0+delta)*aux;//why??
    double C11 = - Outflow/U[0];
    double C21 = -1.0/density*DpDz;
    double C22 = N/U[0];
    double K22 = kinViscosity;
    // K22=0.0;

    double G1 = -Outflow;
    double G2 = N*aux-U[0]/density*DpDz;

    // IV's formulation 01-24-03
    // I am using F,K, and G which is different than Jing's G
    // so will call it GF in its contribution to RHS=-Residual
    // we begin by calculating the residual
    double F1  = U[1];
    double F2  = (1.0+ delta)*U[1]*aux+IntegralpS/density;
    // if(element==0) cout<<"F2"<<" "<<(1.0+ delta)*U[1]*aux<<" "<<IntegralpS/density<<endl;
    // double K22 = kinViscosity;
    double GF1 = -Outflow;
    double GF2 = N*aux+IntegralpD2S/density;

    // evaluate the tau matrix
    double tau[4];
    double h = nodes[1] - nodes[0];
    double A[4] = { 0.0, 1.0, A21, A22};
    double C[4] = { C11, 0.0, C21, C22};
    double modA[4];
    double modC[4];

    if(STABILIZATION==1){
      GetModulus( A, modA);
      if(kinViscosity < SMALL_KINEMATIC_VISCOSITY){
        modC[0] = 0;
        modC[1] = 0;
        modC[2] = 0;
        modC[3] = 0;
      }else{
        GetModulus(C,modC);
      }

      // this is actually the inverse of tau...
      tau[0] = (2.0/deltaTime) + 2.0/h*modA[0] + modC[0];
      tau[1] = 2.0/h*modA[1] + modC[1];
      tau[2] = 2.0/h*modA[2] + modC[2];
      tau[3] = (2.0/deltaTime)+2.0/h*modA[3]+12.0/(h*h)*K22+modC[3];

      double det = tau[0]*tau[3]-tau[1]*tau[2];

      // correcting tau...
      double temp = tau[0];
      tau[0] =  tau[3]/det;
      tau[1] = -tau[1]/det;
      tau[2] = -tau[2]/det;
      tau[3] =  temp/det;
    }

    for( int a = 0; a < numberOfNodes; a++){

      double jw = jacobian * weight[l];

      //
      // DG terms
      //
      double rDG1=0.0;
      double rDG2=0.0;

      if(cvOneDGlobal::CONSERVATION_FORM){

        // IV formulation 01-31-03

        rDG1 = deltaTime*(DxShape[a]*F1+shape[a]*GF1)-shape[a]*(U[0]-Un[0]);
        rDG2 = deltaTime*(DxShape[a]*F2-DxShape[a]*K22*DxU[1]+shape[a]*GF2)-shape[a]*(U[1]-Un[1]);
        // if(element==0 && a==0) cout<<"rDG2"<<" "<<DxShape[a]*F2<<" "<<-DxShape[a]*K22*DxU[1]<<" "<<shape[a]*GF2<<" "<<-shape[a]*(U[1]-Un[1])<<endl;
      }else{

        // Brooke's formulation that I am not using IV 01-31-03

        rDG1 = deltaTime*(shape[a]*(DxU[1])-shape[a]*G1)+shape[a]*(U[0]-Un[0]);
        rDG2 = deltaTime*(shape[a]*(A21*DxU[0]+A22*DxU[1])+DxShape[a]*K22*DxU[1]-shape[a]*G2)+shape[a]*(U[1]-Un[1]);
    }

    double rGLS1 = 0.0;
    double rGLS2 = 0.0;

    if (STABILIZATION == 1){
      //
      // GLS terms
      //
      // create an auxiliary matrix to handle some of the terms
      //
      double auxa[4];
      auxa[0] = -shape[a]*C11;    // A11 = 0.0
      auxa[1] = DxShape[a];    // A12 = 1.0, C12 = 0.0
      auxa[2] = DxShape[a]*A21-shape[a]*C21;
      auxa[3] = DxShape[a]*A22-shape[a]*C22;

      double auxb[2];
      auxb[0] = DxU[1]-G1;
      auxb[1] = A21*DxU[0]+A22*DxU[1]-G2;

      //
      //multiply the matrices to obtain the GLS contribution
      //into auxa
      //
      double auxc[4];  // contains the product: auxa * tau
      auxc[0] = auxa[0]*tau[0]+auxa[1]*tau[2];
      auxc[1] = auxa[0]*tau[1]+auxa[1]*tau[3];
      auxc[2] = auxa[2]*tau[0]+auxa[3]*tau[2];
      auxc[3] = auxa[2]*tau[1]+auxa[3]*tau[3];

      //
      // sum the GLS terms to the DG terms
      //
      rGLS1 = deltaTime*(auxc[0]*auxb[0]+auxc[1]*auxb[1]);
      rGLS2 = deltaTime*(auxc[2]*auxb[0]+auxc[3]*auxb[1]);

    }//end stabilization

      double r1 = rDG1 + rGLS1;
      double r2 = rDG2 + rGLS2;

      // now RHS= -residual , comment added by IV 01-24-03
      elementVector->Add( 2*a  , -r1*jw);
      elementVector->Add( 2*a+1, -r2*jw);
      // if(element==0 && a==0) cout<<"\t"<<"-r2*jw"<<" "<<-r2*jw<<endl;

    }
  }

  if(cvOneDGlobal::CONSERVATION_FORM){
    //Inlet flux term (at z=z_inlet) which is the linearized F-KU
    if(element == 0){
      long node = 0;
      double z = sub->GetNodalCoordinate( node);
      double aux = Q[0]/S[0];

      finiteElement->Evaluate(z,shape,DxShape,&jacobian);
      double dQdz = Q[1]*DxShape[1]+Q[0]*DxShape[0];
      double IntegralpS = material->GetIntegralpS(S[0],z);
      int a = 0;

      double InletR1 = Q[0];
      //double InletR2 = (1.0+delta)*Q[0]*aux + IntegralpS/density- kinViscosity*dQdz;
      double InletR2 = (1.0+delta)*Q[0]*aux + IntegralpS/density;//without viscosity in flux
      elementVector->Add(2*a  , -InletR1*deltaTime);
      elementVector->Add(2*a+1, -InletR2*deltaTime);
      // cout<<"\t"<<"\t"<<"-InletR2*deltaTime"<<" "<<-InletR2*deltaTime<<endl;
      // cout<<"\t"<<"\t"<<"terms of inletR2"<<" "<<(1.0+delta)*Q[0]*aux <<" "<<IntegralpS/density<<endl;
    }// end inlet flux

    if(bound==BoundCondTypeScope::NOBOUND||bound==BoundCondTypeScope::PRESSURE
      ||bound==BoundCondTypeScope::PRESSURE_WAVE||bound==BoundCondTypeScope::AREA||bound==BoundCondTypeScope::FLOW){
      //If no outlet BC or Dirichlet outlet BC, compute the Outlet full flux term (at z=z_outlet) which is the linearized F-KU IV 02-03-03
      if (element == (sub->GetNumberOfElements())-1){
        //cout<<element<<endl;
        double z = sub->GetOutletZ();//checked IV 02-03-03
        double pressure = material->GetPressure( S[1], z);
        finiteElement->Evaluate( z, shape, DxShape, &jacobian);//careful: shape is in the natural coord system (xi)

        int a = 1;

        double dQdz = Q[1]*DxShape[1]+Q[0]*DxShape[0];
        double IntegralpS = material->GetIntegralpS( S[1], z);
        double OutletR1= Q[1];
        // double OutletR2= (1.0+delta)*pow(Q[1],2)/S[1] + IntegralpS/density - kinViscosity*dQdz;

        // without viscosity in flux
        double OutletR2= (1.0+delta)*pow(Q[1],2)/S[1] + IntegralpS/density;
        // without adv term
        // double OutletR2= IntegralpS/density;
        // without M2H2
        // double OutletR2= 0;

        // try linear downstream domain-Hughes
        // double Cp=1.0;
        /*  double Cp = material->GetLinCompliance(z);
        // double Cp = material->GetnonLinCompliance( S[1],z);//tried 02-13-03 worse results
        double OutletR2 = S[1]*S[1]/(2.0*density*Cp) - pow(material->GetArea(material->p1,z),2)/(2*density*Cp);//linear downstream domain-Hughes
        */
        elementVector->Add( 2*a  , OutletR1*deltaTime);
        elementVector->Add( 2*a+1, OutletR2*deltaTime);
      }//end outlet flux term
    }//end if no outletBC or Dirichlet
  }//end   if(CONSERVATION_FORM)
}//end function
