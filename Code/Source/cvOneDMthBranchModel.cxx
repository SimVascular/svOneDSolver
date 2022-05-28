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
//  MthBranchModel.cxx
//  ~~~~~~~~~~~~~~~~~~
//
//  Build up the parts corresponding to the Lagrange nodes in the Jacobian
//  Matrix and the residual vector.
//

# include "cvOneDMthBranchModel.h"
# include "cvOneDUtility.h"

cvOneDMthBranchModel::cvOneDMthBranchModel(const vector<cvOneDSubdomain*>& subdList,
                                           const vector<cvOneDFEAJoint*>& jtList,
                                           const vector<int>& outletList) :
                                           cvOneDMthModelBase(subdList, jtList, outletList){
}

void cvOneDMthBranchModel::GetEquationNumbers(long ele, long* eqNumbers, long ith){
}

// For skyline matrix calculation
long cvOneDMthBranchModel::GetUpmostEqnNumber(long ele, long ithJoint){

  cvOneDFEAJoint* joint = jointList[ithJoint];
  cvOneDSubdomain* subdomain;
  int inletSegs = joint->getNumberOfInletSegments();
  long* eqnNumbers;
  int nonZeros;

  if(ele == 0){
    nonZeros = joint->getNumberOfSegments();
  }else{
    nonZeros = 2;
  }

  eqnNumbers = new long[nonZeros];
  int i, ID;

  if(ele == 0){
    for(i = 0; i < joint->getNumberOfInletSegments(); i++){
      int ID = joint->GetInletID(i);
      subdomain = subdomainList[ID];
      eqnNumbers[i] = (subdomain->GetGlobal1stNodeID()+subdomain->GetNumberOfNodes()-1)*2 + 1;
    }
    for(i = 0; i < joint->getNumberOfOutletSegments(); i++){
      ID = joint->GetOutletID(i);
      subdomain = subdomainList[ID];
      eqnNumbers[i+inletSegs] = subdomain->GetGlobal1stNodeID()*2 + 1;
    }
  }else{
    ID = joint->GetInletID(0);
    subdomain = subdomainList[ID];
    eqnNumbers[0] = (subdomain->GetGlobal1stNodeID()+subdomain->GetNumberOfNodes()-1)*2;
    if(ele/inletSegs == 0){ //they are from inletList
      int ID = joint->GetInletID(ele);
      subdomain = subdomainList[ID];
      eqnNumbers[1] = (subdomain->GetGlobal1stNodeID()+subdomain->GetNumberOfNodes()-1)*2;
    }else{
      int ID = joint->GetOutletID(ele-inletSegs);
      subdomain = subdomainList[ID];
      eqnNumbers[1] = subdomain->GetGlobal1stNodeID()*2;
    }
  }
  long minimum = min(nonZeros, eqnNumbers);
  delete [] eqnNumbers;
  return minimum;
}

void cvOneDMthBranchModel::FormNewton(cvOneDFEAMatrix* lhsMatrix, cvOneDFEAVector* rhsVector){
  for(long i = 0; i < jointList.size(); i++){
    FormLagrangeLHSbyP(i, lhsMatrix);
    FormLagrangeLHSbyQ(i, lhsMatrix);
    FormLagrangeRHSbyQ(i, rhsVector);
    FormLagrangeRHSbyP(i, rhsVector);
  }
}

//default: the first one of the inletSegment array will be the node
//to write mass balance (flow in = flow out), also the pressure of this
//node will be equal to the pressure of every other node
void cvOneDMthBranchModel::FormLagrangeRHSbyQ(long ith, cvOneDFEAVector* rhsVector){
  long eqn[2];
  long i;
  long lagOrigin;
  int currID = 0;

  for(i = 0; i < jointList[ith]->getNumberOfInletSegments(); i++){
    currID = jointList[ith]->GetInletID(i);
    long numNodes = subdomainList[currID]->GetNumberOfNodes();
    GetNodalEquationNumbers(numNodes-1, eqn, currID);
    if(i==0){
      lagOrigin = jointList[ith]->GetGlobal1stLagNodeID();
    }
    rhsVector->Add(eqn[1], -currSolution->Get(lagOrigin));
    double currQ = currSolution->Get(eqn[1]);
    rhsVector->Add(lagOrigin, - currQ);
  }

  for(i = 0; i < jointList[ith]->getNumberOfOutletSegments(); i++){
    currID = jointList[ith]->GetOutletID(i);
    GetNodalEquationNumbers(0, eqn, currID);
    rhsVector->Add(eqn[1], currSolution->Get(lagOrigin));
    double currQ = currSolution->Get(eqn[1]);
    rhsVector->Add(lagOrigin, currQ);
  }
}

//lagrange variables are always ordered as
// 0         1        2           3  ....
// lamda_Q  lamda_pi_in lamda_pi_out     .......
void cvOneDMthBranchModel::FormLagrangeRHSbyP(long ith, cvOneDFEAVector* rhsVector){
  cvOneDMaterial* material;
  cvOneDSubdomain* subdomain;
  long eqn[2];
  int ID = jointList[ith]->GetInletID(0);
  subdomain = subdomainList[ID];
  long numNodes = subdomain->GetNumberOfNodes();
  material = subdomain->GetMaterial();
  GetNodalEquationNumbers(numNodes-1, eqn, ID);
  long origin = eqn[0];
  double currS_orig = currSolution->Get(eqn[0]);
  double z_orig = subdomain->GetOutletZ();
  double pressure_orig = material->GetPressure(currS_orig, z_orig);
  double dpds_orig = material->GetDpDS(currS_orig, z_orig);
  int i;

  for(i=1;i<jointList[ith]->getNumberOfInletSegments();i++){

    int ID = jointList[ith]->GetInletID(i);
    subdomain = subdomainList[ID];
    long numNodes = subdomain->GetNumberOfNodes();
    material = subdomain->GetMaterial();
    GetNodalEquationNumbers(numNodes-1, eqn, ID);
    double lamdapi = currSolution->Get(jointList[ith]->GetGlobal1stLagNodeID()+i);
    double currS = currSolution->Get(eqn[0]);
    double z = subdomain->GetOutletZ();
    double dpidsi = material->GetDpDS(currS, z);
    rhsVector->Add(eqn[0], dpidsi*lamdapi);
    rhsVector->Add(origin, -dpds_orig*lamdapi);
    long lagID = jointList[ith]->GetGlobal1stLagNodeID()+i;
    double pressure = material->GetPressure(currS, z);
    rhsVector->Add(lagID, pressure - pressure_orig);
  }

  for(i=0;i<jointList[ith]->getNumberOfOutletSegments();i++){

    int ID = jointList[ith]->GetOutletID(i);
    subdomain = subdomainList[ID];
    long numNodes = subdomain->GetNumberOfNodes();
    material = subdomain->GetMaterial();
    GetNodalEquationNumbers(0, eqn, ID);
    long lagID = jointList[ith]->GetGlobal1stLagNodeID()+i+
    jointList[ith]->getNumberOfInletSegments();

    double lamdapi = currSolution->Get(lagID);
    double currS = currSolution->Get(eqn[0]);
    double z = subdomain->GetInletZ();
    double dpidsi = material->GetDpDS(currS, z);
    rhsVector->Add(eqn[0], dpidsi*lamdapi);
    rhsVector->Add(origin, -dpds_orig*lamdapi);
    double pressure = material->GetPressure(currS, z);
    rhsVector->Add(lagID, pressure - pressure_orig);
  }
}

void cvOneDMthBranchModel::FormLagrangeLHSbyP(long ith, cvOneDFEAMatrix* lhs){
  cvOneDMaterial* material;
  cvOneDSubdomain* subdomain;
  int currID = 0;
  long eqn[2];
  currID = jointList[ith]->GetInletID(0);
  subdomain = subdomainList[currID];
  material = subdomain->GetMaterial();
  long numNodes = subdomain->GetNumberOfNodes();
  GetNodalEquationNumbers(numNodes-1, eqn, currID);
  int origin = eqn[0];
  double currS_orig = currSolution->Get(eqn[0]);
  double z_orig = subdomain->GetOutletZ();
  double dpds_orig = material->GetDpDS(currS_orig, z_orig);
  double d2pds2_orig = material->GetD2pDS2(currS_orig, z_orig);
  int i;

  for(i = 1; i < jointList[ith]->getNumberOfInletSegments(); i++){

    currID = jointList[ith]->GetInletID(i);
    subdomain = subdomainList[currID];
    long numNodes = subdomain->GetNumberOfNodes();
    material = subdomain->GetMaterial();
    GetNodalEquationNumbers(numNodes-1, eqn, currID);
    long lagID = jointList[ith]->GetGlobal1stLagNodeID()+i;
    double lamdapi = currSolution->Get(lagID);
    double currS = currSolution->Get(eqn[0]);
    double z = subdomain->GetOutletZ();
    double dpds   = material->GetDpDS(currS, z);
    double d2pds2 = material->GetD2pDS2(currS, z);
    //add diagonal part
    lhs->AddValue(eqn[0], eqn[0], -d2pds2_orig*lamdapi);
    lhs->AddValue(origin, origin, d2pds2_orig*lamdapi);
    //add lag nodes part
    lhs->AddValue(eqn[0], lagID, -dpds);
    lhs->AddValue(origin, lagID, dpds_orig);
    lhs->AddValue(lagID, eqn[0], -dpds);
    lhs->AddValue(lagID, origin, dpds_orig);
  }

  for(i = 0; i < jointList[ith]->getNumberOfOutletSegments(); i++){

    currID = jointList[ith]->GetOutletID(i);
    subdomain = subdomainList[currID];
    long numNodes = subdomain->GetNumberOfNodes();
    material = subdomain->GetMaterial();
    GetNodalEquationNumbers(0, eqn, currID);
    long lagID = jointList[ith]->GetGlobal1stLagNodeID()+i+
    jointList[ith]->getNumberOfInletSegments();
    double lamdapi = currSolution->Get(lagID);
    double currS = currSolution->Get(eqn[0]);
    double z = subdomain->GetInletZ();
    double dpds   = material->GetDpDS(currS, z);
    double d2pds2 = material->GetD2pDS2(currS, z);
    lhs->AddValue(eqn[0], eqn[0], -d2pds2*lamdapi);
    lhs->AddValue(origin, origin, d2pds2_orig*lamdapi);
    //add lag nodes part
    lhs->AddValue(eqn[0], lagID, -dpds);
    lhs->AddValue(origin, lagID, dpds_orig);
    lhs->AddValue(lagID, eqn[0], -dpds);
    lhs->AddValue(lagID, origin, dpds_orig);
  }
}

void cvOneDMthBranchModel::FormLagrangeLHSbyQ(long ith, cvOneDFEAMatrix* lhs){
  int i;
  long eqn[2];
  long lagID = jointList[ith]->GetGlobal1stLagNodeID();

  for(i=0;i<jointList[ith]->getNumberOfInletSegments();i++){
    int ID = jointList[ith]->GetInletID(i);
    long numNodes = subdomainList[ID]->GetNumberOfNodes();
    GetNodalEquationNumbers(numNodes-1, eqn, ID);
    lhs->AddValue(eqn[1], lagID, 1);
    lhs->AddValue(lagID, eqn[1], 1);
  }

  for(i = 0; i < jointList[ith]->getNumberOfOutletSegments(); i++){
    int ID = jointList[ith]->GetOutletID(i);
    GetNodalEquationNumbers(0, eqn, ID);
    lhs->AddValue(eqn[1], lagID, -1);
    lhs->AddValue(lagID, eqn[1],  -1);
  }
}
