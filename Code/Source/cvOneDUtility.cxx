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

# include "cvOneDUtility.h"
#include <regex>
#include <string.h>

//
//  Utility.cxx - Source for Some Utility functions 
//  ~~~~~~~~~~~   
//

std::vector<std::string> split_string(std::string& s, const std::string& delims)
{
  std::vector<std::string> sub_strings;
  char* sub_str;
  sub_str = strtok(&s[0], delims.c_str());

  while (sub_str != NULL) {
    auto token = std::string(sub_str);
    if ((token.back() == '\n') || (token.back() == '\r')) {
      token.pop_back();
    } 
    if ((token != "") && (token[0] != '\n') && (token[0] != '\r')) {
      sub_strings.push_back(token);
    }
    sub_str = strtok(NULL, delims.c_str());
  }
  return sub_strings;
}

std::string trim_string(const std::string& s)
{
  return std::regex_replace(s, std::regex("^ +| +$|( ) +"), "$1");
}

std::string upper_string(const std::string& s)
{
  std::string upper_str(s);
  std::transform(upper_str.begin(), upper_str.end(), upper_str.begin(), [](unsigned char c){ return std::toupper(c); });
  return upper_str;
}


long min(long size, long* values){
  long m = values[0];
  long nm;
  for( long i = 0; i < size; i++){
    nm = values[i];
    m = (m < nm) ? m : nm;
  }
  return m;
}

long min(long a, long b){
  return (a < b) ? a : b;
}

long max(long a, long b){
  return (a < b) ? b : a;
}

void clear(long size, long* vector){
  for( long i = 0; i < size; i++){
    vector[i] = 0;
  }
}

long sum(long size, long* values){
  long s = 0;
  for( long i = 0; i < size; i++){
    s += *values++;
  }
  return s;
}

void GetModulus( double* A, double* modulusA){
  double A2[4];	    // A2 = A*A
  A2[0] = A[0] * A[0] + A[1] * A[2];
  A2[1] = A[0] * A[1] + A[1] * A[3];
  A2[2] = A[2] * A[0] + A[3] * A[2];
  A2[3] = A[2] * A[1] + A[3] * A[3];
  
  double traceA2 = A2[0] + A2[3];
  double detA = A[0] * A[3] - A[1] * A[2];
  double b = fabs(detA);
  double a = sqrt(traceA2 + 2.0 * b);

  if(fabs(a)>1.0e-8){
    modulusA[0] = (A2[0] + b) / a;
    modulusA[1] = (A2[1]    ) / a;
    modulusA[2] = (A2[2]    ) / a;
    modulusA[3] = (A2[3] + b) / a;
  }else{
    modulusA[0] = 0.0;
    modulusA[1] = 0.0;
    modulusA[2] = 0.0;
    modulusA[3] = 0.0;
  }
}

/*Quadrature::Quadrature( int n)
{
  double p1, w1, w2;
  
  assert( n >= 1 && n <= 3);
  
  npts = n;
  xi = new double[ npts];
  weights = new double[ npts];
  
  assert( xi != 0 && weights != 0);
  
  switch( npts){
  case 1:
    xi[0] = 0.0;
    weights[0] = 2.0;
    break;
  case 2:
    p1 = 1.0/sqrt(3.0);
    xi[0] = -p1;
    xi[1] = p1;
    weights[0] = weights[1] = 1.0;
    break;
  case 3:
    p1 = sqrt(3.0/5.0);
    w1 = 5.0/9.0;
    w2 = 8.0/9.0;
    xi[0] = -p1;
    xi[1] = 0.0;
    xi[2] = p1;
    weights[0] = weights[2] = w1;
    weights[1] = w2;
    break;
  case 4:
    xi[0] = - 0.861136311594953;
    xi[1] = - 0.339981043584856;
    xi[2] = - xi[0];
    xi[3] = - xi[1];
    weights[0] = weights[2] = 0.347854845137454;
    weights[1] = weights[3] = 0.652145154862546;
    
  }
}
*/
/* Computes the (np+1)/2 non-negative abscissas x[i] & corresponding weights 
   w[i] of the n-point Gauss-Legendre integration rule, normalized to 
   interval [-1,1]. The abscissas appear in ascending order.   
   Ref: Gaussian Quadrature Formulas, Stroud and Secrest, 1966. 
*/ 
cvOneDQuadrature::cvOneDQuadrature(int np){ 
  double pkm1, pk, t1, pkp1, den, d1, dpn, d2pn, d3pn, d4pn, u, v, h, p, dp, fx;
  double tx[25], tw[25]; 
  
  int i, m = (np+1)/2; 
  double e1 = np*(np+1); 
  double t, x0; 
  
  npts = np;
  xi = new double[np];
  weights = new double[np];
  
  assert( xi != 0 && weights != 0);
  
  for (i=1; i<=m; i++){ 
    t = (4*i-1)*3.1415926535897932384626433/(4*np+2); 
    x0 = (1.-(1.-1./np)/(8.*np*np))*cos(t); 
    pkm1 = 1.; 
    pk = x0; 
    for (int k=2; k<=np; k++){ 
      t1 = x0*pk; 
      pkp1 = t1 - pkm1 - (t1-pkm1)/k + t1; 
      pkm1 = pk; 
      pk = pkp1; 
    } 
    den = 1. - x0*x0; 
    d1 = np*(pkm1-x0*pk); 
    dpn = d1/den; 
    d2pn = (2.*x0*dpn-e1*pk)/den; 
    d3pn = (4.*x0*d2pn+(2.-e1)*dpn)/den; 
    d4pn = (6.*x0*d3pn+(6.-e1)*d2pn)/den; 
    u = pk/dpn; 
    v = d2pn/dpn; 
    h = -u*(1.+0.5*u*(v+u*(v*v-d3pn/(3.*dpn)))); 
    p = pk + h*(dpn+0.5*h*(d2pn+h/3.*(d3pn+0.25*h*d4pn))); 
    dp = dpn + h*(d2pn+0.5*h*(d3pn+h*d4pn/3.)); 
    h -= p/dp; 
    tx[i-1] = x0 + h; 
    fx = d1 - h*e1*(pk+0.5*h*(dpn+h/3.*(d2pn+0.25*h*(d3pn+0.2*h*d4pn)))); 
    tw[i-1] = 2.*(1-tx[i-1]*tx[i-1])/(fx*fx); 
  } 
  if (2*m>np) 
    tx[m-1] = 0.0; 
  for (i=0; i<m; i++) { 
    xi[i] = tx[i];        
    weights[i] = tw[i];   
  } 
  for(i = 0; i < m-np%2; i++) {
    xi[np-i-1] = -xi[i];
    weights[np-i-1] = weights[i];
  }
  return; 
} 
 
 cvOneDQuadrature::~cvOneDQuadrature()
{
  delete [] xi;
  delete [] weights;
}


void cvOneDQuadrature::Get(double* w, double* x)const{
  for( int i = 0; i < npts; i++){
    w[i] = weights[i];
    x[i] = xi[i];
  }
}

// GET INDEX FROM STRING VECTOR
int getListIDWithStringKey(string key,cvStringVec list){
  bool found = false;
  size_t count = 0;
  while((!found)&&(count<list.size())){    
    found = upper_string(list[count]) == key;
    //printf("FOUND: %s\n",found ? "True":"False");
    if(!found){
      count++;
    }
  }
  if(found){
    return count;
  }else{
    return -1;
  }
}
