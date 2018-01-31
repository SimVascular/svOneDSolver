#ifndef CVONEDUTILITY_H
#define CVONEDUTILITY_H

//
//  cvOneDUtility.h - Header for Some Utility functions and rules of random point quadratures.
//  ~~~~~~~~~~~   
//  
//  History:
//  May 1999, J.Wan
//      Creation of file,

# include <cmath>
# include <iostream>
# include <cstdarg>
# include <cassert>
# include <string>

# include <boost/algorithm/string.hpp>

# include "cvOneDTypes.h"
# include "cvOneDException.h"

const int MaxChar = 128;

long min(long a, long b);
long max( long a, long b);	
long min( long size, long* values);
long sum( long size, long* values);	
void clear( long size, long* vec);	
// Calculates the modulus of the 2x2 matrix A and put the results
// in modulusA using Cayley-Hamilton theory
void GetModulus(double* A, double* modulusA);

class cvOneDQuadrature{
  public:
    cvOneDQuadrature( int n);
    ~cvOneDQuadrature();
    void Get( double* w, double* x) const;

  private:	
  	// number of integration points in the quadrature
    int npts;  
    double* xi;
    double* weights;
};

int getListIDWithStringKey(string key,cvStringVec list);

#endif // CVONEDUTILITY_H
