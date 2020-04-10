#ifndef CVONEDFEAVECTOR_H
#define CVONEDFEAVECTOR_H

//
//  cvOneDFEAVector.h - Header for a class to handle finite element Vectors
//  ~~~~~~~~~~~~~
//  
//  SYNOPSIS...This is a utility class used to handle the finite element 
//             Vectors (currentSolution, increment, previousSolution, etc.
//

#include <fstream>

#include "cvOneDSolverDefinitions.h"

using namespace std;

enum normType { 
  L2_norm = 1, 
  L1_norm, 
  Linf_norm
};

class cvOneDFEAVector{
  protected:
    char title[MAX_STRING_SIZE + 1];
    long dimension;
    long* equationNumbers;
    double* entries;
    void CreateVector(long dim, const char* tit = "vector");

  public:
    cvOneDFEAVector(long dim, const char* tit = "vector");
    cvOneDFEAVector(long dim, long* eqNumbers, const char* tit = "vector");
    void SetEquationNumbers(long* eqNumbers);
    ~cvOneDFEAVector();
    void Clear();
    void Rename(const char* tit);
    void Add(long component, double value);
    void Add(cvOneDFEAVector& vector);				
    // add the components of vector in the global vector in the 
    // positions specified by its equation numbers
    long GetDimension() const {return dimension;}
    const long* GetEquationNumbers() const {return equationNumbers;}
    double* GetEntries() {return entries;}
    double& operator[](long i);
    cvOneDFEAVector& operator=(const cvOneDFEAVector& rhs);
    cvOneDFEAVector& operator+=(const cvOneDFEAVector& rhs);
    void Set(long i, double value);
    double Get( long i) const;
    double Norm(normType type, int start, int step, int stop_index = -1) const;
    void CheckPositive(int start,int step, int stop);
	
    void print(std::ostream& os);
};

#endif	// CVONEDFEAVECTOR_H
