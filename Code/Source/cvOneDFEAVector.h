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
