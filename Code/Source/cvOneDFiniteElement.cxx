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
//  FiniteElement.cxx - Source for a class to handle individual finite elements
//  ~~~~~~~~~~~~~~~~~
//  
//  This class abstracts the finite element one-dimensional shape function
//

#include <cassert>

#include "cvOneDFiniteElement.h"

// Constuctor
cvOneDFiniteElement::cvOneDFiniteElement(){
  wasSet = false;
}

// Destructor
cvOneDFiniteElement::~cvOneDFiniteElement(){

}

void cvOneDFiniteElement::Set(double* nd, long* conn){
  nodes[0] = nd[0];
  nodes[1] = nd[1];
  
  connectivity[0] = conn[0];
  connectivity[1] = conn[1];
  
  wasSet = true;
}

void cvOneDFiniteElement::Evaluate(double xi, double* shape, double* DxShape, double* jacobian)const{
  assert( wasSet);
  
  shape[0] = 0.5 * (1.0 - xi);
  shape[1] = 0.5 * (1.0 + xi);
  
  *jacobian = 0.5 * (nodes[1] - nodes[0]);
  
  DxShape[0] = -0.5 / (*jacobian);
  DxShape[1] = 0.5 / (*jacobian);
}

double cvOneDFiniteElement::Interpolate(double xi, double* values)const{
  double shape[2];
  double aux1[2]; // redundant pointer
  double aux2;    // redundant value
  
  Evaluate(xi, shape, aux1, &aux2);
  
  return values[0] * shape[0] + values[1] * shape[1];
}
