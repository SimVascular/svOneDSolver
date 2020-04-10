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
