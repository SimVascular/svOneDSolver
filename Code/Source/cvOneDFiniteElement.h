#ifndef CVONEDFINITEELEMENT_H
#define CVONEDFINITEELEMENT_H

//
//  cvOneDFiniteElement.h - Header for a class to handle individual finite elements
//  ~~~~~~~~~~~~~~~~~
//  
//  SYNOPSIS...This class abstracts the finite element and provides utilities
//             for maintaining a solution.
//

class cvOneDFiniteElement{

  public:
    
    cvOneDFiniteElement();
    virtual ~cvOneDFiniteElement();
    void Set( double* nd, long* conn);
    void Evaluate( double xi, double* shape,double* DxShape, double* jacobian)const;

    double Interpolate( double xi, double* values)const;
 
  private:
    
    bool wasSet;
    double nodes[2];
    long connectivity[2];
};

#endif // CVONEDFINITEELEMENT_H
