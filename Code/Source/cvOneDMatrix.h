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

#ifndef CVONEDMATRIX_H
#define CVONEDMATRIX_H

//
//  cvOneDMatrix.h: A Custom templated matrix class
//
//  SYNOPSIS...Templated matrix class, extends cvOneDVector.h to two 
//             dimensional "safe" (range-checked) matrices.  
//
//             Originally Written 4.17.95
//
//  USAGE......For a matrix of Items use Matrix<Item>, e.g., 
//             Matrix <int> intmatrix; note: Item must have a 
//             default constructor.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// constructors:
//   Matrix()                   -- default, matrix of size 0x0
//   Matrix(int rows,int cols)  -- matrix with dimensions rows x cols
//   Matrix(int rows,int cols
//          Item fillValue)     -- matrix w/ entries all == fillValue
//   Matrix(const Matrix & vec) -- copy constructor
//    
//   int Rows()                  -- returns # of rows (capacity)
//   int Cols()                  -- returns # of columns (capacity
//                                                        
//   void SetSize(int rows, int cols) -- resizes the matrix to rows x cols
//                                   can lose entries if size is smaller
//                                   in either dimension
//
//   void Fill(Item fillValue)    -- set all entries == fillValue
//
//   operator =                   -- assignment operator works properly
//   operator []                  -- indexes const and non-const matrixs
//
//
//  examples of use:
//        Matrix<double> dmat(100,80);      // 100 x 80 matrix of doubles
//        Matrix<double> dzmat(100,80,0.0); // initialized to 0.0
//
//        Matrix<String> smat(300,1);       // 300 strings
//
//        Matrix<int> imat;                 // has room for 0 ints
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#include <cstdlib>
#include <cassert>
#include <iostream>

#include "cvOneDVector.h"

template <class Item> class cvOneDMatrix{

  public:

    // postcondition: matrix of zero items constructed  
    cvOneDMatrix(): myMatrix(0),myRows(0),myCols(0){}
    
    // postcondition: matrix w/ dimensions rows x cols constructed   
    cvOneDMatrix(int rows,int cols): myMatrix(rows),myRows(rows),myCols(cols){
      int k;
      for(k=0;k<rows;k++){
        myMatrix[k].SetSize(cols);
       }
    }

    // Postcondition: matrix w/ dimensions rows x cols,
    //                (each entry initialized to fillValue) constructed
    cvOneDMatrix(int rows, int cols, Item fillValue):myMatrix(rows),myRows(rows),myCols(cols){
      int k;
      for(k=0;k<rows;k++){
        myMatrix[k].SetSize(cols);
        myMatrix[k].Fill(fillValue);
      }
    }

    // Copy constructor
    // Precondition: Item supports assignment
    // Postcondition: make copy of mat
    cvOneDMatrix(const cvOneDMatrix<Item> & mat):myMatrix(mat.myMatrix),myRows(mat.Rows()),myCols(mat.Cols()){
      // initialier list takes care of things
    }
    
    // free new'd storage
    // postcondition: dynamically allocated storage freed
    ~cvOneDMatrix (){
      // leave in "empty" state
      myRows = 0;
      myCols = 0;                  
    }

    // Postcondition: all entries == fillvalue  
    void Fill(Item fillValue){
      int k;
      for(k=0;k<myRows;k++){
        myMatrix[k].Fill(fillValue);
      }
    }

    // Overload assignment
    // Precondition: Item supports assignment     
    // Postcondition: self is assigned mat
    cvOneDMatrix & operator = (const cvOneDMatrix<Item> & mat){
      // don't assign to self!
      if (this != &mat){
        myMatrix = mat.myMatrix;
        myRows = mat.myRows;
        myCols = mat.myCols;
      }
      return *this;           
    }
    
    // Row capacity of matrix
    int Rows()const{return myRows;}
    // Column capacity of matrix
    int Cols()const{return myCols;}

    // precondition: 0 < rows and 0 < cols
    // postcondition: matrix has dimensions rows x cols
    //                the first rows of which are copies of original
    //                unless rows < original rows, whence truncation
    //                same with cols (copy as much as possible)    
    void SetSize(int rows, int cols){
      int k;
      myMatrix.SetSize(rows);
        
      for(k=0;k<rows;k++){
        myMatrix[k].SetSize(cols);
      }
      myRows = rows;
      myCols = cols;
    }

    // safe indexing, returning reference
    // precondition: 0 <= index < myRows
    // postcondition: return index-th item
    // exception: aborts if index is out-of-bounds
    cvOneDVector<Item> & operator [] (int index){
      if ((unsigned) index >= (unsigned) myRows || index < 0){
        // cerr << "Illegal matrix index: " << index
        //      << " (max = " << myRows-1 << ")" << endl;
        assert(index >= 0);
        assert(index < myRows);         
      }
      return myMatrix[index];         
    }

    // safe indexing, returning const reference to avoid modification
    // precondition: 0 <= index < myRows
    // postcondition: return index-th item
    // exception: aborts if index is out-of-bounds    
    // const index 
    const cvOneDVector<Item> & operator [] (int index)const{
      if((unsigned) index >= (unsigned) myRows || index < 0){
        cerr << "Illegal matrix index: " << (unsigned)index
             << " (max = " << myRows << ")" << endl;
        assert(index >= 0);
        assert(index < myRows);         
      }
      return myMatrix[index]; 
    }

    Item ** getPointer(void){  
      cvOneDVector<Item> *temp = myMatrix.getPointer();
      Item **K = new Item* [myCols];

      for (int i=0;i<myRows;i++){
        K[i] = temp[i].getPointer();
      }
      return K;
    }

  private:
    
    cvOneDVector<cvOneDVector<Item> > myMatrix; // the matrix of items
    int myRows;                     // # of rows (capacity)
    int myCols;                     // # of cols (capacity)
};

#endif // CVONEDMATRIX_H

