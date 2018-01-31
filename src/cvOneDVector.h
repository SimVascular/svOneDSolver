#ifndef CVONEDVECTOR_H
#define CVONEDVECTOR_H

//
//  cvOneDVector.h: This is a custom templated vector class
//
//  SYNOPSIS...templated vector class, partially based on Budd, Classic 
//             Data Structures in C++. Originally written 11/5/93, 
//             modified 3/23/94.  Changed on 3/9/95 to make all methods 
//             inline (defined in class decl)
//
//  USAGE......for a vector of Items use Vector<Item>, e.g., 
//             Vector <int> intvector.  Note: Item must have a default 
//             constructor.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// constructors:
//   Vector()                     -- default, vector of size 0 (no entries)
//   Vector(int size)             -- vector with size entries
//   Vector(int size,
//          Item fillValue)       -- vector w/ size entries all == fillValue
//   Vector(const Vector & vec)   -- copy constructor
//    
//   int Length()                 -- returns size of vector (capacity)
//   void SetSize(int newSize)    -- resizes the vector to newSize elements
//                                  (can result in losing elements if
//                                   new size < old size)
//   void resize(int newSize)     -- synonym for SetSize
//
//   void Fill(Item fillValue)    -- set all entries == fillValue
//
//   operator =                   -- assignment operator works properly
//   operator []                  -- indexes both const and non-const vectors
//    
//
//  examples of use:
//             Vector<double> dlist(100);       // a list of 100 doubles
//             Vector<double> dzlist(100,0.0);  // initialized to 0.0
//
//             Vector<String> slist(300);       // 300 strings
//
//             Vector<int> ilist;               // has room for 0 ints
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <cstdlib>
#include <cassert>
#include <iostream>

template <class Item> class cvOneDVector{
  public:
    cvOneDVector();                               // default constructor 0 elts
    cvOneDVector(int size);                       // specify size of vector
    cvOneDVector(int size, Item fillValue);      // specify size and fill value
    cvOneDVector(const cvOneDVector<Item> & vec);             // copy constructor
    ~cvOneDVector ();                            // free new'd storage
    cvOneDVector & operator = (const cvOneDVector<Item> & vec); // overload assignment
    int Length() const;                    // capacity of vector
    int length() const;
    void Fill(Item fillValue);
    void SetSize(int newSize);             // change size dynamically
    void resize(int newSize);
    Item & operator [] (int index);      
    const Item & operator [] (int index) const; // const index 
    Item * getPointer(void);
  private:
    Item * myList;   // the array of items
    int myLength;   // # things in vector (array), 0,1,...,(myLength-1)
};

template<class Item> cvOneDVector<Item>::cvOneDVector(){
  myLength = 0; myList = 0;
}

template<class Item> cvOneDVector<Item>::cvOneDVector(int size){
  myLength = size;
  myList = new Item [size];
  assert(myList != 0);
}

template<class Item> cvOneDVector<Item>::cvOneDVector(int size, Item fillValue){
  myLength = size;
  myList = new Item [size];
  assert(myList != 0);
  for(int k = 0; k < size; k++){
    myList[k] = fillValue;
  }       
}

// copy constructor
// precondition: Item supports assignment
// postcondition: return copy of vec        
template<class Item>cvOneDVector<Item>::cvOneDVector(const cvOneDVector<Item> & vec){
  // allocate storage  
  myList = new Item [myLength = vec.myLength];
  assert(myList != 0);  
  // copy elements
  for(int k=0;k<vec.myLength;k++){
    myList[k] = vec.myList[k];
  }       
}

template<class Item> cvOneDVector<Item>::~cvOneDVector(){
  delete [] myList;
  myList = 0;
  // leave in "empty" state
  myLength = 0;
}

template<class Item> cvOneDVector<Item>& cvOneDVector<Item>::operator= (const cvOneDVector<Item>& vec){
  // don't assign to self!
  if(this != &vec){       
    delete [] myList;            // out with old list, in with new
    myList = new Item [myLength = vec.myLength]; 
    assert(myList != 0);
      
    // copy vec
    myLength = vec.myLength;
    
    for(int k=0; k < myLength; k++){
	    myList[k] = vec.myList[k];
	  }
  }
  return *this;        
}

template<class Item> int cvOneDVector<Item>::Length() const{
  return myLength;
}

template<class Item> int cvOneDVector<Item>::length() const{
  return myLength;
}

// Postcondition: all entries == fillvalue
template<class Item> void cvOneDVector<Item>::Fill(Item fillValue){
  int k;
  for(k=0; k < myLength; k++){
    myList[k] = fillValue;
  }
}

// change size dynamically
// precondition: vector has room for myLength entries     
// postcondition: vector has room for newSize entries
//                the first myLength of which are copies of original
//                unless newSize < myLength, then truncated copy occurs       
template<class Item> void cvOneDVector<Item>::SetSize(int newSize){
  int numToCopy = newSize < myLength ? newSize : myLength;
  
  // allocate new storage
  Item * newList = new Item[newSize];
  assert(newList != 0);   // be sure storage allocated
  
  int k;
  for(k=0;k<numToCopy;k++){
    newList[k] = myList[k];
  }
  
  delete [] myList;         // de-allocate old storage
  myLength = newSize;
  myList = newList;
}


template<class Item> void cvOneDVector<Item>::resize(int newSize){
  SetSize(newSize);
}

// safe indexing, returning reference
// precondition: 0 <= index < myLength
// postcondition: return index-th item
// exception: aborts if index is out-of-bounds
template<class Item> Item & cvOneDVector<Item>::operator [] (int index){
  if ((unsigned) index >= (unsigned)myLength || index < 0){
    //nmw      cerr << "Illegal vector index: " << index
    //nmw    << " (max = " << (myLength-1) << ")" << endl;
    assert(index >= 0);
    assert(index < myLength);
  }
  return myList[index];     
}

// const index 
// safe indexing, returning const reference to avoid modification
// precondition: 0 <= index < myLength
// postcondition: return index-th item
// exception: aborts if index is out-of-bounds
template<class Item> const Item & cvOneDVector<Item>::operator [] (int index)const{  
  if ((unsigned) index >= (unsigned)myLength || index < 0){
    cerr << "Illegal vector index: " << (unsigned)index << " (max = " << myLength << ")" << endl;
    assert(index >= 0);
    assert(index < myLength);
  }  
  return myList[index]; 
}

template<class Item> Item * cvOneDVector<Item>::getPointer(void){
  return myList;
}
    
#endif // CVONEDVECTOR_H
