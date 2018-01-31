//
//  bnsString.cxx - Source for a simple string class.
//

# include <string.h>
# include <assert.h>

# include "cvOneDString.h"

cvOneDString::cvOneDString(){
  len = 0;
  theData = 0;
}

cvOneDString::cvOneDString(const char* str){
  len = strlen( str);
  theData = new char[len + 1];
  assert( theData != 0);
  strcpy( theData, str);
}

cvOneDString::cvOneDString(const cvOneDString& str){
  len = str.len;
  theData = new char[len + 1];
  assert( theData != 0);
  strcpy( theData, str.theData);
}

cvOneDString::~cvOneDString(){
  delete [] theData;
}

const char* cvOneDString::data(){
  return theData;
}

cvOneDString cvOneDString::operator+(const char* rhs){
  cvOneDString str = *this;
  str += rhs;
  return str;
}

cvOneDString& cvOneDString::operator+=(const char* rhs){
  len += strlen( rhs);
  char *p = new char[ len + 1];
  assert( p != 0);
  strcpy( p, theData);
  strcat( p, rhs);
  delete [] theData;
  theData = p;

  return *this;
}

const cvOneDString& cvOneDString::operator=(const cvOneDString& rhs){
  len = rhs.len;
  delete [] theData;
  theData = new char[len + 1];
  strcpy( theData, rhs.theData);

  return *this;
}
