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
//  cvOneDString.cpp - Source for a simple string class.
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
