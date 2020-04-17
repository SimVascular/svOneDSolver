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

#ifndef CVONEDERROR_H
#define CVONEDERROR_H

//
//  cvOneDError.h: Header file for cvOneDError.handling/Debugging
//  ~~~~~~~~~
//
//  NOTES:
//
//  If you want to use this class to handle errors in your code
//  (and you should becase it provides a consistant way to do 
//  this throughout the code and any applications derived from it),
//  Simply set the appropriate Error Number with setErrorNumber(...)
//  and a *descriptive* error string using setErrorString(...)
//  If you want to handle the error yourself (ie. you don't think
//  the application programmer should have to do this explicitly, 
//  such as a max iterations exceeded, etc...), call the error 
//  handler -- which will be defined by the programmer in his/her 
//  code.
//  
 
#include "cvOneDEnums.h"

typedef void (*ErrorHandlerFunctionPtr)();

class cvOneDError{
  public: // these are used by applications
    static ErrorType getErrorNumber(void);
    static char*     getErrorString(void);

    static ErrorHandlerFunctionPtr setErrorHandler(ErrorHandlerFunctionPtr fcn);

    static void setDebugLevel(int level);
    static int  getDebugLevel(void);

 public: // these should be used predominantly by API

    static void setErrorNumber(ErrorType errorNumber);
    static void setErrorString(char *string);
    static void CallErrorHandler(void);

 private:

    cvOneDError(void); // prohibit instancing

    // error handling state info
    static ErrorType ErrorNo;
    static char* ErrorString;
    static ErrorHandlerFunctionPtr errorHandler;

    // debugging state info
    static int debugLevel;
}; 

#endif // CVONEDERROR_H
