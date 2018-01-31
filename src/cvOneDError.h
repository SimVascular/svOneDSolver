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
//  History:
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
