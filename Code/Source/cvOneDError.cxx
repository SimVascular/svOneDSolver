//
//  Error.cxx: Source file for cvOneDError.handling/Debugging
//  ~~~~~~~~~
//

#include "cvOneDError.h"

using namespace std;

// Declare the static globals
ErrorType  cvOneDError::ErrorNo;
char*  cvOneDError::ErrorString;
ErrorHandlerFunctionPtr cvOneDError::errorHandler;
int cvOneDError::debugLevel;

// Functions to be used by applications
ErrorType cvOneDError::getErrorNumber(void){
  return ErrorNo;
}

char* cvOneDError::getErrorString(void){
  return ErrorString;
}

ErrorHandlerFunctionPtr cvOneDError::setErrorHandler(ErrorHandlerFunctionPtr fcn){
  errorHandler = fcn;
  return errorHandler;
}

void cvOneDError::setDebugLevel(int level){
  debugLevel = level;
}

int cvOneDError::getDebugLevel(void){
  return debugLevel;
}

// Functions to be used by API
void cvOneDError::setErrorNumber(ErrorType errorNumber){
  ErrorNo = errorNumber;
}

void cvOneDError::setErrorString(char *string){
  ErrorString = string;
}

void cvOneDError::CallErrorHandler(void){
  if(errorHandler == NULL){
    cerr << "There's an error, but No error handler!" << endl;
    return;
  }
  errorHandler();
}
