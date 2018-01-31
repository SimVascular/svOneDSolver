#ifndef CVONEDEXCEPTION_H
#define CVONEDEXCEPTION_H

#include <string>

using namespace std;

class cvException: public exception{
  public:
    // Constructor and Destructor
    cvException(const char* m):msg(m){}
    virtual ~cvException() throw(){}
    // Member Functions
	virtual const char* what() const throw() {return msg.c_str();}
  protected:
    string msg;
};

#endif // CVONEDEXCEPTION_H
