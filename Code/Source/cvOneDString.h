#ifndef CVONEDSTRING_H
#define CVONEDSTRING_H

//
//  cvOneDString.h - Header for a simple string class.
//

class cvOneDString{

  public:

    cvOneDString();
    cvOneDString(const char* str);
    cvOneDString(const cvOneDString& str);
    ~cvOneDString();
    const char* data();

    const cvOneDString& operator=(const cvOneDString& rhs);
    cvOneDString operator+(const char* rhs);
    cvOneDString& operator+=(const char* rhs);

  private:

    long len;
    char* theData;
};

#endif // CVONEDSTRING_H

