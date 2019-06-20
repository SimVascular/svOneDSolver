#ifndef CVONEDTYPES_H
#define CVONEDTYPES_H

# include <vector>
# include <string>

using namespace std;

// CONSTANTS
const int CV_STRLEN = 256;

// RETURN VALUES
const int CV_ERROR = -1;
const int CV_OK = 0;

// TYPES

// STRUCTS
struct cvOneDKentry{
  int row;
  int col;
  double value;
};

// Vectors
typedef vector<string> cvStringVec;
typedef vector<long>   cvLongVec;
typedef vector<double> cvDoubleVec;

// Matrices
typedef vector<vector<string> > cvStringMat;
typedef vector<vector<long> >   cvLongMat;
typedef vector<vector<double> > cvDoubleMat;

#endif // CVONEDTYPES_H

