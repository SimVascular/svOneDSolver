%module oneDSolver

 %{
 /* Includes the header in the wrapper code */
 #include "cvOneDOptions.h"
 #include "main.h"
 using namespace std;
 %}

 /* Parse the header file to generate wrappers */
 %include <typemaps.i>
 %include <std_string.i>
 %include <std_vector.i>
 %include <std_map.i>

 %apply const int & { int & }; 
 %apply const double & { double & }; 
 
 %include "cvOneDOptions.h"
 %include "main.h"

 namespace std{
   typedef std::string String;
   // Vectors
   %template(cvStringVec) vector< string >;
   %template(cvLongVec)   vector< long >;
   %template(cvDoubleVec) vector< double >;
   // Matrices
   %template(cvStringMat) vector< vector< string > >;
   %template(cvLongMat)   vector< vector< long > >;
   %template(cvDoubleMat) vector< vector< double > >;
 }  

