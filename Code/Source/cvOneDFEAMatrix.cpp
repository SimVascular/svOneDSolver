//
//  SkylineMatrix.cxx - Solve Skyline Matrices
//

# include <cassert>
# include <iostream>
# include "cvOneDSkylineMatrix.h"
# include "cvOneDDenseMatrix.h"

using namespace std;

cvOneDFEAMatrix::cvOneDFEAMatrix(const char* tit){
  int i = 0;
  while( i < MAX_STRING_SIZE && tit[i] != '\0'){
    title[i] = tit[i];
    i++;
  }
  title[i] = '\0';
}

cvOneDFEAMatrix::~cvOneDFEAMatrix(){
}


