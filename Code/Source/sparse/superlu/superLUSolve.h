#ifndef SUPERLUSOLVE_H
#define SUPERLUSOLVE_H

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "cvOneDTypes.h"

int superLUSolve(cvOneDKentry* Kentries,  double *b,int numEntries, int NNZ,int nunknown,double u[]);

#endif // SUPERLUSOLVE_H
