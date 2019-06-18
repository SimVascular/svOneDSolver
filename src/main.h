# include <string.h>
# include <iostream>
# include <fstream>

# include <boost/algorithm/string.hpp>

# include "cvOneDGlobal.h"
# include "cvOneDUtility.h"
# include "cvOneDModelManager.h"
# include "cvOneDOptions.h"
# include "cvOneDDataTable.h"
# include "cvOneDException.h"

using namespace std;

void WriteHeader();
int  getDataTableIDFromStringKey(string key);
void createAndRunModel(cvOneDOptions* opts);
void readModelFile(string inputFile, cvOneDOptions* opts, cvStringVec includedFiles);
void readModel(string inputFile, cvOneDOptions* opts);
void runOneDSolver(string inputFile);

