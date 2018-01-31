
#include <iostream>
#include <fstream>
#include "GenBC0D1D.h"


using namespace std;



int main(int argc, char** argv){


	cGenBC* genBC = new cGenBC();

	genBC->GenBC0D1D();


	delete genBC;

	return 0;

}