// CRUNCH_bins.cpp
// implementation of the CRUNCH_bins object
// it is responsible for aggregating geochemical data
// across multiple columns and also includes printing 
// capabilities


#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include "mathutil.hpp"
#include "CRN_parameters.hpp"
#include "CRN_tParticle_bins.hpp"
#include "CRUNCH_bins.hpp"
using namespace std;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::create()
{
  cout << "You need to at least provide the number of bins and cells. " << endl;
}


