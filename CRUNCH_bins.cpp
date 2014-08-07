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
#include "VolumeParticleInfo.hpp"
#include "CRN_tParticle_bins.hpp"
using namespace std;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::create()
{
  cout << "You need to at least provide the number of bins and cells. " << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function:
// this just has some information about the cell locations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::create(int n_pdz, int n_caz, VolumeParticleInfo this_vpi)
{
  n_pdz_cells = n_pdz;
  n_caz_cells = n_caz;
  vpi = ths_vpi;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function:
// this just has some information about the cell locations
// the data is stored in vectors, you need to get the index into
// the vector of from the bin and cell number
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int CRUNCH_bins::retrieve_index_of_cell(int bin, int cell)
{
  int index_into_cell = (n_caz_cells+n_pdz_cell)*bin+cell;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function:
// this just has some information about the cell locations
// the data is stored in vectors, you need to get the index into
// the vector of from the bin and cell number
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::populate_cells_with_geochemical_data_from_CRNtPb(
                                                CRN_tParticle_bins& CRN_tPb);
{


  // get the number of bins
  

}                                                
