// CRUNCH_bins.cpp
// implementation of the CRUNCH_bins object
// it is responsible for aggregating geochemical data
// across multiple columns and also includes printing 
// capabilities


#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include "mathutil.hpp"
#include "flowtube.hpp"
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
void CRUNCH_bins::create(int n_pdz_per_bin, int n_caz_per_bin, 
                 double b_depth, VolumeParticleInfo this_vpi)
{
  n_pdz_cells_per_bin = n_pdz_per_bin;
  n_caz_cells_per_bin = n_caz_per_bin;
  bottom_depth = b_depth;
  vpi = this_vpi;
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
  int index_into_cell = (n_caz_cells_per_bin+n_pdz_cells_per_bin)*bin+cell;

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function:
// this just has some information about the cell locations
// the data is stored in vectors, you need to get the index into
// the vector of from the bin and cell number
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::populate_cells_with_geochemical_data_from_CRNtPb(flowtube& ft,
                                                CRN_tParticle_bins& CRN_tPb)
{


  // get the number of bins from the CRN_tParticle_bin_object
  n_bins = CRN_tPb.get_n_bins();
  
  // calculate the total number of cells 
  total_cells = (n_caz_cells_per_bin+n_pdz_cells_per_bin)*n_bins;
  
  // the bottom and top locations of all the cells
  // the reason why there is not one master vector for all the bins
  // is because PDZ thickness can change spatially. 
  vector<double> d_bottom_locs(total_cells);
  vector<double> d_top_locs(total_cells);   
  
  // the bottom and top locations of individual bins
  // these vectors will get replaced in each bin and are temporary
  vector<double> bin_d_bottom_locs;
  vector<double> bin_d_top_locs;

  // some temporary vectors for holding cell information
  vector<double> verts_s; 
  vector<double> verts_z;
  vector<double> verts_d;
  
  // some temporary vectors for holding cell indices
  vector<int> cell_node1;
  vector<int> cell_node2;
  vector<int> cell_node3;
  vector<int> cell_node4;
  
  // First partition the data into cells
  CRN_tPb.update_particles_cell_index(ft,
								n_pdz_cells_per_bin, n_caz_cells_per_bin,
								bottom_depth,
								verts_s, verts_z, verts_d,
							  cell_node1, cell_node2, cell_node3, cell_node4); 			

   cell_data_map["verts_s"] = verts_s;
   cell_data_map["verts_z"] = verts_z;
   cell_data_map["verts_d"] = verts_d;
   
   cell_index_map["cell_node1"] = cell_node1;
   cell_index_map["cell_node2"] = cell_node2;
   cell_index_map["cell_node3"] = cell_node3;
   cell_index_map["cell_node4"] = cell_node4;


/* 
  // get the locations of the cells and the cell indices. 
  // this function also updates the cell indices of the particles
  for (int bn = 0; bn<n_bins; bn++)
  {
    // first partition the bins into cells. This just gives the top and bottom
    // depths in the centrepoint of each cell. 
    CRN_tPb.partition_bins_into_cells(bn, ft, n_PDZ_cells_per_bin, n_CAZ_cells_per_bin,
										bottom_depth,bin_d_top_locs,bin_d_bottom_locs);
										
					
  
    // now collect the data from from the particles
    // this collects data from each cell in this bin
	  list< vector<double> > mineral_vfracs_old;
 	  list< vector<double> > mineral_ssa_old;
 	  list< vector<double> > mineral_mass_old;
 	  list< vector<double> > mineral_surface_area_old;
	  CRN_tPb.get_data_by_cell_volumetric_for_CRUNCH(bn,n_PDZ_cells_per_bin, 
                n_CAZ_cells_per_bin,bottom_depth,
								verts_s, verts_z, verts_d,
								cell_node1, cell_node2, cell_node3, cell_node4, vpi,
								mineral_vfracs_old,mineral_ssa_old,
								mineral_surface_area_old,mineral_mass_old);    
								
    // now you need to distribute data to the main vectors								
  
  
  }
*/  
  
  
}                                                
