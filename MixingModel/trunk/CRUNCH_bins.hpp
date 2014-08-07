// CRUNCH_bins.hpp
// header file for the CRUNCH binning object
// it is responsible for aggregating geochemical data
// across multiple columns and also includes printing 
// capabilities

#include <iostream>
#include <vector>
#include <map>
#include "CRN_tParticle_bins.hpp"
#include "VolumeParticleInfo.hpp"
#include "CRN_tParticle_bins.hpp"
#include "flowtube.hpp"
using namespace std;

#ifndef CRUNCH_bins_H
#define CRUNCH_bins_H

class CRUNCH_bins
{
	public:
	  /// @brief the default constructor which doens't do anything
	  CRUNCH_bins()			{ create(); }
	  
	  /// @brief constructor just initialises the volume particle info
	  /// object, and the number of cells
	  /// @author SMM
	  /// @date 07/08/2014
	  CRUNCH_bins( int n_pdz_per_bin, int n_caz_per_bin, double b_depth, 
                                 VolumeParticleInfo start_vpi)
	                    { create(n_pdz_per_bin, n_caz_per_bin, b_depth, start_vpi); }

    /// @brief the data is stored in vectors, you need to get the index into
    /// the vector of from the bin and cell number
    /// @author SMM
    /// @date 07/08/2014
    int retrieve_index_of_cell(int bin, int cell);
    
    /// @brief this gathers all the information from the CRN tParticle bin
    /// object including the locations of the cell boundaries, 
    /// and the relevant geochemical information
    /// @author SMM
    /// @date 07/08/2014
    void populate_cells_with_geochemical_data_from_CRNtPb(flowtube& ft,
                                  CRN_tParticle_bins& CRN_tPb);
        

  protected:
  
    /// The number of bins
    int n_bins;
    
    /// The number of particle types
    int n_types;
    
    /// the volume particle info
    VolumeParticleInfo vpi;
    
    /// The number of PDZ cells
    int n_pdz_cells_per_bin;
    
    /// The number of CAZ cells
    int n_caz_cells_per_bin;  
    
    /// The total number of cells
    int total_cells;
    
    /// The bottom depth of the cells you want
    double bottom_depth;
    
    /// This stores all the coordinates of the cell edges
    vector<double> cell_corners;
    
    /// This map stores all the data. Each vector in the map 
    /// is associated with a specific data type. This data can be found
    /// using a series of strings 
    map< string, vector<double> > cell_data_map; 
    
    /// this map stores indices, such as the node corners
    map< string, vector<int> > cell_index_map; 
    
       

	private:
	  void create();
	  
	  void create(int n_pdz_per_bin, int n_caz_per_bin, 
                double bottom_depth, VolumeParticleInfo this_vpi);
  
};

#endif
