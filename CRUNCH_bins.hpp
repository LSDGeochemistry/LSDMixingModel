// CRUNCH_bins.hpp
// header file for the CRUNCH binning object
// it is responsible for aggregating geochemical data
// across multiple columns and also includes printing 
// capabilities

#include <iostream>
#include <vector>
#include "CRN_tParticle_bins.hpp"
using namespace std;

#ifndef CRUNCH_bins_H
#define CRUNCH_bins_H

class CRUNCH_bins
{
	public:
	  CRUNCH_binss()			{ create(); }

  protected:
  
    /// The number of bins
    int n_bins;
    
    /// The number of PDZ cells
    int n_pdz_cells;
    
    /// The number of CAZ cells
    int n_caz_cells;
    
    /// This stores all the coordinates of the cell edges
    vector<double> cell_corners;
    
    /// This map stores all the data. Each vector in the map 
    /// is associated with a specific data type. This data can be found
    /// using a series of strings 
    map< string, vector<double> > cell_data_map; 
   

	private:
	  void create();
  
};

#endif
