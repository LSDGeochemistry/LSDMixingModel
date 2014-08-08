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
#include "CRUNCH_engine.hpp"
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

    /// @brief This function prints the corners of the cell to screen
    /// @author SMM
    /// @date 07/08/2014
    void cell_location_to_screen(int bin, int cell);

    /// @brief this gathers all the information from the CRN tParticle bin
    /// object including the locations of the cell boundaries, 
    /// and the relevant geochemical information
    /// @author SMM
    /// @date 07/08/2014
    void populate_cells_with_geochemical_data_from_CRNtPb(flowtube& ft,
                                  CRN_tParticle_bins& CRN_tPb);
        
    /// @brief this function prints data members to vtk
    /// @author SMM
    /// @date 08/08/2014
    void vtk_print_cell_header(int reference_frame_switch, ofstream& vtk_cell_out);
 
    /// @brief this takes a map of vectors and prints them to a VTK cell
    /// file
    /// @author SMM
    /// @date 08/08/2014
    void vtk_print_cell_from_map_of_vectors(ofstream& vtk_cell_out, 
                                       map<string, vector<double> >& data_map);
                                       
    /// @brief This function prints the vtk files for the solid state
    /// properties of the minerals. Before you run this you need to 
    /// print the vtk_print_cell_header function (although this
    /// will be bundled in a wrapper later)
    /// @author SMM
    /// @date 08/08/2014
    void vtk_print_cell_mineral_solid_state(ofstream& vtk_cell_out);                                  
    
    
    /// @brief this function takes vectors of list vecs and reorganises them
    /// into vectors so that it is easier to plot the data to cells
    /// it also places them in a map container so the key to the map
    /// is used to name the variable in the vtk files
    /// @authors SMM
    /// @date 08/08/2014
    map< string,vector<double> > parse_vec_list_vec_to_vec_map(string master_name, 
                            list<string> element_list,
                            vector< list < vector<double> > >& vlv);
                            
    /// @brief this function takes a VolumeParticleInfo file and 
    /// gets a list of the minerals as a list<string>. This is then used
    /// for printing the parameters to vtk files
    list<string> get_names_of_minerals();                        
    
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
    
    // vec list vecs to hold mineral information
	  vector< list< vector<double> > > vec_mineral_vfracs_old;
 	  vector< list< vector<double> > > vec_mineral_ssa_old;
 	  vector< list< vector<double> > > vec_mineral_mass_old;
 	  vector< list< vector<double> > > vec_mineral_surface_area_old;       

    // vecvec to hold details about columns (that is, these are not mineral
    // or species specific)
	  vector< vector<double> > vec_spacings;
	  vector< vector<double> > vec_CRUNCH_tdepths;
	  vector< vector<double> > vec_CRUNCH_bdepths;
	  vector< vector<double> > vec_pH_vec;
	  
	  // these are vlvs to hold information from crunch
    vector< list< vector<double> > > vec_new_conc;
	  vector< list< vector<double> > > vec_mineral_vfracs_new;
	  vector< list< vector<double> > > vec_new_min_ssa;
	  vector< list< vector<double> > > vec_new_rxn_rates;


	private:
	  void create();
	  
	  void create(int n_pdz_per_bin, int n_caz_per_bin, 
                double bottom_depth, VolumeParticleInfo this_vpi);
  
};

#endif
