//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// CRUNCH_bins.hpp
// Object responsible for aggregating geochemical data
// across multiple columns and also includes printing 
// capabilities
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group mixing model
//  for exploring hillslope mixing and particle weathering
//
// Developed by:
//  Simon M. Mudd
//
// Copyright (C) 2018 Simon M. Mudd 2018
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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

    /// @brief this gets the mass fractions and the depletions and collects 
    /// the data into cells
    /// @author SMM
    /// @date 14/08/2014
    void populate_cells_with_mfracs_and_depletion(CRN_tParticle_bins& CRN_tPb);

    /// @brief This loops through the bins generating a CrunchFlow .in
    /// file for each bin
    /// @author SMM
    /// @date 10/08/2014
    void generate_CRUNCH_in_files(CRUNCH_engine& Ceng, flowtube& ft, 
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
 
    /// @brief this takes raw data and calcualtes the fraction from each mineral
    /// for example if the cell data has mass in kg in each cell, this will 
    /// print the mass fraction in each cell
    /// @author SMM
    /// @date 13/08/2014
    void vtk_print_cell_mineral_fractions(ofstream& vtk_cell_out, 
                    map< string, vector<double> >& raw_data);
                                       
    /// @brief This function prints the vtk files for the solid state
    /// properties of the minerals. Before you run this you need to 
    /// print the vtk_print_cell_header function (although this
    /// will be bundled in a wrapper later)
    /// @author SMM
    /// @date 08/08/2014
    void vtk_print_cell_mineral_solid_state(ofstream& vtk_cell_out);                                  

    /// @brief this function prints the mass fractions and the mineral
    /// depletion ratios to vtk
    /// @author SMM
    /// @date 14/08/2014
    void vtk_print_cell_mineral_mfracs_and_depletion(ofstream& vtk_cell_out);

    
    /// @brief This function prints the concentrations and reaction rates
    /// to a vtk file after crunch has run
    /// @author SMM
    /// @date 11/08/2014
    void vtk_print_cell_CRUNCH_data(ofstream& vtk_cell_out, 
                     CRUNCH_engine& Ceng);


    /// @brief This function bundles several of the vtk printing functions so
    /// you get the correct header, an appropriate filename, and both solute and
    /// mineral data
    /// @author SMM
    /// @date 13/08/2014
    void vtk_cell_bundler(double t_ime, int reference_switch, 
                                   string vtk_cell_fname, CRUNCH_engine& Ceng, 
                                   CRN_tParticle_bins& CRN_tPb); 
    
    
    /// @brief this function takes vectors of list vecs and reorganises them
    /// into vectors so that it is easier to plot the data to cells
    /// it also places them in a map container so the key to the map
    /// is used to name the variable in the vtk files
    /// @author SMM
    /// @date 08/08/2014
    map< string,vector<double> > parse_vec_list_vec_to_vec_map(string master_name, 
                            list<string> element_list,
                            vector< list < vector<double> > >& vlv);

    /// @brief This function gets data from a CRUNCH formatted list vec. 
    /// This style of list vec only has data fro the cells in the bin. 
    /// @author SMM
    /// @date 10/08/2014
    map< string, vector<double> > parse_CRUNCH_vec_list_vec_to_vec_map(string master_name, 
                            list<string> element_list,
                            vector< list < vector<double> > >& vlv);

    /// @brief this function call CRUNCH and then parses the data. It goes into
    /// vector list vectors, which can be read by the vtk routine
    /// @author SMM
    /// @date 11/08/2014
    void call_CRUNCH_and_parse_data(CRUNCH_engine& Ceng);
 
    /// @brief this weathers the particles based on the elements in the 
    /// vec list vecs that chave come from previous mineral properties
    /// or from CRUNCHFLOW
    /// @author SMM
    /// @date 11/08/2014
    void weather_particles_based_on_CRUNCH(CRN_tParticle_bins& CRN_tPb);
 
    /// @brief This bundles several CRUNCH_bins functions into one
    /// function that does a timestep of CRUNCH weathering
    /// @details IMPORTANT: the timestep is set my the MasterCrunch file
    /// you can change this with the CRUNCH_engine interface
    /// @author SMM
    /// @date 12/08/2014
    void run_CRUNCH_timestep(CRN_tParticle_bins& CRN_tPb, 
                  flowtube& ft, CRUNCH_engine& Ceng); 
                            
    /// @brief gets a list of the minerals as a list<string>. This is then used
    /// for printing the parameters to vtk files
    /// @author SMM
    /// @date 08/08/2014
    list<string> get_names_of_minerals();     
    
    /// @brief this function gets the names of the primary species in the 
    /// CRUNCH simulation. These are taken from the CRUNCH_engine object
    /// @author SMM
    /// @date 10/08/2014
    list<string> get_names_of_primary_species(CRUNCH_engine& cEng);                    
    
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
    vector< list< vector<double> > > vec_mineral_vpercents_old;
    vector< list< vector<double> > > vec_mineral_ssa_old;
    vector< list< vector<double> > > vec_mineral_mass_old;
    vector< list< vector<double> > > vec_mineral_surface_area_old;  
    vector< list< vector<double> > > vec_mineral_mfracs_old;  
    vector< list< vector<double> > > vec_mineral_depletion_old;        
          

    // vecvec to hold details about columns (that is, these are not mineral
    // or species specific)
    vector< vector<double> > vec_spacings;
    vector< vector<double> > vec_CRUNCH_tdepths;
    vector< vector<double> > vec_CRUNCH_bdepths;
    vector< vector<double> > vec_pH_vec;
    
    // these are vlvs to hold information from crunch
    vector< list< vector<double> > > vec_new_conc;
    vector< list< vector<double> > > vec_mineral_vpercents_new;
    vector< list< vector<double> > > vec_new_min_ssa;
    vector< list< vector<double> > > vec_new_rxn_rates;


  private:
    void create();
    
    void create(int n_pdz_per_bin, int n_caz_per_bin, 
                double bottom_depth, VolumeParticleInfo this_vpi);
  
};

#endif
