//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// flowtub
// An object that controls a 1-D hillslope, keeps track of hillslope
// evolution, sediment transport and interfaces with particles
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
#include <fstream>
#include "tParticle.hpp"
//#include "CRN_tParticle_bins.hpp"
using namespace std;

#ifndef flowtube_H
#define flowtube_H

/// This object controls sediment flux along a "flowtube". This is a pseudo 1-D
/// hillslope profile where the width of the profile changes in order
/// to account for divergence or convergence of the terrain
class flowtube
{
  public:
  
    /// @brief default constructor. Doesn't do anything
    /// @author SMM
    /// @date 01/01/2008
    flowtube()				{ create(); }
    
    /// @brief Constructor that loads in input file
    ///  *add input file format here*
    /// @param infile an ifstream object
    /// @author SMM
    /// @date 01/01/2008
    flowtube(ifstream& infile)         { create(infile); }
    
    /// @brief A constructor that assigns all the parameters (used for copying)
    /// @author SMM
    /// @date 01/01/2008
    flowtube(int sn_nodes, double sS_c, double sK_h, double sW_0, double sgamma,
             double srho_s, double srho_r, double sN, double sN_0, double sN_m,
             double sbeta, double sK_g, vector<double> sA, vector<double> sA_bins,
             vector<double> sb, vector<double> szeta, vector<double> seta,
             vector<double> sh, vector<double> ss_h, vector<double> ss_b,
             vector<double> sDeltaXh, vector<double> sDeltaXb,
             vector<double> sbin_edge_loc,vector<double> sold_h,
             vector<double> sold_eta, vector<double> sold_zeta,
             vector<double> sintermediate_zeta, vector<double> sintermediate_h,
             vector<double> spre_surface_zeta, vector<double> spre_surface_h,
             vector<double> sMass_Flux, vector<double> sfluff)
             { create(sn_nodes, sS_c, sK_h, sW_0, sgamma,
                      srho_s, srho_r, sN, sN_0, sN_m,
                      sbeta, sK_g, sA, sA_bins, sb,szeta, seta,
                      sh, ss_h,  ss_b, sDeltaXh,  sDeltaXb,
                      sbin_edge_loc, sold_h, sold_eta, sold_zeta,
                      sintermediate_zeta, sintermediate_h,
                      spre_surface_zeta, spre_surface_h,
                      sMass_Flux, sfluff); }

    /// @brief The copy constructor
    /// @author SMM
    /// @date 01/01/2008
    flowtube& operator=(flowtube& ft);

    /// @brief Prints flowtube properties to file
    /// @param outfile An ofstream object
    /// @author SMM
    /// @date 01/01/2008
    void print_ft_properties(ofstream& outfile);
  
    void initialize_interpolation_nodes(vector<double> interp_x_loc,
                                        vector<int>& ds_interp_node_num,
                                        vector<int>& us_interp_node_num,
                                        vector<double>& interpolated_fractional_distance);
            
    vector<double> interpolated_modeled_h(vector<int> ds_interp_node_num,
                                          vector<int> us_interp_node_num,
                                          vector<double> interpolated_fractional_distance );
            
    vector<double> interpolated_modeled_zeta(vector<int> ds_interp_node_num,
                                             vector<int> us_interp_node_num,
                                             vector<double> interpolated_fractional_distance );
    
    /// @brief This loads a profile. The file has five columns but only the last two
    ///   are read. These are the zeta (elevation of surface) and eta (elevation of saprolite surface)
    ///   The soil thickness is calculated from these data. The data are reported at the nodes and 
    ///   not the node boundaries 
    /// @author SMM
    /// @date 01/01/2011
    void load_profile(ifstream& infile);

    /// @brief Exports a profile. It has data at centre nodes but not boundaries
    ///  It has 5 columns: 
    ///   distance, elevation, soil thickness, elevation, and soil thickness
    ///  yes the elevation and soil thicknesses are repeated. I can't remember
    ///  why I did that. I think maybe before there were old and new thicknesses
    ///  and elevations. 
    /// @author SMM
    /// @date 01/01/2011    
    void export_input_profile(ofstream& outfile);
    
    /// @brief This raises all elevation points relative to the downslope boundary
    /// @param ds_elev The elevation at the downslope boundary (in m) 
    /// @author SMM
    /// @date 01/01/2011
    void raise_zeta_eta_ds_bound(double ds_elev);
    
    /// @brief This raises all elevation points relative to the mean elevation
    /// @param mean_elev The mean elevation of the hillslope (in m)
    /// @author SMM
    /// @date 01/01/2011    
    void raise_zeta_eta_mean(double mean_elev);
    
    /// @brief This creates a flat surface with constant soil thickness
    /// @param zeta_flat The elevation of the surface (in m)
    /// @param h_flat The constant soil thickness (in m)
    /// @author SMM
    /// @date 01/01/2011   
    void set_const_zeta_eta_h(double zeta_flat, double h_flat);

    /// @brief Runs a timestep with the boundary condition to be a set flux
    /// @paramater dt The timestep (in years)
    /// @parameter flux_us The flux from upslope in kg/yr
    /// @parameter flux_ds The flux at the downslope boundary in kg/yr
    /// @parameter flux_switch A switch for the type of flux (need to look up)
    /// @parameter prod_switch A switch for the trpe of production function
    /// @parameter surface_change_rate A vector holding the change in surface elevation
    ///      in m/yr to simulate, for exampe erosion from overland flow
    /// @author SMM
    /// @date 01/01/2011    
    void flux_timestep_flux_bc(double dt,
							double flux_us, double flux_ds,
							int flux_switch, int prod_switch,
							vector<double> surface_change_rate);
    
    /// @brief Runs a timestep with the boundary condition to be a set elevation
    /// @paramater dt The timestep (in years)
    /// @parameter flux_us The flux from upslope in kg/yr
    /// @parameter ds_elevation the elevation of the downslope boundary (in m)
    /// @parameter flux_switch A switch for the type of flux (need to look up)
    /// @parameter prod_switch A switch for the trpe of production function
    /// @parameter surface_change_rate A vector holding the change in surface elevation
    ///      in m/yr to simulate, for exampe erosion from overland flow
    /// @author SMM
    /// @date 01/01/2011
    void flux_timestep_elev_bc(double dt,
							double flux_us, double ds_elevation,
							int flux_switch, int prod_switch,
							vector<double> surface_change_rate);

    // Getter functions
    int get_n_nodes()     { return n_nodes; }
    double get_S_c()      { return S_c; }
    double get_K_h()      { return K_h; }
    double get_W_0()      { return W_0; }
    double get_gamma()    { return gamma; }
    double get_rho_s()    { return rho_s; }
    double get_rho_r()    { return rho_r; }
    double get_N()        { return N; }
    double get_N_0()      { return N_0; }
    double get_N_m()      { return N_m; }
    double get_beta()     { return beta; }
    double get_K_g()      { return K_g; }

    vector<double> get_A()		{ return A; }
    vector<double> get_A_bins()	{ return A_bins; }
    vector<double> get_b()		{ return b; }
    vector<double> get_zeta()	{ return zeta; }
    vector<double> get_eta()	{ return eta; }
    vector<double> get_h()		{ return h; }
    vector<double> get_s_h()	{ return s_h; }
    vector<double> get_s_b()	{ return s_b; }
    vector<double> get_DeltaXh()
								{ return DeltaXh; }
	vector<double> get_DeltaXb()
								{ return DeltaXb; }
	vector<double> get_bin_edge_loc()
								{ return bin_edge_loc; }


    // Getter functions, these get old data members
    vector<double> get_old_h()	{ return old_h; }

    vector<double> get_old_zeta()
								{ return old_zeta; }
    vector<double> get_old_eta()
								{ return old_eta; }
    vector<double> get_intermediate_zeta()
								{ return intermediate_zeta; }
    vector<double> get_intermediate_h()
								{ return intermediate_h; }
    vector<double> get_pre_surface_zeta()
								{ return pre_surface_zeta; }
    vector<double> get_pre_surface_h()
								{ return pre_surface_h; }
    vector<double> get_Mass_Flux()
								{ return Mass_Flux; }
    vector<double> get_fluff()	{ return fluff; }

    /// This gets the ds_elevation for the elevation boundary condition
    double get_zeta_ds()		{ return zeta[n_nodes-1]; }
	
    
    ///20/09/18 LK
    ///need to update for what N, N_0, N_m and K_g are for using cases 5, 6 and 7 in the flux laws
    ///beta describes the rate of exponential decline with depth? (see Roering 2008 equations 10 through 16).
    ///The rest I'm not sure of but will need to make sure they're read in the sed_trans_param file to use these.
    
    
    /// @brief This sets the transport parameters for the soil transport 
    /// @param temp_S_c Critical slope (dimensionless)
    /// @param temp_K_h The sediment transport coefficient (units depend on flux law)
    /// @param temp_W_0 Soil production at zero thickness in m/yr
    /// @paarm temp_gamma E folding distance of soil production (m)
    /// @param temp_rho_s Soil density (kg/m^3)
    /// @param temp_rho_r Rock density (kg/m^3)
    /// @author SMM
    /// @date 01/01/2011
    void set_transport_params(double temp_S_c, double temp_K_h,
								double temp_W_0, double temp_gamma,
							  	double temp_rho_s, double temp_rho_r);
    
	void set_transport_params(double temp_S_c, double temp_K_h,
								double W_0, double temp_gamma,
							  	double temp_rho_s, double temp_rho_r,
							  	double temp_N, double temp_N_0,
							  	double temp_N_m,
							  	double temp_beta, double temp_K_g);
	double calculate_steady_flux(double steady_erosion_rate, double flux_us);

	void print_zeta(double t_ime, ofstream& zeta_out);
	void print_relative_zeta(double t_ime, ofstream& zeta_out);
	void print_eta(double t_ime, ofstream& eta_out);
	void print_h(double t_ime, ofstream& h_out);
	void print_ft_vecs_to_screen();
	void print_h_s(ofstream& outfile);

	protected:
	int n_nodes;			// number of nodes in the flowtube
	double S_c;				// critical slope (dimensionless)
	double K_h;				// diffusivity
	double W_0;
	double gamma;
	double rho_s;
	double rho_r;
	double N;
	double N_0;
	double N_m;
	double beta;
	double K_g;
	vector<double> A;		// area (in horizontal plane) of the nodes (m^2)
	vector<double> A_bins;	// area of the particle bins
	vector<double> b;		// width of the flowtube (m)
	vector<double> zeta;	// surface elevation (m)
	vector<double> eta;		// elevation of soil-bedrock boundary (m)
	vector<double> h;		// soil thickness (m)
	vector<double> s_h;		// distance downslope (m)
	vector<double> s_b;		// distance downslope (m)
	vector<double> DeltaXh;	// node spacing (m)
	vector<double> DeltaXb;	// node spacing (m)
	vector<double> bin_edge_loc;
							// location of bins for particle tracking

	// a number of vectors store information from timestapes and are used
	// for particle tracking.
	vector<double> old_h;
	vector<double> old_zeta;
	vector<double> old_eta;
	vector<double> intermediate_zeta;
	vector<double> intermediate_h;
	vector<double> pre_surface_zeta;
	vector<double> pre_surface_h;
	vector<double> Mass_Flux;
	vector<double> fluff;

	private:
	void create();
	void create(ifstream& infile);
	void create(int sn_nodes, double sS_c, double sK_h, double sW_0, double sgamma,
			 double srho_s, double srho_r, double sN, double sN_0, double sN_m,
			 double sbeta, double sK_g, vector<double> sA, vector<double> sA_bins,
			 vector<double> sb, vector<double> szeta, vector<double> seta,
			 vector<double> sh, vector<double> ss_h, vector<double> ss_b,
			 vector<double> sDeltaXh, vector<double> sDeltaXb,
			 vector<double> sbin_edge_loc,vector<double> sold_h,
			 vector<double> sold_eta, vector<double> sold_zeta,
			 vector<double> sintermediate_zeta, vector<double> sintermediate_h,
			 vector<double> spre_surface_zeta, vector<double> spre_surface_h,
			 vector<double> sMass_Flux, vector<double> sfluff);
};





#endif


