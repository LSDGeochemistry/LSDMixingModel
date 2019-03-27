//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// CRN_tParticle_bins
// Cosmogenic nuclide particle bins object
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

#ifndef CRN_tParticle_bins_H
#define CRN_tParticle_bins_H

#include <iostream>
#include <vector>
#include <list>
#include "tParticle.hpp"
#include "flowtube.hpp"
#include "chronos_particle_info.hpp"
#include "CRN_parameters.hpp"
#include "VolumeParticleInfo.hpp"
using namespace std;




class CRN_tParticle_bins
{
  public:
    
    /// @brief the default constructor which doens't do anything
    CRN_tParticle_bins()				{ create(); }
    
    /// @brief Construct a particle bins object from a flowtube
    /// @param ft A flowtube oject
    /// @author SMM
    /// @date 01/01/2011
    CRN_tParticle_bins(flowtube ft)		{ create(ft); }
    
    /// @brief Construct a particle bins object from raw data of a flowtube
    /// @param bin_edge_locations Vector with downslope locations of bin edges (in metres)
    /// @param s_h vector of downslope locations of bin centrepoints
    /// @param s_b vector of downslope locations of bin edges
    /// @param b vector of widths of each bin
    /// @param A_bins vector of basal area of each bin
    /// @author SMM
    /// @date 01/01/2011
    CRN_tParticle_bins(vector<double> bin_edge_locations,vector<double> s_h,
                       vector<double> s_b, vector<double> b, vector<double> A_bins)
                  { create(bin_edge_locations,s_h,s_b,b,A_bins); }

    // insert functions
    int insert_particles(flowtube ft,vector<double> Delta_lowered, vector<double>& old_bottom_depth,
                         double part_conc, vector<int> starting_pID, vector<double> starting_p_mfrac);

	int insert_particles(flowtube ft, vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double part_conc, vector<int> starting_pID, vector<double> starting_p_mfrac,
									double C_10Be, double C_f10Be, double C_26Al, double C_36Cl, double C_14C, double C_21Ne, double C_3He);
	int insert_particles(flowtube ft, vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double part_conc, vector<int> starting_pID, vector<double> starting_p_mfrac,
									CRN_parameters& CRNp,
									double erosion_rate_in_mass_per_time_per_area);
	int insert_particles_volumetric(flowtube ft,
									vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double C_10Be, double C_f10Be, double C_26Al, double C_36Cl, double C_14C, double C_21Ne, double C_3He,
									VolumeParticleInfo vpi);

    vector< list<CRN_tParticle> > particle_motion(double dt, flowtube ft,
										double Omega,double vert_vel_fluctuating,
										double horiz_vel_fluctuating,const int CRN_switch,
										CRN_parameters& CRNp);

	// fallout functions
	void update_fallout_10Be_bins(double dt, double M_supply_surface,
					double rho_s, double k_f10Be, double deltad, CRN_parameters& CRNp);
	void update_fallout_10Be_bins(double dt, double M_supply_surface,
					double rho_s, double k1_f10Be, double k2_f10Be, double chi_f10Be,
					double deltad, CRN_parameters& CRNp);

	// functions related to sampling
	int get_bin_number_of_sloc(double s_loc);
	void update_particles_cell_index(flowtube ft,
									int n_depthintervals_soil, int n_depthintervals_parent,
									double bottom_depth, vector<double>& verts_s, vector<double>& verts_z,
									vector<double>& verts_d,
									vector<int>& cell_node1, vector<int>& cell_node2,
									vector<int>& cell_node3, vector<int>& cell_node4);

	void check_particles_in_cells(int bn, vector<double>& verts_s, vector<double>& verts_z,
								vector<double>& verts_d,
								vector<int>& cell_node1, vector<int>& cell_node2,
								vector<int>& cell_node3, vector<int>& cell_node4);

	void get_data_by_cell(int bn, int n_PDZ_intervals, int n_CAZ_intervals,
									double bottom_depth, vector<double> verts_s, vector<double> verts_z,
									vector<double>& verts_d,
									vector<int>& cell_node1, vector<int>& cell_node2,
									vector<int>& cell_node3, vector<int>& cell_node4);
	
    /// @brief This function assumes the cells have already been identified for the
    /// particles. It then calcualted the mass fractions in each cell of the
    /// particle types, as well as the mass depletion, as measured by the
    /// total mass of a type in a cell divided by the total starting type of 
    /// mass in a cell
    /// @author SMM
    /// @date 14/08/2014								
    void get_mineral_mass_loss_and_mfracs_volumetric(int bn, int n_PDZ_intervals, 
                  int n_CAZ_intervals,
									VolumeParticleInfo& vpi,
									list< vector<double> >& mineral_mfracs,
									list< vector<double> >& mineral_depletion);								
									
	void get_data_by_cell_volumetric_for_CRUNCH(int bn, int n_PDZ_intervals, int n_CAZ_intervals,
									double bottom_depth, vector<double> verts_s, vector<double> verts_z,
									vector<double>& verts_d,
									vector<int>& cell_node1, vector<int>& cell_node2,
									vector<int>& cell_node3, vector<int>& cell_node4,
									VolumeParticleInfo vpi,
									list< vector<double> >& mineral_vpercents,
									list< vector<double> >& mineral_ssa,
									list< vector<double> >& mineral_surface_area,
									list< vector<double> >& mineral_mass);
	list<CRN_tParticle> get_particles_in_depth_interval(int bn, double d_top, double d_bottom);

	void calculate_sample_averages(int bn, vector<double>& d_top_locs, vector<double>& d_bottom_locs,
													vector<double>& sample_mean_age,
													vector<double>& sample_mean_C10Be,
													vector<double>& sample_mean_Cf10Be);
	void calculate_sample_averages(vector<double>& s_locs,
								vector<double>& d_top_locs, vector<double>& d_bottom_locs,
								Particle_info& pi,
								vector<double>& sample_mean_age, vector<double>& sample_mean_C10Be,
								vector<double>& sample_enrich, vector<double>& sample_frac_clay,
								vector<double>& sample_frac_t1, vector<double>& sample_frac_t2,
								vector<double>& sample_frac_t3, vector<double>& sample_frac_t4,
								vector<double>& sample_frac_t5);

	void weather_particles_from_CRUNCH(int bn, int n_PDZ_intervals, int n_CAZ_intervals,
									double bottom_depth, vector<double> verts_s, vector<double> verts_z,
									vector<double>& verts_d,
									vector<int>& cell_node1, vector<int>& cell_node2,
									vector<int>& cell_node3, vector<int>& cell_node4,
									VolumeParticleInfo vpi,
									list< vector<double> >& mineral_vpercents_old,
                                    list< vector<double> >& mineral_vpercents_new,
                                    list< vector<double> >& surface_area_in_cell_old,
                                    list< vector<double> >& mass_in_cell_old);

	// this function takes a PDZ and CAZ depth interval: it seperates the soil into
	// n_PDZ_intervals and the CAZ into n_CAZ_intervals. The bottom of the CAZ interval
	// is depth_to_bottom
	void partition_bins_into_cells(int bn, flowtube ft,
							int n_PDZ_intervals, int n_CAZ_intervals,
									double depth_to_bottom,
								   vector<double>& d_top_locs,
								   vector<double>& d_bottom_locs);


	// functions for updating CRN concentration
	void update_CRN_conc_const(double C_10Be, double C_f10Be, double C_26Al,
										double C_36Cl, double C_14C,
										 double C_21Ne, double C_3He);

	// printing functions
	void print_particle_stats(double t_ime, ofstream& particle_out);
	void print_particle_stats(double t_ime, flowtube ft,ofstream& particle_out);
	void print_particle_stats_soil(double t_ime, flowtube ft,
								ofstream& particle_out);
	void print_eroded_stats(double t_ime,
							vector< list<CRN_tParticle> > eroded_list_vec,
							flowtube ft, ofstream& particle_out);
	void print_particle_stats_vtk(double t_ime, flowtube ft,
								string vtk_fname);

  /// @brief This function is a basic vtk printing routine for volume 
  /// particles
  /// @author SMM
  /// @date 13/08/2014
  void vtk_print_basic_volume_particles(double t_ime, 
								 string vtk_particle_fname, int reference_frame_switch);


	void print_age_cdfpdf_bulk(double t_ime, double max_age, int n_spacings,
							double K_times_D, double D, double sigma,
						  flowtube& ft,
						  ofstream& cdf_out,
						  ofstream& pdf_out);
	void print_age_cdfpdf_bins(double t_ime, double max_age, int n_spacings,
						  double K_times_D, double D, double sigma,
						  flowtube& ft,
						  ofstream& cdf_out,
						  ofstream& pdf_out);
	void print_age_cdfpdf_eroded_bulk(double t_ime, double max_age, int n_spacings,
						  vector< list<CRN_tParticle> >& eroded_particle_bins,
						  double K_times_D, double D, double sigma,
						  ofstream& cdf_out,
						  ofstream& pdf_out);
	void print_age_cdfpdf_eroded_bins(double t_ime, double max_age, int n_spacings,
						  double K_times_D, double D, double sigma,
						  vector< list<CRN_tParticle> >& particle_bins,
				     	  ofstream& cdf_out,
						  ofstream& pdf_out);

	void cell_printing_vtk(double t_ime, flowtube ft,string vtk_fname,
						   int n_depthintervals_soil, int n_depthintervals_parent,
						   double bottom_depth);

	void cell_and_particle_printing_vtk(double t_ime, flowtube ft,
								string vtk_particle_fname, string vtk_cell_fname,
								int n_depthintervals_soil, int n_depthintervals_parent,
								double bottom_depth);

	void cell_and_particle_chemistry_printing_vtk(double t_ime, flowtube ft,
									 	Particle_info& pi,
									 	string vtk_particle_fname, string vtk_cell_fname,
										int n_depthintervals_soil, int n_depthintervals_parent,
										double bottom_depth, int reference_frame_switch);


	int get_n_bins()			{ return n_bins; }
	vector< list<CRN_tParticle> > get_particle_bins()
								{ return particle_bins; }
	vector<double> get_dx_h()   { return dx_h; }
	vector<int>    get_h_node_us()
								{ return h_node_us; }
	vector<int>    get_h_node_ds()
								{ return h_node_ds; }
	vector<double> get_s_us_h()	{ return s_us_h; }
    void create(flowtube ft);
	protected:
	int n_bins;				// number of bins
	vector<double> bin_edge_loc;
	vector< list<CRN_tParticle> > particle_bins;

	vector<int> h_node_us;
	vector<int> h_node_ds;
	vector<int> b_node_us;			// these are node reference for the width of
	vector<int> b_node_ds;			// upslope and downslope locations of the
									// bin edges

	vector<double> dx_h;
	vector<double> dx_b;
	vector<double> s_us_h;
	vector<double> s_us_b;
	vector<double> b_us;
	vector<double> Slope_b;
	vector<double> bin_width;
	vector<double> A_bins;

	private:
	void create();
	void create(vector<double> bel,vector<double> s_h,
						vector<double> s_b, vector<double> b,
						vector<double> A_bin_temp);
	//void create(flowtube ft);


};





#endif
