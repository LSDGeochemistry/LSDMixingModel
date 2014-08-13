// CRN_tParticle_bins.h
//

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
	CRN_tParticle_bins()				{ create(); }
	CRN_tParticle_bins(flowtube ft)		{ create(ft); }
	CRN_tParticle_bins(vector<double> bin_edge_locations,vector<double> s_h,
						vector<double> s_b, vector<double> b, vector<double> A_bins)
								{ create(bin_edge_locations,s_h,s_b,b,A_bins); }

	// insert functions
	int insert_particles(flowtube ft,vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double part_conc, vector<int> starting_pID, vector<double> starting_p_mfrac);

	int insert_particles(flowtube ft, vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double part_conc, vector<int> starting_pID, vector<double> starting_p_mfrac,
									double C_10Be, double C_26Al, double C_36Cl, double C_14C, double C_21Ne, double C_3He);
	int insert_particles(flowtube ft, vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double part_conc, vector<int> starting_pID, vector<double> starting_p_mfrac,
									CRN_parameters& CRNp,
									double erosion_rate_in_mass_per_time_per_area);
	int insert_particles_volumetric(flowtube ft,
									vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double C_10Be, double C_26Al, double C_36Cl, double C_14C, double C_21Ne, double C_3He,
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
	void update_CRN_conc_const(double C_10Be, double C_26Al,
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
	void create(flowtube ft);


};





#endif
