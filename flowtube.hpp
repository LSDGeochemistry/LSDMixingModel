// flowtube.hpp
// the flowtube object

#include <iostream>
#include <vector>
#include <fstream>
#include "tParticle.hpp"
//#include "CRN_tParticle_bins.hpp"
using namespace std;

#ifndef flowtube_H
#define flowtube_H

class flowtube
{
	public:
	flowtube()				{ create(); }
	flowtube(ifstream& infile)
							{ create(infile); }
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

	flowtube& operator=(flowtube& ft);

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
	void load_profile(ifstream& infile);
	void export_input_profile(ofstream& outfile);
	void raise_zeta_eta_ds_bound(double ds_elev);
	void raise_zeta_eta_mean(double mean_elev);
	void set_const_zeta_eta_h(double zeta_flat, double h_flat);

	void flux_timestep_flux_bc(double dt,
							double flux_us, double flux_ds,
							int flux_switch, int prod_switch,
							vector<double> surface_change_rate);
	void flux_timestep_elev_bc(double dt,
							double flux_us, double ds_elevation,
							int flux_switch, int prod_switch,
							vector<double> surface_change_rate);

	int get_n_nodes() 			{ return n_nodes; }
	double get_S_c()			{ return S_c; }
	double get_K_h()			{ return K_h; }
	double get_W_0()			{ return W_0; }
	double get_gamma()			{ return gamma; }
	double get_rho_s()			{ return rho_s; }
	double get_rho_r()			{ return rho_r; }
	double get_N()				{ return N; }
	double get_N_0()			{ return N_0; }
	double get_N_m()			{ return N_m; }
	double get_beta()			{ return beta; }
	double get_K_g()			{ return K_g; }

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

	double get_zeta_ds()		{ return zeta[n_nodes-1]; }
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


