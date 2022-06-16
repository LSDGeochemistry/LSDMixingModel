// Part_util

#ifndef Part_UTIL_H
#define Part_UTIL_H

#include "flowtube.hpp"
#include "CRN_tParticle_bins.hpp"
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>

using namespace std;

void part_ft_steady_erate(double SS_erate, double end_time,
					double ds_elev,
					vector<double> surf_erate,
					double insert_interval_time,
					double particle_print_interval,
					double zhe_print_interval,
					double age_print_interval,
					double eroded_catch_window,
					double max_age,
					int n_spacings,
					const int CRN_switch,
					string vtk_particle_fname,
					string vtk_cell_fname,
					string particle_list_fname,
					string particle_mfracs_fname,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,
					string CRN_parameter_fname,
					string ft_initial_profile_fname,
					string profile_out_fname,
					string particle_out_fname,
					string eroded_particle_out_fname,
					string age_pdf_out_fname,
					string age_cdf_out_fname,
					string zeta_out_fname, string eta_out_fname, string h_out_fname);

void part_ft_steady_flux(double SS_flux, double end_time,
					vector<double> surf_erate,
					double insert_interval_time,
					double particle_print_interval,
					double zhe_print_interval,
					double age_print_interval,
					double eroded_catch_window,
					double max_age,
					int n_spacings,
					const int CRN_switch,
					string vtk_particle_fname,
					string vtk_cell_fname,
					string particle_list_fname,
					string particle_mfracs_fname,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,
					string CRN_parameter_fname,
					string ft_initial_profile_fname,
					string profile_out_fname,
					string particle_out_fname,
					string eroded_particle_out_fname,
					string age_pdf_out_fname,
					string age_cdf_out_fname,
					string zeta_out_fname, string eta_out_fname, string h_out_fname);

void part_ft_erate_from_erh(string run_name);

void part_ft_erate_from_erh(string run_name, vector<double>& sample_s_locs, vector<double>& d_top, vector<double>& d_bottom,
							vector<double>& meas, vector<double>& unc, vector<double>& modelled, double& MLE,
							vector<double>& c_frac, vector<double>& pf1, vector<double>& pf2,
							vector<double>& pf3, vector<double>& pf4, vector<double>& pf5, vector<double>& C10Be);

void part_ft_steady_flux_schaller(double SS_flux, double end_time,
					vector<double> surf_erate,
					double insert_interval_time,
					double particle_print_interval,
					double zhe_print_interval,
					double age_print_interval,
					double eroded_catch_window,
					double max_age,
					int n_spacings,
					const int CRN_switch,
					string vtk_particle_fname,
					string vtk_cell_fname,
					string particle_list_fname,
					string particle_mfracs_fname,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,
					string CRN_parameter_fname,
					string ft_initial_profile_fname,
					string profile_out_fname,
					string particle_out_fname,
					string eroded_particle_out_fname,
					string age_pdf_out_fname,
					string age_cdf_out_fname,
					string zeta_out_fname, string eta_out_fname, string h_out_fname);

vector<double> part_ft_moraine_fit(string run_name, vector<double>& d_top_insitu, vector<double>& d_bottom_insitu,
							vector<double>& meas_insitu, vector<double>& modelled_insitu, double& MLE_insitu,
							vector<double>& d_top_meteor, vector<double>& d_bottom_meteor,
							vector<double>& meas_meteor, vector<double>& modelled_meteor, double& MLE_meteor);

void modify_parameter_files(string run_name, double start_depth, double vert_mix_vel, double part_conc,
					double single_scaling, double M_supply_surface,
					double k_f10Be, double deltad, double k2_f10Be, double chi_f10Be,
					double W_0, double gamma);

void age_calculations_bulk(double t_ime, double max_age, int n_spacings,
							double K_times_D, double D, double sigma,
						  flowtube& ft,
						  CRN_tParticle_bins& CRN_tpb,
						  vector< list<LSDCRNParticle> >& particle_bins,
						  ofstream& cdf_out,
						  ofstream& pdf_out);

void age_calculations_bins(double t_ime, double max_age, int n_spacings,
						double K_times_D, double D, double sigma,
						  flowtube& ft,
						  CRN_tParticle_bins& CRN_tpb,
						  vector< list<LSDCRNParticle> >& particle_bins,
						  ofstream& cdf_out,
						  ofstream& pdf_out);

void age_calculations_eroded_bulk(double t_ime, double max_age, int n_spacings,
							double K_times_D, double D, double sigma,
						  vector< list<CRN_tParticle> >& eroded_particle_bins,
						  ofstream& cdf_out,
						  ofstream& pdf_out);

void age_calculations_eroded_bins(double t_ime, double max_age, int n_spacings,
						double K_times_D, double D, double sigma,
						  vector< list<LSDCRNParticle> >& particle_bins,
						  ofstream& cdf_out,
						  ofstream& pdf_out);

double fraction_remaining(double pAge, double K_times_D, double D, double sigma);

#endif
