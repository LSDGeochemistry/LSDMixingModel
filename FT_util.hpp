// FT_util

#ifndef FT_UTIL_H
#define FT_UTIL_H

#include "flowtube.hpp"
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>

using namespace std;


double RMSE(vector<double> interp_meas, vector<double> modeled);
double RMSE_mean(vector<double> interp_meas, vector<double> modeled);
double RMSE(vector<double> interp_meas, vector<double> modeled,
					vector<double> independant_var, double ind_var_start,
					double efold_weight_dist);
double RMSE_local(vector<double> interp_meas, vector<double> modeled,
			vector<double> independant_var, double ind_var_before_end);
void RMSE_zeta_h_twofiles(int n_runs,string zeta_name1, string zeta_name2,
					string h_name1, string h_name2,
					string steady_modeled_fname,
					string RMSE_compare_out_fname);
double fraction_in_window(double window_time,
					vector<double> time_vec,
					vector<double> erate_vec,
					vector<double> upper_window,
					vector<double> lower_window);
double load_erate_history(ifstream& er_in,vector<double>& erate_change_times,
					vector<double>& erate_change_rates);
double erosion_rate_from_erate_history_file(int& erate_change_index,
					double& erate_change_time,
					vector<double> erate_change_times,
					vector<double> erate_change_rates);
double erate_history_arbitrary_start(int& erate_change_index,
					double& erate_change_time,
					vector<double> erate_change_times,
					vector<double> erate_change_rates,
					double start_time);
void create_erate_history(long seed,ofstream& erh_out,
					double start_time, double start_rate,
					double end_time, double dt,
					double max_eros, double min_eros,
					double mean_change, double mean_interval,
                    vector<double>& erate_change_times,
					vector<double>& erate_change_rates);
void create_erate_history(long seed,
					double start_time, double start_rate,
					double end_time, double dt,
					double max_eros, double min_eros,
					double mean_change, double mean_interval,
                    vector<double>& erate_change_times,
					vector<double>& erate_change_rates);
void create_erate_hist_datfile(int start_runnumber, int end_runnumber,
					double start_time,
					double end_time, double max_eros,
					double min_eros, double mean_change,
					double mean_interval, double start_rate,
					string model_run_params_fname,
					string erf_datfile_fname);
void print_erate_hist_datfile(ofstream& erh_datfile,
					int run_number,
					double RMSE_zeta,double RMSE_h,
					vector<double> erate_change_times,
					vector<double> erate_change_rates);
void generate_erate_hist_from_datfile(string erf_datfile_fname,
					string erf_out_fname, int run_number);
void generate_erate_hist_from_datfile(ifstream& erf_datfile,
					ofstream& erf_out, int run_number,
					vector<double>& erate_change_times,
					vector<double>& erate_change_rates);
void generate_erate_hist_from_datfile(ifstream& erf_datfile,
					int run_number,
					vector<double>& erate_change_times,
					vector<double>& erate_change_rates);
void generate_erate_hist_from_datfile(string erf_datfile_name,
					int run_number,
					vector<double>& erate_change_times,
					vector<double>& erate_change_rates);
void generate_regularly_spaced_erate_hist(
					double t_spacing,
					int run_number,
					ifstream& erf_datfile,
					vector<double>& er_times,
					vector<double>& erate_hist_regular);
void compare_single_erate_hist_with_others(
					double t_spacing,
					int run_number_comparison_erh,
					int start_compare_rn,
					int end_compare_rn,
					double efold_weight_dist,
					string erh_fname,
					string RMSE_erate_fname);
void compare_single_erate_hist_with_others(
					double t_spacing,
					int run_number_comparison_erh,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> efold_weight_dist_vec,
					string erh_fname,
					string RMSE_erate_fname);
void compare_single_erate_hist_with_others(
					double t_spacing,
					vector<double> compare_reg_times,
					vector<double> compare_reg_rates,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> efold_weight_dist_vec,
					string erh_fname,
					string RMSE_erate_fname);
void compare_single_erate_hist_block_RMSE(
					double t_spacing,
					vector<double> compare_reg_times,
					vector<double> compare_reg_rates,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> last_block_vec,
					string erh_fname,
					string RMSE_erate_fname);
void compare_erate_windows(double t_spacing,
					int run_number_comparison_erh,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> window_times_vec,
					string window_2p_out_fname,
					string window_5p_out_fname,
					string window_10p_out_fname,
					string erh_fname);
void compare_erate_windows(double t_spacing,
					vector<double> compare_reg_times,
					vector<double> compare_reg_rates,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> window_times_vec,
					string window_2p_out_fname,
					string window_5p_out_fname,
					string window_10p_out_fname,
					string erh_fname);
void ft_steady_erate(double SS_erate, double end_time,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,
					string ft_initial_profile_fname,
					string profile_out_fname);
void ft_steady_erate(double SS_erate, double end_time,
					double ds_elev,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,
					string ft_initial_profile_fname,
					string profile_out_fname);
double ft_steady_erate_flat(double SS_erate,
					double start_zeta, double start_h,
					double end_time,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,
					string profile_out_fname);

//void run_ft_from_simple_erh_print_final_profile(
//				    string model_run_params_fname,
//					string sed_trans_param_fname,
//					string ft_parameter_fname,string ft_initial_profile_fname,
//					string erh_fname, string profile_out_fname);

void run_ft_from_erh_list_print_zh_final(int start_erh_runnumber,
					int end_erh_runnumber, string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,string ft_initial_profile_fname,
					string erh_datfile_name,
					string zeta_out_fname, string h_out_fname);
void run_ft_from_erh_print_zhe(double print_interval,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,string ft_initial_profile_fname,
					string erh_fname,
					string zeta_out_fname, string eta_out_fname, string h_out_fname);
flowtube flowtube_initializer(string sed_trans_param_fname,
					string ft_parameter_fname,
					string ft_initial_profile_fname);
flowtube flowtube_initializer(string sed_trans_param_fname,
					string ft_parameter_fname);
void create_meas_zeta_h_from_modeled(int run_number,
					string steady_modeled_fname,
					string zeta_final_fname,
					string h_final_fname,
					string false_zeta_meas_fname,
					string false_h_meas_fname);
void compare_model_with_measured(string measured_zeta_fname,
					string measured_h_fname,
					string final_zeta_fname,
					string final_h_fname,
					string steady_modeled_fname,
					string RMSE_fname);
void linear_interpolation_setup(vector<double> interp_x_loc,
					vector<double> modeled_x_loc,
					vector<int>& ds_interp_node_num,
					vector<int>& us_interp_node_num,
					vector<double>& interpolated_distance_fraction);
vector<double> linear_interpolate_vector(vector<int> ds_interp_node_num,
					vector<int> us_interp_node_num,
					vector<double> interpolated_fractional_distance,
					vector<double> modeled_data);
void virtual_erate_history_analysis(int start_datafile,
					int end_datafile, int n_runs_per_datafile,
					int virtual_hs_datfile,
					int virtual_hs_runnumber,
					double t_spacing,
					vector<double> erate_RMSE_efolding_values,
					vector<double> window_times,
					string steady_in_fname,
					string virtual_analysis_out_fname);
void virtual_erate_history_analysis(int start_datafile,
					int end_datafile, int n_runs_per_datafile,
					int virtual_hs_datfile,
					int virtual_hs_runnumber,
					double t_spacing,
					vector<double> erate_RMSE_efolding_values,
					vector<double> window_times,
					string steady_in_fname,
					string virtual_analysis_out_fname,
					string path);
void virtual_erate_history_analysis(int start_datafile,
					int end_datafile, int n_runs_per_datafile,
					int virtual_hs_datfile,
					int virtual_hs_runnumber,
					double t_spacing,
					vector<double> erate_RMSE_efolding_values,
					string steady_in_fname,
					string virtual_analysis_out_fname,
					string pathname);
#endif
