// FT_util.cpp
// utility fie for the flowtube model
#ifndef FT_UTIL_CPP
#define FT_UTIL_CPP


#include "FT_util.hpp"
#include "LSDStatsTools.hpp"
#include "flowtube.hpp"
#include <iostream>
#include <time.h>
#include <cstdlib>
#include <string>
using namespace std;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Gets the RMSE of a measured and modeled vector
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double RMSE(vector<double> interp_meas, vector<double> modeled)
 {
   int sz = interp_meas.size();
   double err =0;
   double n_counter = 0;
   for (int i = 0; i< sz; i++)
    {
     //cout << "LINE 1039 modeled["<<i<<"]: " << modeled[i]
     //	<< " interp h: " << interp_meas[i] << endl;
     if (interp_meas[i]>0)
      {
       n_counter ++;

       err += (modeled[i]-interp_meas[i])*(modeled[i]-interp_meas[i]);
      }

    }
   err = sqrt(err/n_counter);
   if (err<1e-10)
   {
	  err=0;
   }

   //cout << "LINE 1049 h_err: " << err << " n_counter: " << n_counter << endl;
   return err;
 }
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Gets the RMSE of a measured and modeled vector
// the rmse values are weighted.
// the sum of the weights must equal 1
// for each node unnormalized weighting is
// exp((ind_var-ind_var_start)/efold_weight_dist)
// you then add up the unnormalized weightings
// and then divide the weightings by this number
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double RMSE(vector<double> interp_meas, vector<double> modeled,
			vector<double> independant_var, double ind_var_start,
			double efold_weight_dist)
{
	int sz = interp_meas.size();
	double err =0;
	vector<double> err_local(sz);
	vector<double> weight(sz);
	double weight_total=0;
	double n_counter = 0;
	for (int i = 0; i< sz; i++)
	{
		n_counter++;
		//cout << "LINE 68 i: "<< i << " mod: " <<modeled[i] << " and meas: " << interp_meas[i] << endl;
		//err_local[i] = sqrt((modeled[i]-interp_meas[i])*(modeled[i]-interp_meas[i]));
		err_local[i] = (modeled[i]-interp_meas[i])*(modeled[i]-interp_meas[i]);
		weight[i] = exp((independant_var[i]-ind_var_start)/efold_weight_dist);
		weight_total+=weight[i];
		//cout << "err_loc: " << err_local[i] << endl;
		//cout << "weight: " << weight[i] << endl;
	}
	//cout << "weight total: " << weight_total <<endl;

	//cout << "err last: " << sqrt(err_local[sz-1]) << endl;
	double RMSE_test = 0;
	for (int i = 0; i< sz; i++)
	{
		RMSE_test+=err_local[i];
		weight[i]=weight[i]/weight_total;
		err+= weight[i]*err_local[i];
	}
	RMSE_test = sqrt(RMSE_test/n_counter);
	err = sqrt(err);

	//cout << "LINE 88 RMSE_test: " << RMSE_test << endl;
	//err = err/double(sz);
	if (err<1e-10)
	{
	err=0;
	}

	// bug check
	//cout << endl << endl << endl;
	//weight_total = 0;
	//for (int i = 0; i< sz; i++)
	//{
	//
	//	weight_total+=weight[i];
	//	cout << weight[i] << endl;
	//}
	//
	//cout << " line 90 FT_util, weight_total = " << weight_total << endl;
	//cout << endl << endl << endl;
	//cout << "LINE 1049 h_err: " << err << " n_counter: " << n_counter << endl;
	return err;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Gets the RMSE of a measured and modeled vector
// the rmse values are weighted.
// the sum of the weights must equal 1
// for each node unnormalized weighting is
// exp((ind_var-ind_var_start)/efold_weight_dist)
// you then add up the unnormalized weightings
// and then divide the weightings by this number
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double RMSE_local(vector<double> interp_meas, vector<double> modeled,
			vector<double> independant_var, double ind_var_before_end)
{
	int sz = interp_meas.size();
	double last_ind_var = independant_var[sz-1];
	double RMSE_window_start = last_ind_var-ind_var_before_end;
	double err =0;
	vector<double> err_local(sz);
	vector<double> weight(sz);

	double n_counter = 0;
	for (int i = 0; i< sz; i++)
	{

		if (independant_var[i]>=RMSE_window_start)
		{
			n_counter++;
			err+= (modeled[i]-interp_meas[i])*(modeled[i]-interp_meas[i]);
		}
	}
	err = sqrt(err/n_counter);

	//cout << "LINE 88 RMSE_test: " << RMSE_test << endl;
	//err = err/double(sz);
	if (err<1e-10)
	{
	err=0;
	}

	return err;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this runction calculates the root square mean error, but the data points
// that are used are normalized by teh mean of each dataset. This is
// used for topography and means a low RMSE hillslope needs to have similar
// relief to measured values  for a low RMSE value
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double RMSE_mean(vector<double> interp_meas, vector<double> modeled)
 {
   int sz = interp_meas.size();
   double sum_interp = 0;
   double sum_modeled = 0;
   double mean_interp;
   double mean_modeled;
   vector<double> normalized_interp;
   vector<double> normalized_modeled;

   for (int i = 0; i<sz; i++)
   {
	   sum_interp+=interp_meas[i];
	   sum_modeled+=modeled[i];
   }
   mean_interp = sum_interp/double(sz);
   mean_modeled = sum_modeled/double(sz);

   for (int i = 0; i<sz; i++)
   {
	   normalized_interp.push_back(interp_meas[i]-mean_interp);
	   normalized_modeled.push_back(modeled[i]-mean_modeled);
   }

   double err =0;
   double n_counter = 0;
   for (int i = 0; i< sz; i++)
    {
     //cout << "LINE 1039 modeled["<<i<<"]: " << normalized_modeled[i]
     //		<< " interp h: " << normalized_interp[i]<< endl;
     if (interp_meas[i]>0)
      {
       n_counter ++;

       err += (normalized_modeled[i]-normalized_interp[i])*
              (normalized_modeled[i]-normalized_interp[i]);
      }
    }
   err = sqrt(err/n_counter);
   if (err<1e-10)
   {
   	  err=0;
   }
   //cout << "LINE 1049 h_err: " << err << " n_counter: " << n_counter << endl;
   return err;
 }
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// RMSE from two files
// calucaltes the RMSE from two files of h and zeta
// used to test influence of initial conditions
// IT IS IMPORTANT TO NOTE THAT THE TWO SETS OF FILES SHOULD HAVE THE
// SAME NUMBER OF RUNS
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void RMSE_zeta_h_twofiles(int n_runs,string zeta_name1, string zeta_name2,
					string h_name1, string h_name2,
					string steady_modeled_fname,
					string RMSE_compare_out_fname)
{

	// open the data files
	ifstream zeta1_in, zeta2_in, h1_in, h2_in;
	zeta1_in.open(zeta_name1.c_str());
	zeta2_in.open(zeta_name2.c_str());
	h1_in.open(h_name1.c_str());
	h2_in.open(h_name2.c_str());

	ofstream RMSE_compare_out;
	RMSE_compare_out.open(RMSE_compare_out_fname.c_str());

	// get the number of nodes in teh flowtubes
	int n_nodes = 0;
	string temp;
	ifstream sm_in;
	sm_in.open(steady_modeled_fname.c_str());
	while (sm_in >> temp >> temp >> temp >> temp >> temp)
	{
		n_nodes++;
	}
	sm_in.close();

	// initialize the paired vectors
	vector<double> zeta1(n_nodes);
	vector<double> zeta2(n_nodes);
	vector<double> h1(n_nodes);
	vector<double> h2(n_nodes);

	double temp_double;
	double RMSE_zeta;
	double RMSE_h;
	// now loop through the runs, getting the RMSE between vectors
	for (int rn = 0; rn<n_runs; rn++)
	{
		zeta1_in >> temp >> temp;
		zeta2_in >> temp >> temp;
		h1_in >> temp >> temp;
		h2_in >> temp >> temp;

		for (int n_num = 0; n_num < n_nodes; n_num++)
		{
			zeta1_in >> temp_double;
			zeta1[n_num] = temp_double;
			zeta2_in >> temp_double;
			zeta2[n_num] = temp_double;
			h1_in >> temp_double;
			h1[n_num] = temp_double;
			h2_in >> temp_double;
			h2[n_num] = temp_double;
		}


		RMSE_zeta = RMSE_mean(zeta1, zeta2);
		RMSE_h = RMSE(h1, h2);

		RMSE_compare_out << RMSE_zeta << " " << RMSE_h << endl;
	}

	RMSE_compare_out.close();
	zeta1_in.close();
	zeta2_in.close();
	h1_in.close();
	h2_in.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// fraction in window: this function returns the frasction of time
// the simulated erosion rate is within a window of the
// 'true' erosion rate
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double fraction_in_window(double window_time,
				vector<double> time_vec,
				vector<double> erate_vec,
				vector<double> upper_window,
				vector<double> lower_window)
{
	int n_ts = time_vec.size();
	double final_time = time_vec[n_ts-1];

	int tot_counter =0;
	int inwindow_counter = 0;

	// loop through all time (this could be more efficient)
	for(int i = 0; i<n_ts; i++)
	{
		// if the time is within the appropriate time frame
		if (time_vec[i] >= final_time-window_time)
		{
			tot_counter++;				// increment the counter
			// and if the erate is within the bounds of the
			// 'true' erate
			if (erate_vec[i] >= lower_window[i] &&
			    erate_vec[i] <= upper_window[i])
			{
				inwindow_counter++;	// increment a counter for
									// being within the bounds
			}
		}
	}
	//  calcualte the fraction of time in the window
	double fraction = double(double(inwindow_counter)/double(tot_counter));

	return fraction;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function loads a flux history file
//The files have two columns, the first is the time and the
//second is the erosion rate. Each period of constant erosion is
//recorded with two points for the beginning and the end of that
//erosion rate, in order to plot a step function in matlab. So it looks like
//
//start_t1 r1
//end_t1 r1
//start_t2 r2
//end_t2 r2
//
// and so on. The start time of tn is the same as the end time of t_n-1.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double load_erate_history(ifstream& er_in,vector<double>& erate_change_times,
						vector<double>& erate_change_rates)
{
	vector<double> erct;
	vector<double> ercr;

	double temp_time;
	double temp_rate;
	double end_time;
	while (er_in >> temp_time >> temp_rate)
	{
		erct.push_back(temp_time);
		ercr.push_back(temp_rate);
	}
	erate_change_times = erct;
	erate_change_rates = ercr;
	end_time = temp_time;
	return end_time;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function updates the erosion rate using the flux history
// teh flux history file is loaded into two vectors.
// the strategy is to only call this function when needed
// it is called when the time hits some threshold.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double erosion_rate_from_erate_history_file(int& erate_change_index,
						double& erate_change_time,
						vector<double> erate_change_times,
						vector<double> erate_change_rates)
{
	// the changing rates come in pairs. so if the change time is triggered
	// this function is called. The new rate has the erate_change_index+1
	// and the new time has the erate_change_index+2
	double newrate = erate_change_rates[erate_change_index+1];
	erate_change_time = erate_change_times[erate_change_index+1];
	erate_change_index = erate_change_index+2;

	return newrate;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function advances the flux history to an arbitrary start time
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double erate_history_arbitrary_start(int& erate_change_index,
						double& erate_change_time,
						vector<double> erate_change_times,
						vector<double> erate_change_rates,
						double start_time)
{
	int eri = 0;
	double newrate =erate_change_rates[eri];
	int n_changes = erate_change_times.size();

	while (eri<n_changes)
	{
		//cout << "eri: " << eri << " and n_changes: " << n_changes << endl;
		if (erate_change_times[eri] <= start_time &&
		    erate_change_times[eri+1] > start_time)
		{
			erate_change_index = eri;
			erate_change_time = erate_change_times[eri+1];
			newrate = erate_change_rates[eri];
			eri = n_changes;
		}
		else
		{
			eri++;
		}
	}
	return newrate;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// generate_regularly_spaced_erate
// this function takes the name of an erate history and
// places it in a vector of regular time intervals.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void generate_regularly_spaced_erate_hist(
						double t_spacing,
						int run_number,
						ifstream& erf_datfile,
						vector<double>& er_times,
						vector<double>& erate_hist_regular)
{
	double t_ime = 0;
	double start_time;

	// vectors for holding individual erate histories
	vector<double> erate_change_times;
	vector<double> erate_change_rates;

	vector<double> reg_times;
	vector<double> reg_rates;

	generate_erate_hist_from_datfile(erf_datfile,
					run_number, erate_change_times,
					erate_change_rates);


	// initialize some variables
	double erate_change_time;	// the time that erate changes
	int erate_change_index;	// an index into the erate change history vector
	double erate;			// the erosion rate
	int n_changes;			// the number of nodes in the erate_change_times vector
	double end_time;		// the time at which the model run ends


	// reset the erate change time and erate change index
	erate_change_time = 0;
	erate_change_index = 0;

	// reset the time
	t_ime = 0;
	start_time = 0;

	// get the number of time changes
	n_changes = erate_change_times.size();

	// get the end_time
	end_time = erate_change_times[n_changes-1];

	// get the first change time
	erate =  erate_history_arbitrary_start(erate_change_index,
							erate_change_time,erate_change_times,
							erate_change_rates,start_time);

	//cout << "LINE 432 erate change time: " << erate_change_time << endl;

	// record the first time and rate
	reg_times.push_back(start_time);
	reg_rates.push_back(erate);

	// now loop through the erosion history
	double small_spacing = 1;
	int spacing_trigger = int(double(double(t_spacing)/double(small_spacing)));
	int st = 0;
	while (t_ime < end_time)
	{
		//t_ime += t_spacing;
		t_ime += small_spacing;
		st++;

		if (t_ime >= erate_change_time && t_ime < end_time)
		{
			erate = erosion_rate_from_erate_history_file(erate_change_index,
							erate_change_time,
							erate_change_times,
							erate_change_rates);
			//cout << "LINE 450 new_rate: " << erate << " and change_time: " << erate_change_time << endl;
		}

		if (st%spacing_trigger == 0)
		{
			// record the first time and rate
			reg_times.push_back(t_ime);
			reg_rates.push_back(erate);
		}

	}

	// change vectors for return
	er_times = reg_times;
	erate_hist_regular= reg_rates;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// compare_single_erate_hist_with_others
// this function compares a single erate history with a number of erate
// histories
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void compare_single_erate_hist_with_others(
					double t_spacing,
					int run_number_comparison_erh,
					int start_compare_rn,
					int end_compare_rn,
					double efold_weight_dist,
					string erh_fname,
					string RMSE_erate_fname)
{

	// open the erate history file
	ifstream erf_datfile;
	erf_datfile.open(erh_fname.c_str());

	vector<double> compare_reg_times;
	vector<double> compare_reg_rates;

	// now get the erate history that will be compared
	// with all the others
	generate_regularly_spaced_erate_hist(
						t_spacing,
						run_number_comparison_erh,
						erf_datfile,
						compare_reg_times,
						compare_reg_rates);
	erf_datfile.close();

	int n_ts = compare_reg_times.size();
	double ind_var_start =compare_reg_times[n_ts-1];

	// now reopen the erf_datfile and loop through
	// the erates
	ifstream erh_datfile;
	erh_datfile.open(erh_fname.c_str());

	vector<double> reg_times;
	vector<double> reg_rates;


	vector<double> RMSE_erates;

	// loop through the erate histories
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		generate_regularly_spaced_erate_hist(
						t_spacing,i,erh_datfile,
						reg_times, reg_rates);

		if (efold_weight_dist==0)
		{
			RMSE_erates.push_back(RMSE(compare_reg_rates, reg_rates));
		}
		else
		{
			RMSE_erates.push_back(RMSE(compare_reg_rates, reg_rates,
					reg_times, ind_var_start,
					efold_weight_dist));
		}
		//
		//cout << endl << endl << reg_rates[30] << "and RMSE: "
		//     << RMSE(compare_reg_rates, reg_rates) << endl << endl;
	}
	erh_datfile.close();

	// now print the RMSE to file
	ofstream RMSE_erate_out;
	RMSE_erate_out.open(RMSE_erate_fname.c_str());
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		RMSE_erate_out << RMSE_erates[i-start_compare_rn] << endl;
	}
	RMSE_erate_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// compare_single_erate_hist_with_others
// this function compares a single erate history with a number of erate
// histories
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void compare_single_erate_hist_with_others(
					double t_spacing,
					int run_number_comparison_erh,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> efold_weight_dist_vec,
					string erh_fname,
					string RMSE_erate_fname)
{
	// open the erate history file
	ifstream erf_datfile;
	erf_datfile.open(erh_fname.c_str());

	vector<double> compare_reg_times;
	vector<double> compare_reg_rates;

	// now get the erate history that will be compared
	// with all the others
	generate_regularly_spaced_erate_hist(
						t_spacing,
						run_number_comparison_erh,
						erf_datfile,
						compare_reg_times,
						compare_reg_rates);
	erf_datfile.close();

	int n_ts = compare_reg_times.size();
	double ind_var_start =compare_reg_times[n_ts-1];

	// now reopen the erf_datfile and loop through
	// the erates
	ifstream erh_datfile;
	erh_datfile.open(erh_fname.c_str());

	vector<double> reg_times;
	vector<double> reg_rates;

	int n_efold = efold_weight_dist_vec.size();
	double efold_weight_dist;

	vector< vector<double> > RMSE_erates_multiple_ef;
	vector< vector<double> >::iterator v_iter;
	vector<double> temp_vec;
	for (int i = 0; i<n_efold; i++)
	{
		RMSE_erates_multiple_ef.push_back(temp_vec);
	}

	// loop through the erate histories
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		generate_regularly_spaced_erate_hist(
						t_spacing,i,erh_datfile,
						reg_times, reg_rates);

		v_iter = RMSE_erates_multiple_ef.begin();
		// loop through the e_folding times
		for (int efn = 0; efn<n_efold; efn++)
		{
			efold_weight_dist = efold_weight_dist_vec[efn];
			if (efold_weight_dist==0)
			{
				(*v_iter).push_back(RMSE(compare_reg_rates, reg_rates));
			}
			else
			{
				(*v_iter).push_back(RMSE(compare_reg_rates, reg_rates,
						reg_times, ind_var_start,
						efold_weight_dist));
			}
			v_iter++;
		}

	}
	erh_datfile.close();

	// now print the RMSE to file
	ofstream RMSE_erate_out;
	RMSE_erate_out.open(RMSE_erate_fname.c_str());
	RMSE_erate_out << "-99 ";
	for (int efn = 0; efn<n_efold; efn++)
	{
		RMSE_erate_out <<  efold_weight_dist_vec[efn] << " ";
	}
	RMSE_erate_out << endl;
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		v_iter = RMSE_erates_multiple_ef.begin();
		RMSE_erate_out << i << " ";
		for (int efn = 0; efn<n_efold; efn++)
		{
			RMSE_erate_out << (*v_iter)[i-start_compare_rn] << " ";
			v_iter++;
		}
		RMSE_erate_out << endl;
	}
	RMSE_erate_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// compare_single_erate_hist_with_others
// this function compares a single erate history with a number of erate
// histories
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void compare_single_erate_hist_with_others(
					double t_spacing,
					vector<double> compare_reg_times,
					vector<double> compare_reg_rates,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> efold_weight_dist_vec,
					string erh_fname,
					string RMSE_erate_fname)
{
	int n_ts = compare_reg_times.size();
	double ind_var_start =compare_reg_times[n_ts-1];

	// now reopen the erf_datfile and loop through
	// the erates
	ifstream erh_datfile;
	erh_datfile.open(erh_fname.c_str());

	vector<double> reg_times;
	vector<double> reg_rates;

	int n_efold = efold_weight_dist_vec.size();
	double efold_weight_dist;

	vector< vector<double> > RMSE_erates_multiple_ef;
	vector< vector<double> >::iterator v_iter;
	vector<double> temp_vec;
	for (int i = 0; i<n_efold; i++)
	{
		RMSE_erates_multiple_ef.push_back(temp_vec);
	}

	// loop through the erate histories
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		generate_regularly_spaced_erate_hist(
						t_spacing,i,erh_datfile,
						reg_times, reg_rates);

		int n_times = reg_times.size();
		for (int nn = 0; nn<n_times; nn++)
		{
			//cout << "LINE 697 i: " << nn << " " << reg_times[nn] << " "
			//     << reg_rates[nn] << endl;
		}

		v_iter = RMSE_erates_multiple_ef.begin();
		// loop through the e_folding times
		for (int efn = 0; efn<n_efold; efn++)
		{
			efold_weight_dist = efold_weight_dist_vec[efn];
			if (efold_weight_dist==0)
			{
				(*v_iter).push_back(RMSE(compare_reg_rates, reg_rates));
			}
			else
			{
				(*v_iter).push_back(RMSE(compare_reg_rates, reg_rates,
						reg_times, ind_var_start,
						efold_weight_dist));
			}
			v_iter++;
		}

	}
	erh_datfile.close();

	// now print the RMSE to file
	ofstream RMSE_erate_out;
	RMSE_erate_out.open(RMSE_erate_fname.c_str());
	RMSE_erate_out << "-99 ";
	for (int efn = 0; efn<n_efold; efn++)
	{
		RMSE_erate_out <<  efold_weight_dist_vec[efn] << " ";
	}
	RMSE_erate_out << endl;
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		v_iter = RMSE_erates_multiple_ef.begin();
		RMSE_erate_out << i << " ";
		for (int efn = 0; efn<n_efold; efn++)
		{
			RMSE_erate_out << (*v_iter)[i-start_compare_rn] << " ";
			v_iter++;
		}
		RMSE_erate_out << endl;
	}
	RMSE_erate_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// compare_single_erate_hist_with_others
// this function compares a single erate history with a number of erate
// histories
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void compare_single_erate_hist_block_RMSE(
					double t_spacing,
					vector<double> compare_reg_times,
					vector<double> compare_reg_rates,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> last_block_vec,
					string erh_fname,
					string RMSE_erate_fname)
{
	int n_ts = compare_reg_times.size();
	//double ind_var_start =compare_reg_times[n_ts-1];

	// now reopen the erf_datfile and loop through
	// the erates
	ifstream erh_datfile;
	erh_datfile.open(erh_fname.c_str());

	vector<double> reg_times;
	vector<double> reg_rates;

	int n_efold = last_block_vec.size();
	double last_block;

	vector< vector<double> > RMSE_erates_multiple_ef;
	vector< vector<double> >::iterator v_iter;
	vector<double> temp_vec;
	for (int i = 0; i<n_efold; i++)
	{
		RMSE_erates_multiple_ef.push_back(temp_vec);
	}

	// loop through the erate histories
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		generate_regularly_spaced_erate_hist(
						t_spacing,i,erh_datfile,
						reg_times, reg_rates);

		int n_times = reg_times.size();
		for (int nn = 0; nn<n_times; nn++)
		{
			//cout << "LINE 697 i: " << nn << " " << reg_times[nn] << " "
			//     << reg_rates[nn] << endl;
		}

		v_iter = RMSE_erates_multiple_ef.begin();
		// loop through the e_folding times
		for (int efn = 0; efn<n_efold; efn++)
		{
			last_block = last_block_vec[efn];
			if (last_block==0)
			{
				(*v_iter).push_back(RMSE(compare_reg_rates, reg_rates));
			}
			else
			{
				(*v_iter).push_back(RMSE_local(compare_reg_rates, reg_rates,
						reg_times, last_block));
			}
			v_iter++;
		}

	}
	erh_datfile.close();

	// now print the RMSE to file
	ofstream RMSE_erate_out;
	RMSE_erate_out.open(RMSE_erate_fname.c_str());
	RMSE_erate_out << "-99 ";
	for (int efn = 0; efn<n_efold; efn++)
	{
		RMSE_erate_out <<  last_block_vec[efn] << " ";
	}
	RMSE_erate_out << endl;
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		v_iter = RMSE_erates_multiple_ef.begin();
		RMSE_erate_out << i << " ";
		for (int efn = 0; efn<n_efold; efn++)
		{
			RMSE_erate_out << (*v_iter)[i-start_compare_rn] << " ";
			v_iter++;
		}
		RMSE_erate_out << endl;
	}
	RMSE_erate_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// compare_erate_window
// this function generates several windows of flux about the
// 'true' erostion rate history. It then tests if simulated
// erosion rates lie within this window.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void compare_erate_windows(double t_spacing,
					int run_number_comparison_erh,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> window_times_vec,
					string window_2p_out_fname,
					string window_5p_out_fname,
					string window_10p_out_fname,
					string erh_fname)
{
	// open the erate history file
	ifstream erf_datfile;
	erf_datfile.open(erh_fname.c_str());

	vector<double> compare_reg_times;
	vector<double> compare_reg_rates;

	// now get the erate history that will be compared
	// with all the others
	generate_regularly_spaced_erate_hist(
						t_spacing,
						run_number_comparison_erh,
						erf_datfile,
						compare_reg_times,
						compare_reg_rates);
	erf_datfile.close();

	int n_ts = compare_reg_times.size();
	//double ind_var_start =compare_reg_times[n_ts-1];


	// now generate the windows
	vector<double> upper_2percent(n_ts);
	vector<double> lower_2percent(n_ts);
	vector<double> upper_5percent(n_ts);
	vector<double> lower_5percent(n_ts);
	vector<double> upper_10percent(n_ts);
	vector<double> lower_10percent(n_ts);

	// now loop through the 'true' erate history compiling the
	// erate windows
	for (int i = 0; i< n_ts; i++)
	{
		upper_2percent[i] = compare_reg_rates[i]+compare_reg_rates[i]*0.02;
		lower_2percent[i] = compare_reg_rates[i]-compare_reg_rates[i]*0.02;
		upper_5percent[i] = compare_reg_rates[i]+compare_reg_rates[i]*0.05;
		lower_5percent[i] = compare_reg_rates[i]-compare_reg_rates[i]*0.05;
		upper_10percent[i] = compare_reg_rates[i]+compare_reg_rates[i]*0.1;
		lower_10percent[i] = compare_reg_rates[i]-compare_reg_rates[i]*0.1;
	}

	// now reopen the erf_datfile and loop through
	// the erates
	ifstream erh_datfile;
	erh_datfile.open(erh_fname.c_str());

	vector<double> reg_times;
	vector<double> reg_rates;

	int n_windows = window_times_vec.size();
	double window_time;

	vector< vector<double> > windows_multiple_2percent;
	vector< vector<double> > windows_multiple_5percent;
	vector< vector<double> > windows_multiple_10percent;
	vector< vector<double> >::iterator v_iter_2p;
	vector< vector<double> >::iterator v_iter_5p;
	vector< vector<double> >::iterator v_iter_10p;
	vector<double> temp_vec;
	for (int i = 0; i<n_windows; i++)
	{
		windows_multiple_2percent.push_back(temp_vec);
		windows_multiple_5percent.push_back(temp_vec);
		windows_multiple_10percent.push_back(temp_vec);
	}

	// loop through the erate histories
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		generate_regularly_spaced_erate_hist(
						t_spacing,i,erh_datfile,
						reg_times, reg_rates);

		v_iter_2p = windows_multiple_2percent.begin();
		v_iter_5p = windows_multiple_5percent.begin();
		v_iter_10p = windows_multiple_10percent.begin();
		// loop through the e_folding times
		for (int wn = 0; wn<n_windows; wn++)
		{
			window_time = window_times_vec[wn];
			(*v_iter_2p).push_back(fraction_in_window(window_time,
								reg_times,reg_rates,
								upper_2percent,lower_2percent));
			(*v_iter_5p).push_back(fraction_in_window(window_time,
								reg_times,reg_rates,
								upper_5percent,lower_5percent));
			(*v_iter_10p).push_back(fraction_in_window(window_time,
								reg_times,reg_rates,
								upper_10percent,lower_10percent));
			v_iter_2p++;
			v_iter_5p++;
			v_iter_10p++;
		}

	}
	erh_datfile.close();

	// now print window fractions to file
	ofstream window_2p_out,window_5p_out,window_10p_out;
	window_2p_out.open(window_2p_out_fname.c_str());
	window_5p_out.open(window_5p_out_fname.c_str());
	window_10p_out.open(window_10p_out_fname.c_str());

	// print the window_times
	window_2p_out << "-99 ";
	window_5p_out << "-99 ";
	window_10p_out <<"-99 ";
	for (int wn = 0; wn<n_windows; wn++)
	{
		window_2p_out << window_times_vec[wn] << " ";
		window_5p_out << window_times_vec[wn] << " ";
		window_10p_out << window_times_vec[wn] << " ";
	}
	window_2p_out << endl;
	window_5p_out << endl;
	window_10p_out << endl;

	// now print the fraction that lies within the bounds
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		v_iter_2p = windows_multiple_2percent.begin();
		v_iter_5p = windows_multiple_5percent.begin();
		v_iter_10p = windows_multiple_10percent.begin();
		window_2p_out << i<<" ";
		window_5p_out << i<<" ";
		window_10p_out<<i<<" ";
		for (int wn = 0; wn<n_windows; wn++)
		{
			window_2p_out << (*v_iter_2p)[i-start_compare_rn] << " ";
			window_5p_out << (*v_iter_5p)[i-start_compare_rn] << " ";
			window_10p_out << (*v_iter_10p)[i-start_compare_rn] << " ";
			v_iter_2p++;
			v_iter_5p++;
			v_iter_10p++;
		}
		window_2p_out << endl;
		window_5p_out << endl;
		window_10p_out << endl;
	}
	window_2p_out.close();
	window_5p_out.close();
	window_10p_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// compare_erate_window
// this function generates several windows of flux about the
// 'true' erostion rate history. It then tests if simulated
// erosion rates lie within this window.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void compare_erate_windows(double t_spacing,
					vector<double> compare_reg_times,
					vector<double> compare_reg_rates,
					int start_compare_rn,
					int end_compare_rn,
					vector<double> window_times_vec,
					string window_2p_out_fname,
					string window_5p_out_fname,
					string window_10p_out_fname,
					string erh_fname)
{
	int n_ts = compare_reg_times.size();
	//double ind_var_start =compare_reg_times[n_ts-1];


	// now generate the windows
	vector<double> upper_2percent(n_ts);
	vector<double> lower_2percent(n_ts);
	vector<double> upper_5percent(n_ts);
	vector<double> lower_5percent(n_ts);
	vector<double> upper_10percent(n_ts);
	vector<double> lower_10percent(n_ts);

	// now loop through the 'true' erate history compiling the
	// erate windows
	for (int i = 0; i< n_ts; i++)
	{
		upper_2percent[i] = compare_reg_rates[i]+compare_reg_rates[i]*0.02;
		lower_2percent[i] = compare_reg_rates[i]-compare_reg_rates[i]*0.02;
		upper_5percent[i] = compare_reg_rates[i]+compare_reg_rates[i]*0.05;
		lower_5percent[i] = compare_reg_rates[i]-compare_reg_rates[i]*0.05;
		upper_10percent[i] = compare_reg_rates[i]+compare_reg_rates[i]*0.1;
		lower_10percent[i] = compare_reg_rates[i]-compare_reg_rates[i]*0.1;
	}

	// now reopen the erf_datfile and loop through
	// the erates
	ifstream erh_datfile;
	erh_datfile.open(erh_fname.c_str());

	vector<double> reg_times;
	vector<double> reg_rates;

	int n_windows = window_times_vec.size();
	double window_time;

	vector< vector<double> > windows_multiple_2percent;
	vector< vector<double> > windows_multiple_5percent;
	vector< vector<double> > windows_multiple_10percent;
	vector< vector<double> >::iterator v_iter_2p;
	vector< vector<double> >::iterator v_iter_5p;
	vector< vector<double> >::iterator v_iter_10p;
	vector<double> temp_vec;
	for (int i = 0; i<n_windows; i++)
	{
		windows_multiple_2percent.push_back(temp_vec);
		windows_multiple_5percent.push_back(temp_vec);
		windows_multiple_10percent.push_back(temp_vec);
	}

	// loop through the erate histories
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		generate_regularly_spaced_erate_hist(
						t_spacing,i,erh_datfile,
						reg_times, reg_rates);

		v_iter_2p = windows_multiple_2percent.begin();
		v_iter_5p = windows_multiple_5percent.begin();
		v_iter_10p = windows_multiple_10percent.begin();
		// loop through the e_folding times
		for (int wn = 0; wn<n_windows; wn++)
		{
			window_time = window_times_vec[wn];
			(*v_iter_2p).push_back(fraction_in_window(window_time,
								reg_times,reg_rates,
								upper_2percent,lower_2percent));
			(*v_iter_5p).push_back(fraction_in_window(window_time,
								reg_times,reg_rates,
								upper_5percent,lower_5percent));
			(*v_iter_10p).push_back(fraction_in_window(window_time,
								reg_times,reg_rates,
								upper_10percent,lower_10percent));
			v_iter_2p++;
			v_iter_5p++;
			v_iter_10p++;
		}

	}
	erh_datfile.close();

	// now print window fractions to file
	ofstream window_2p_out,window_5p_out,window_10p_out;
	window_2p_out.open(window_2p_out_fname.c_str());
	window_5p_out.open(window_5p_out_fname.c_str());
	window_10p_out.open(window_10p_out_fname.c_str());

	// print the window_times
	window_2p_out << "-99 ";
	window_5p_out << "-99 ";
	window_10p_out <<"-99 ";
	for (int wn = 0; wn<n_windows; wn++)
	{
		window_2p_out << window_times_vec[wn] << " ";
		window_5p_out << window_times_vec[wn] << " ";
		window_10p_out << window_times_vec[wn] << " ";
	}
	window_2p_out << endl;
	window_5p_out << endl;
	window_10p_out << endl;

	// now print the fraction that lies within the bounds
	for (int i = start_compare_rn; i<=end_compare_rn; i++)
	{
		v_iter_2p = windows_multiple_2percent.begin();
		v_iter_5p = windows_multiple_5percent.begin();
		v_iter_10p = windows_multiple_10percent.begin();
		window_2p_out << i<<" ";
		window_5p_out << i<<" ";
		window_10p_out<<i<<" ";
		for (int wn = 0; wn<n_windows; wn++)
		{
			window_2p_out << (*v_iter_2p)[i-start_compare_rn] << " ";
			window_5p_out << (*v_iter_5p)[i-start_compare_rn] << " ";
			window_10p_out << (*v_iter_10p)[i-start_compare_rn] << " ";
			v_iter_2p++;
			v_iter_5p++;
			v_iter_10p++;
		}
		window_2p_out << endl;
		window_5p_out << endl;
		window_10p_out << endl;
	}
	window_2p_out.close();
	window_5p_out.close();
	window_10p_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// create erate history
// this function has two components.
// the first is the change in erosion rate
// it can be either positive or negative
// there are several switched
// first one can have a negative exponential
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void create_erate_history(long seed,ofstream& erh_out,
						double start_time, double start_rate,
						double end_time, double dt,
						double max_eros, double min_eros,
						double mean_change, double mean_interval,
                        vector<double>& erate_change_times,
						vector<double>& erate_change_rates)
{
	vector<double> ert;
	vector<double> err;

	double t_ime = start_time;
	double old_time = start_time;
	double old_rate = start_rate;
	double new_rate;
	double new_time;
	double ran_time;
	int n_ts_ran_time;
	double pos_or_neg;
	double erate_change;
	double beyond_threshold_erate;
	// loop through the history. with each step you add
	// additional erosion rates.
	while (t_ime < end_time)
	{
		// get the new time
		// in this function the time is exponentially distributed
		ran_time = expdev(&seed)*mean_interval;
		n_ts_ran_time = int (double(ran_time/dt+0.5));
		//cout << "raw ran_time: " << ran_time << " and increments: "
		//		<< n_ts_ran_time << endl;
		ran_time = double(n_ts_ran_time)*dt;
		//cout << "fixed ran_time: " << ran_time << endl;
		new_time += ran_time;

		// now get the change in the erosion rate
		// first determine if the change is positive or negative
		pos_or_neg = ran3(&seed);
		if (pos_or_neg >= 0.5)
		{
			pos_or_neg = 1.0;
		}
		else
		{
			pos_or_neg = -1.0;
		}

		// now get the cahnge. This is also exponentially distributed
		erate_change = mean_change*expdev(&seed)*pos_or_neg;

		// now change the erosion rate. This gets reflected if it lies
		// outside the minimum or maximum bounds
		new_rate = old_rate+erate_change;

		if (new_rate > max_eros)
		{
			beyond_threshold_erate = new_rate-max_eros;
			new_rate = max_eros-beyond_threshold_erate;
		}
		if (new_rate < min_eros)
		{
			beyond_threshold_erate = min_eros-new_rate;
			new_rate = min_eros+beyond_threshold_erate;
		}

		if (new_time > end_time)
		{
			new_time = end_time;
		}

		erh_out << old_time << " " << old_rate << endl;
		erh_out << new_time << " " << old_rate << endl;
		ert.push_back(old_time);
		ert.push_back(new_time);
		err.push_back(old_rate);
		err.push_back(old_rate);

		t_ime = new_time;
		old_time = new_time;
		old_rate = new_rate;
	}
	erate_change_times = ert;
	erate_change_rates = err;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// create erate history
// this function has two components.
// the first is the change in erosion rate
// it can be either positive or negative
// there are several switched
// first one can have a negative exponential
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void create_erate_history(long seed,
						double start_time, double start_rate,
						double end_time, double dt,
						double max_eros, double min_eros,
						double mean_change, double mean_interval,
                        vector<double>& erate_change_times,
						vector<double>& erate_change_rates)
{
	vector<double> ert;
	vector<double> err;

	double t_ime = start_time;
	double old_time = start_time;
	double old_rate = start_rate;
	double new_rate;
	double new_time;
	double ran_time;
	int n_ts_ran_time;
	double pos_or_neg;
	double erate_change;
	double beyond_threshold_erate;
	// loop through the history. with each step you add
	// additional erosion rates.
	while (t_ime < end_time)
	{
		// get the new time
		// in this function the time is exponentially distributed
		ran_time = expdev(&seed)*mean_interval;
		n_ts_ran_time = int (double(ran_time/dt+0.5));
		ran_time = double(n_ts_ran_time)*dt;
		new_time += ran_time;

		// now get the change in the erosion rate
		// first determine if the change is positive or negative
		pos_or_neg = ran3(&seed);
		if (pos_or_neg >= 0.5)
		{
			pos_or_neg = 1.0;
		}
		else
		{
			pos_or_neg = -1.0;
		}

		// now get the cahnge. This is also exponentially distributed
		erate_change = mean_change*expdev(&seed)*pos_or_neg;

		// now change the erosion rate. This gets reflected if it lies
		// outside the minimum or maximum bounds
		new_rate = old_rate+erate_change;

		if (new_rate > max_eros)
		{
			beyond_threshold_erate = new_rate-max_eros;
			new_rate = max_eros-beyond_threshold_erate;
		}
		if (new_rate < min_eros)
		{
			beyond_threshold_erate = min_eros-new_rate;
			new_rate = min_eros+beyond_threshold_erate;
		}

		if (new_time > end_time)
		{
			new_time = end_time;
		}

		ert.push_back(old_time);
		ert.push_back(new_time);
		err.push_back(old_rate);
		err.push_back(old_rate);

		t_ime = new_time;
		old_time = new_time;
		old_rate = new_rate;
	}
	erate_change_times = ert;
	erate_change_rates = err;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function creates an erate histfile of arbitrary length
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void create_erate_hist_datfile(int start_runnumber, int end_runnumber,
					double start_time,
					double end_time, double max_eros,
					double min_eros, double mean_change,
					double mean_interval, double start_rate,
					string model_run_params_fname,
					string erf_datfile_fname)
{

	int flux_switch;
	int prod_switch;
	double flux_us;
	double dt;

	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
							>> temp >> flux_us >> temp >> dt;
	model_run_params_in.close();


	long seed = time(NULL);

	// open the erf_datfile. This contains multiple erate histories
	ofstream erf_datfile_out;
	erf_datfile_out.open(erf_datfile_fname.c_str());

	// these two vectors hold the erate change information
	vector<double> erate_change_times;
	vector<double> erate_change_rates;
	vector<double> empty_vec;

	// these are placeholders
	double RMSE_zeta = -99;
	double RMSE_h = -99;

	// loop through the number of histories
	for (int i = start_runnumber; i<= end_runnumber; i++)
	{
		// create an erate history, this function copies
		// the hisotry to the two erate_change vectors
		create_erate_history(seed,
						start_time, start_rate,
						end_time, dt,
						max_eros, min_eros,
						mean_change, mean_interval,
                        erate_change_times,
						erate_change_rates);

		// now print the erate history to file
		print_erate_hist_datfile(erf_datfile_out, i,
								RMSE_zeta,RMSE_h,
								erate_change_times,
								erate_change_rates);
	}
	erf_datfile_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//	As you loop through simulations, you only want 1 output file:
// so each simulation has two rows of data. The first number in
// both rows is the simulation number. The second number in the
// top row is the RMSE of zeta and the second number in the bottom
// row is the RMSE of h. The third number in both rows is the number
// of points in the erate history, and the remaining numbers are the
// erate history, with the top row being the times and the
// bottom row being the rates.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void print_erate_hist_datfile(ofstream& erh_datfile,
							int run_number,
							double RMSE_zeta,double RMSE_h,
							vector<double> erate_change_times,
							vector<double> erate_change_rates)
{
	int sz = erate_change_times.size();
	erh_datfile << run_number << " " << RMSE_zeta << " " << sz;
	for (int i = 0; i<sz; i++)
	{
		erh_datfile << " " << erate_change_times[i];
	}
	erh_datfile << endl;

	erh_datfile << run_number << " " << RMSE_h << " " << sz;
	for (int i = 0; i<sz; i++)
	{
		erh_datfile << " " << erate_change_rates[i];
	}
	erh_datfile << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes and erate_hist_datfile and spits out an erate
// hist file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void generate_erate_hist_from_datfile(string erf_datfile_fname,
					string erf_out_fname, int run_number)
{
	int thisline_run_number;
	double temp;
	int continue_trigger = 0;
	int n_hist_steps;

	ifstream erf_datfile;
	erf_datfile.open(erf_datfile_fname.c_str());
	ofstream erf_out;
	erf_out.open(erf_out_fname.c_str());

	vector<double> ert;
	vector<double> err;

	while (continue_trigger == 0)
	{
		erf_datfile >> thisline_run_number;

		if (thisline_run_number == run_number)
		{
			erf_datfile >> temp >> n_hist_steps;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
				ert.push_back(temp);
			}
			erf_datfile >> temp >> temp >> temp;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
				err.push_back(temp);
			}
			continue_trigger = 1;
		}
		else
		{
			erf_datfile >> temp >> n_hist_steps;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
			}
			erf_datfile >> temp >> temp >> temp;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
			}
		}
	}

	for (int i = 0; i<n_hist_steps; i++)
	{
		erf_out << ert[i] << " " << err[i] << endl;
	}

	erf_datfile.close();
	erf_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes and erate_hist_datfile and spits out an erate
// hist pair of vectors as well as a erate hist file
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void generate_erate_hist_from_datfile(ifstream& erf_datfile,
					ofstream& erf_out, int run_number,
					vector<double>& erate_change_times,
					vector<double>& erate_change_rates)
{
	int thisline_run_number;
	double temp;
	int continue_trigger = 0;
	int n_hist_steps;

	vector<double> ert;
	vector<double> err;

	while (continue_trigger == 0)
	{
		erf_datfile >> thisline_run_number;

		if (thisline_run_number == run_number)
		{
			erf_datfile >> temp >> n_hist_steps;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
				ert.push_back(temp);
			}
			erf_datfile >> temp >> temp >> temp;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
				err.push_back(temp);
			}
			continue_trigger = 1;
		}
		else
		{
			erf_datfile >> temp >> n_hist_steps;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
			}
			erf_datfile >> temp >> temp >> temp;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
			}
		}
	}

	for (int i = 0; i<n_hist_steps; i++)
	{
		erf_out << ert[i] << " " << err[i] << endl;
	}
	erate_change_times = ert;
	erate_change_rates = err;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes and erate_hist_datfile and spits out an erate
// history vectors
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void generate_erate_hist_from_datfile(ifstream& erf_datfile,
					int run_number,
					vector<double>& erate_change_times,
					vector<double>& erate_change_rates)
{
	int thisline_run_number;
	double temp;
	int continue_trigger = 0;
	int n_hist_steps;

	vector<double> ert;
	vector<double> err;

	while (continue_trigger == 0)
	{
		erf_datfile >> thisline_run_number;

		if (thisline_run_number == run_number)
		{
			erf_datfile >> temp >> n_hist_steps;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
				ert.push_back(temp);
			}
			erf_datfile >> temp >> temp >> temp;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
				err.push_back(temp);
			}
			continue_trigger = 1;
		}
		else
		{
			erf_datfile >> temp >> n_hist_steps;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
			}
			erf_datfile >> temp >> temp >> temp;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
			}
		}
	}

	erate_change_times = ert;
	erate_change_rates = err;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes an erate_hist_datfile
// name and spits out an erate
// history vectors
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void generate_erate_hist_from_datfile(string erf_datfile_name,
					int run_number,
					vector<double>& erate_change_times,
					vector<double>& erate_change_rates)
{
	int thisline_run_number;
	double temp;
	int continue_trigger = 0;
	int n_hist_steps;

	vector<double> ert;
	vector<double> err;

	ifstream erf_datfile;
	erf_datfile.open(erf_datfile_name.c_str());

	while (continue_trigger == 0)
	{
		erf_datfile >> thisline_run_number;

		if (thisline_run_number == run_number)
		{
			erf_datfile >> temp >> n_hist_steps;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
				ert.push_back(temp);
			}
			erf_datfile >> temp >> temp >> temp;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
				err.push_back(temp);
			}
			continue_trigger = 1;
		}
		else
		{
			erf_datfile >> temp >> n_hist_steps;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
			}
			erf_datfile >> temp >> temp >> temp;
			for (int i = 0; i<n_hist_steps; i++)
			{
				erf_datfile >> temp;
			}
		}
	}

	erate_change_times = ert;
	erate_change_rates = err;
	erf_datfile.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=








/*

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes a single erate history and
// runs the flowtube model, printing the final profile
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void run_ft_from_simple_erh_print_final_profile(
				    string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,string ft_initial_profile_fname,
					string erh_fname, string profile_out_fname)
{
	// initialize some variables
	double erate;			// the erosion rate
	double old_erate;
	double old_erate_time;
	double new_erate;
	double new_erate_time;
	double erate_slope;
	double time_diff;
	int n_changes;			// the number of nodes in the
							// erate_change_times vector
	double end_time;		// the time at which the model run ends
	double ds_elev;			// elevation at teh downslope node

	double t_ime;
	double start_time;
	int tt;

	int flux_switch;
	int prod_switch;
	double flux_us;
	double dt;

	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt;
	model_run_params_in.close();

	// vectors for holding individual erate histories
	vector<double> erate_times;
	vector<double> erate_rates;

	// first open the erh_file
	ifstream erh_file_in;
	erh_file_in.open(erh_fname.c_str());

	// now load the erosion rates into a vector
	double temp_time, temp_erate;
	while (erh_file_in >> temp_time >> temp_erate)
	{
		erate_times.push_back(temp_time);
		erate_rates.push_back(temp_erate);
	}
	// close the file
	erh_file_in.close();

	//cout << " run number: " << endl;
	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname,
						  ft_initial_profile_fname);
	//cout << "flowtube initialized!" << endl;

	// reset the erate change time and erate change index
	double erate_change_time = 0;
	double erate_change_index = 0;


	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;

	// variables for the sediment transport loop
	vector<double> surface_change_rate(ft_test.get_n_nodes(),0.0);
										// the erosion rate of deposition rate
										// at the surface

	// get the number of time changes
	n_changes = erate_change_times.size();

	// get the starting downslope elevation
	ds_elev = ft_test.get_zeta_ds();

	// now get information about erate
	int n_ts = erate_times.size();
	end_time = estart_times[n_ts-1];

	int curr_erate_index = 0;
	old_erate = erate_rates[curr_erate_index];
	old_erate_time = erate_times[curr_erate_index];
	new_erate = erate_rates[curr_erate_index+1];
	new_erate_time = erate_times[curr_erate_index+1];
	erate_slope = (new_erate-old_erate)/(new_erate_time-old_erate_time);


	// now loop through the erosion history
	while (t_ime < end_time)
	{
		t_ime += dt;
		tt++;


		if (t_ime < new_erate_time)
		{
			time_diff = t_ime-old_erate_time;
			erate = time_diff*erate_slope+old_erate;
		}
		else
		{
			curr_erate_index++;
			old_erate = erate_rates[curr_erate_index];
			old_erate_time = erate_times[curr_erate_index];
			new_erate = erate_rates[curr_erate_index+1];
			new_erate_time = erate_times[curr_erate_index+1];
			erate_slope = (new_erate-old_erate)/(new_erate_time-old_erate_time);
			time_diff = t_ime-old_erate_time;
			erate = time_diff*erate_slope+old_erate;
		}

		ds_elev -= erate*dt;
		ft_test.flux_timestep_elev_bc(dt, flux_us, ds_elev,
						flux_switch, prod_switch,
						surface_change_rate);



	}

	// now print the result to a profile file
	ofstream profile_out;
	profile_out.precision(10);
	profile_out.open(profile_out_fname.c_str());
	ft_test.export_input_profile(profile_out);
	profile_out.close();

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
*/



















//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes a list of erate histories and
// runs the flowtube model, printing the final h and zeta
// data.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void run_ft_from_erh_list_print_zh_final(int start_erh_runnumber,
					int end_erh_runnumber,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,string ft_initial_profile_fname,
					string erh_datfile_name,
					string zeta_out_fname, string h_out_fname)
{
	// initialize some variables
	double erate_change_time;	// the time that erate changes
	int erate_change_index;	// an index into the erate change history vector
	double erate;			// the erosion rate
	int n_changes;			// the number of nodes in the erate_change_times vector
	double end_time;		// the time at which the model run ends
	double ds_elev;			// elevation at teh downslope node

	double t_ime;
	double start_time;

	int flux_switch;
	int prod_switch;
	double flux_us;
	double dt;

	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt;
	model_run_params_in.close();

	// vectors for holding individual erate histories
	vector<double> erate_change_times;
	vector<double> erate_change_rates;

	// first open the erh_file
	ifstream erh_datfile_in;
	erh_datfile_in.open(erh_datfile_name.c_str());

	// now open the zeta and h outfiles
	ofstream zeta_out,h_out;
	zeta_out.open(zeta_out_fname.c_str());
	h_out.open(h_out_fname.c_str());

	// now loop through the erate_histories
	// we initialize a flowtube during each loop
	for (int rn = start_erh_runnumber; rn<= end_erh_runnumber; rn++)
	{
		long run_start_time = time(NULL);
		//cout << " run number: " << endl;
		// initialize a flowtube
		flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
							  ft_parameter_fname,
							  ft_initial_profile_fname);
		//cout << "flowtube initialized!" << endl;

		// reset the erate change time and erate change index
		erate_change_time = 0;
		erate_change_index = 0;

		// reset the time
		t_ime = 0;
		start_time = 0;

		// variables for the sediemnt transport loop
		vector<double> surface_change_rate(ft_test.get_n_nodes(),0.0);
											// the erosion rate of deposition rate
											// at the surface

		// get the erate history of this run number
		generate_erate_hist_from_datfile(erh_datfile_in, rn,
					erate_change_times, erate_change_rates);

		// get the number of time changes
		n_changes = erate_change_times.size();

		// get the end_time
		end_time = erate_change_times[n_changes-1];

		// bug checking: print out the run number and the history
		//cout << endl << endl << endl;
		//cout << "run number: " << rn << endl;
		//cout << "end time: " << end_time << endl;
		//cout << "time_change    erosion rate" << endl;
		//for (int i = 0; i<n_changes; i++)
		//{
		//	cout << erate_change_times[i] << " " << erate_change_rates[i] << endl;
		//}

		// get the first change time
		erate =  erate_history_arbitrary_start(erate_change_index,
								erate_change_time,erate_change_times,
								erate_change_rates,start_time);

		// get the starting downslope elevation
		ds_elev = ft_test.get_zeta_ds();

		// now loop through the erosion history
		while (t_ime < end_time)
		{
			t_ime += dt;

			if (t_ime >= erate_change_time && t_ime < end_time)
			{
				erate = erosion_rate_from_erate_history_file(erate_change_index,
											erate_change_time,
											erate_change_times,
								erate_change_rates);
				//cout << "t_ime: " << t_ime << " and rate: " << erate
				//	 << " and t_change is: " << erate_change_time << endl;
			}
			ds_elev -= erate*dt;
			ft_test.flux_timestep_elev_bc(dt, flux_us, ds_elev,
							flux_switch, prod_switch,
							surface_change_rate);
		}

		// now print the final h and zeta (using the relative zeta function
		zeta_out << rn << " ";
		ft_test.print_relative_zeta(t_ime,zeta_out);
		h_out << rn << " ";
		ft_test.print_h(t_ime,h_out);

		long run_end_time = time(NULL);
		cout << "run number: " << rn
		     << " and runtime: " << run_end_time-run_start_time << endl;
	}

	// close the files
	erh_datfile_in.close();
	zeta_out.close();
	h_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes a single erate history and
// runs the flowtube model, printing the h zeta, and eta
// data at specified intervals
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void run_ft_from_erh_print_zhe(double print_interval,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,string ft_initial_profile_fname,
					string erh_fname,
					string zeta_out_fname, string eta_out_fname, string h_out_fname)
{
	// initialize some variables
	double erate_change_time;	// the time that erate changes
	int erate_change_index;	// an index into the erate change history vector
	double erate;			// the erosion rate
	int n_changes;			// the number of nodes in the
							// erate_change_times vector
	double end_time;		// the time at which the model run ends
	double ds_elev;			// elevation at teh downslope node

	double t_ime;
	double start_time;
	int tt;

	int flux_switch;
	int prod_switch;
	double flux_us;
	double dt;

	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt;
	model_run_params_in.close();

	const int p_i = int(double(print_interval/dt+0.5));

	// vectors for holding individual erate histories
	vector<double> erate_change_times;
	vector<double> erate_change_rates;

	// first open the erh_file
	ifstream erh_file_in;
	erh_file_in.open(erh_fname.c_str());

	// now open the zeta, eta, and h outfiles
	ofstream zeta_out,h_out,eta_out;
	zeta_out.open(zeta_out_fname.c_str());
	h_out.open(h_out_fname.c_str());
	eta_out.open(eta_out_fname.c_str());

	//cout << " run number: " << endl;
	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname,
						  ft_initial_profile_fname);
	//cout << "flowtube initialized!" << endl;

	// reset the erate change time and erate change index
	erate_change_time = 0;
	erate_change_index = 0;

	// load in the erate history
	ifstream erh_in;
	erh_in.open(erh_fname.c_str());
	end_time = load_erate_history(erh_in,erate_change_times,
						erate_change_rates);
	erh_in.close();

	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;

	// variables for the sediemnt transport loop
	vector<double> surface_change_rate(ft_test.get_n_nodes(),0.0);
										// the erosion rate of deposition rate
										// at the surface

	// get the number of time changes
	n_changes = erate_change_times.size();

	// get the first change time
	erate =  erate_history_arbitrary_start(erate_change_index,
							erate_change_time,erate_change_times,
							erate_change_rates,start_time);

	// get the starting downslope elevation
	ds_elev = ft_test.get_zeta_ds();

	// print the initial conditions
	ft_test.print_h_s(zeta_out);
	ft_test.print_h_s(eta_out);
	ft_test.print_h_s(h_out);

	ft_test.print_zeta(t_ime, zeta_out);
	ft_test.print_eta(t_ime, eta_out);
	ft_test.print_h(t_ime, h_out);

	// now loop through the erosion history
	while (t_ime < end_time)
	{
		t_ime += dt;
		tt++;

		if (t_ime >= erate_change_time && t_ime < end_time)
		{
			erate = erosion_rate_from_erate_history_file(erate_change_index,
										erate_change_time,
										erate_change_times,
							erate_change_rates);
			cout << "t_ime: " << t_ime << " and rate: " << erate
				 << " and t_change is: " << erate_change_time << endl;
		}
		ds_elev -= erate*dt;
		ft_test.flux_timestep_elev_bc(dt, flux_us, ds_elev,
						flux_switch, prod_switch,
						surface_change_rate);


		// print if the time is appropriate
		if (tt%p_i == 0)
		{
			cout << "FT_util, run_ft_from_erh_print_zhe: printing, time: "
				<< t_ime << endl;
			ft_test.print_zeta(t_ime, zeta_out);
			ft_test.print_eta(t_ime, eta_out);
			ft_test.print_h(t_ime, h_out);
		}
	}

	// close the files
	zeta_out.close();
	h_out.close();
	eta_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function runs a flowtube at a steady erosion rate
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void ft_steady_erate(double SS_erate, double end_time,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,
					string ft_initial_profile_fname,
					string profile_out_fname)
{
	// initialize some variables
	//double erate;			// the erosion rate
	double ds_elev;			// elevation at teh downslope node

	double t_ime;
	double start_time;
	int tt;

	int flux_switch;
	int prod_switch;
	double flux_us;
	double dt;

	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt;
	model_run_params_in.close();


	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname,
						  ft_initial_profile_fname);

	// raise the flowtube so the downstream boundary is at elevation
	// 100
	//ds_elev = 100.0;
	//ft_test.raise_zeta_eta_ds_bound(ds_elev);

	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;

	// variables for the sediment transport loop
	vector<double> surface_change_rate(ft_test.get_n_nodes(),0.0);
										// the erosion rate of deposition rate

	// get the starting downslope elevation
	ds_elev = ft_test.get_zeta_ds();
										// at the surface
	//ds_elev = 164;

	// now loop through time
	while(t_ime < end_time)
	{
		tt++;
		t_ime += dt;		// increment the time

		ds_elev -= dt*SS_erate;

		// run flux
		ft_test.flux_timestep_elev_bc(dt, flux_us, ds_elev,flux_switch, prod_switch,
							surface_change_rate);

		if (tt%100000 == 0)
		{
			cout << "time is: " << t_ime << endl;
		}
	}

	// now print the result to a profile file
	ofstream profile_out;
	profile_out.precision(10);
	profile_out.open(profile_out_fname.c_str());
	ft_test.export_input_profile(profile_out);
	profile_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function runs a flowtube at a steady erosion rate
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void ft_steady_erate(double SS_erate, double end_time,
					double ds_elev,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,
					string ft_initial_profile_fname,
					string profile_out_fname)
{
	// initialize some variables
	//double erate;			// the erosion rate

	double t_ime;
	double start_time;
	int tt;

	int flux_switch;
	int prod_switch;
	double flux_us;
	double dt;

	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt;
	model_run_params_in.close();


	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname,
						  ft_initial_profile_fname);

	// raise the flowtube so the downstream boundary is at elevation
	// 100
	//vector<double> zeta_old = ft_test.get_zeta();
	//double zeta_bdry = zeta_old[ ft_test.get_n_nodes()-1 ];
	//double increment = zeta_bdry-ds_elev;
	//double zeta_zero = 100.0;
	//ft_test.raise_zeta_eta_ds_bound(zeta_zero);
	//ds_elev = zeta_zero-increment;
	//cout << "ds_elev: " << ds_elev << endl;

	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;

	// variables for the sediment transport loop
	vector<double> surface_change_rate(ft_test.get_n_nodes(),0.0);
										// the erosion rate of deposition rate


	ft_test.print_ft_vecs_to_screen();
	cout << endl;
	cout << "ds_elev: " << ds_elev << endl;
	vector<double> eta_old;
	vector<double> eta_new;

	// now loop through time
	while(t_ime < end_time)
	{
		tt++;
		t_ime += dt;		// increment the time

		ds_elev -= dt*SS_erate;

		// run flux
		eta_old = ft_test.get_eta();
		ft_test.flux_timestep_elev_bc(dt, flux_us, ds_elev,flux_switch, prod_switch,
							surface_change_rate);
		eta_new = ft_test.get_eta();

		if (tt%100000 == 0)
		{
			cout << "time is: " << t_ime << endl;
		}
	}
	for (int i = 0; i<ft_test.get_n_nodes(); i++)
	{
		cout << "prod["<<i<<"]: " << (eta_old[i]-eta_new[i])/dt << endl;
	}
	cout << endl;


	// now print the result to a profile file
	ofstream profile_out;
	profile_out.precision(10);
	profile_out.open(profile_out_fname.c_str());
	ft_test.export_input_profile(profile_out);
	profile_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function runs a flowtube at a steady erosion rate
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
double ft_steady_erate_flat(double SS_erate,
					double start_zeta, double start_h,
					double end_time,
					string model_run_params_fname,
					string sed_trans_param_fname,
					string ft_parameter_fname,
					string profile_out_fname)
{
	// initialize some variables
	//double erate;			// the erosion rate
	double ds_elev;			// elevation at teh downslope node

	double t_ime;
	double start_time;
	int tt;

	int flux_switch;
	int prod_switch;
	double flux_us;
	double dt;

	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt;
	model_run_params_in.close();

	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname);

	// start with a flat plateau
	ft_test.set_const_zeta_eta_h(start_zeta, start_h);

	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;

	// variables for the sediment transport loop
	vector<double> surface_change_rate(ft_test.get_n_nodes(),0.0);
										// the erosion rate of deposition rate

	// get the starting downslope elevation
	ds_elev = ft_test.get_zeta_ds();
										// at the surface
	vector<double> eta_old;
	vector<double> eta_new;
	// now loop through time
	while(t_ime < end_time)
	{
		tt++;
		t_ime += dt;		// increment the time


		ds_elev -= dt*SS_erate;

		// run flux
		eta_old = ft_test.get_eta();
		ft_test.flux_timestep_elev_bc(dt, flux_us, ds_elev,flux_switch, prod_switch,
							surface_change_rate);
		eta_new = ft_test.get_eta();

		if (tt%100000 == 0)
		{
			cout << "time is: " << t_ime << endl;
		}
	}

	for (int i = 0; i<ft_test.get_n_nodes(); i++)
	{
		cout << "prod["<<i<<"]: " << (eta_old[i]-eta_new[i])/dt << endl;
	}
	cout << endl;

	ft_test.print_ft_vecs_to_screen();


	// now print the result to a profile file
	ofstream profile_out;
	profile_out.open(profile_out_fname.c_str());
	profile_out.precision(10);
	ft_test.export_input_profile(profile_out);
	profile_out.close();

	// do one more timestep
	tt++;
	t_ime += dt;		// increment the time
	ds_elev -= dt*SS_erate;
	// run flux
	eta_old = ft_test.get_eta();
	ft_test.flux_timestep_elev_bc(dt, flux_us, ds_elev,flux_switch, prod_switch,
						surface_change_rate);
	eta_new = ft_test.get_eta();
	for (int i = 0; i<ft_test.get_n_nodes(); i++)
	{
		cout << "prod["<<i<<"]: " << (eta_old[i]-eta_new[i])/dt << endl;
	}
	cout << endl;
	ft_test.print_ft_vecs_to_screen();



	return ds_elev;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function initializes a flowtube
// that is used only for sediment transport
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
flowtube flowtube_initializer(string sed_trans_param_fname,
							  string ft_parameter_fname,
							  string ft_initial_profile_fname)
{
	// first the sediment transport parameters for the flowtube#
	double rho_s;			// soil density in kg/m^3
	double rho_r;			// rock density in kg/m^3
	double K_h;
	double S_c;
	double W_0;
	double gamma;
	string temp;

	// load the parameters from a parameter file
	ifstream sed_trans_param_in;
	sed_trans_param_in.open(sed_trans_param_fname.c_str());
	sed_trans_param_in >> temp >> rho_s >> temp >> rho_r
					   >> temp >> K_h >> temp >> S_c
					   >> temp >> W_0 >> temp >>gamma;
	sed_trans_param_in.close();

	cout << "LINE 2548 loaded sed trans param; gamma = " << gamma << endl;

	// now open the flowtube parameter file
	ifstream ft_param_in;
	ft_param_in.open(ft_parameter_fname.c_str());

	// now open the flowtube initial profile file
	ifstream ft_initial_profile_in;
	ft_initial_profile_in.open(ft_initial_profile_fname.c_str());

	cout << "LINE 2558 loaded ft parameters as well as initial profile" << endl;

	// initialize the flowtube
	flowtube ft_test(ft_param_in);	// this sets up the s locations
									// and the areas
	// load the transport parameters
	ft_test.set_transport_params(S_c,K_h,W_0,gamma,rho_s,rho_r);

	cout << "LINE 2566 set the parameters " << endl;

	// now load the profile
	ft_test.load_profile(ft_initial_profile_in);
							// and this loads a profile of zeta and eta

	ft_param_in.close();
	ft_initial_profile_in.close();

	return ft_test;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function initializes a flowtube
// that is used only for sediment transport
// in his function the zeta, eta, and h data members are left empty
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
flowtube flowtube_initializer(string sed_trans_param_fname,
							  string ft_parameter_fname)
{
	// first the sediment transport parameters for the flowtube#
	double rho_s;			// soil density in kg/m^3
	double rho_r;			// rock density in kg/m^3
	double K_h;
	double S_c;
	double W_0;
	double gamma;
	string temp;

	// load the parameters from a parameter file
	ifstream sed_trans_param_in;
	sed_trans_param_in.open(sed_trans_param_fname.c_str());
	sed_trans_param_in >> temp >> rho_s >> temp >> rho_r
					   >> temp >> K_h >> temp >> S_c
					   >> temp >> W_0 >> temp >>gamma;
	sed_trans_param_in.close();

	// now open the flowtube parameter file
	ifstream ft_param_in;
	ft_param_in.open(ft_parameter_fname.c_str());

	// initialize the flowtube
	flowtube ft_test(ft_param_in);	// this sets up the s locations
									// and the areas
	// load the transport parameters
	ft_test.set_transport_params(S_c,K_h,W_0,gamma,rho_s,rho_r);

	ft_param_in.close();

	return ft_test;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes in the names of five files:
// 1)a measured zeta file and
// 2)a measured h file
// 3)a list of final zeta from a model
// 4) a list of final h from a model
// 5) a steady profile file
// it then loops through all of the modeled values calcualting the
// RMSE between the measured topography and soil thickness and the
// predicted topography and soil thickness
// the steady profile file is just used to get teh s locations
// of the nodes
// the measured zeta and h files are in the format:
// s_loc1 zeta1
// s_loc2 zeta2
// etc
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void compare_model_with_measured(string measured_zeta_fname,
					string measured_h_fname,
					string final_zeta_fname,
					string final_h_fname,
					string steady_modeled_fname,
					string RMSE_fname)
{
	int n_nodes;		// the number of nodes in the modeled
						// h and zeta vectors
	// root mean square error
	vector<double> RMSE_zeta;
	vector<double> RMSE_h;

	// open the measured zeta file
	ifstream zeta_meas_in;
	zeta_meas_in.open(measured_zeta_fname.c_str());

	// read in the data
	vector<double> meas_zeta_s_loc;
	vector<double> meas_zeta;
	double temp_s_loc;
	double temp_zeta;
	while(zeta_meas_in >> temp_s_loc >> temp_zeta)
	{
		meas_zeta_s_loc.push_back(temp_s_loc);
		meas_zeta.push_back(temp_zeta);
	}

	// close the measured zeta file
	zeta_meas_in.close();

	// now do the same for the h file
	ifstream h_meas_in;
	h_meas_in.open(measured_h_fname.c_str());

	// read in the data
	vector<double> meas_h_s_loc;
	vector<double> meas_h;
	double temp_h;
	while(h_meas_in >> temp_s_loc >> temp_h)
	{
		meas_h_s_loc.push_back(temp_s_loc);
		meas_h.push_back(temp_h);
	}

	// close the measured h file
	h_meas_in.close();

	// now get the modeled s locations from the initial file
	ifstream sm_in;
	sm_in.open(steady_modeled_fname.c_str());
	vector<double> modeled_s_loc;

	while (sm_in >> temp_s_loc >> temp_h >> temp_h >> temp_h >> temp_h)
	{
		modeled_s_loc.push_back(temp_s_loc);
	}
	sm_in.close();
	n_nodes = modeled_s_loc.size();



	// first we complete all steps for zeta
	// first set up the interpolation
	vector<int> zeta_ds_interp_node_num;
	vector<int> zeta_us_interp_node_num;
	vector<double> zeta_interpolated_distance_fraction;
	linear_interpolation_setup(meas_zeta_s_loc,
						modeled_s_loc,
						zeta_ds_interp_node_num,
						zeta_us_interp_node_num,
						zeta_interpolated_distance_fraction);

	// next open the modeled_zeta file
	ifstream modeled_zeta_in;
	modeled_zeta_in.open(final_zeta_fname.c_str());

	// the final zeta and h files have the format:
	// run_number t_ime zeta_1 zeta_2 ... zeta_n_nodes
	// so each model run has a row
	// reading the files we ignore the first two data elements
	double temp1,temp2;
	vector<double> modeled_zeta(n_nodes);
	vector<double> interpolated_modeled_zeta;
	while (modeled_zeta_in >> temp1 >> temp2)
	{
		//cout << "run number is: " << temp1 << endl;

		// loop through a row collecting zeta values
		for(int i = 0; i<n_nodes; i++)
		{
			modeled_zeta_in >> temp2;
			modeled_zeta[i] = temp2;
		}
		//cout << modeled_zeta.size();

		// now get the interpolated modeled vector
		interpolated_modeled_zeta =
				linear_interpolate_vector(
							zeta_ds_interp_node_num,
							zeta_us_interp_node_num,
							zeta_interpolated_distance_fraction,
							modeled_zeta);

		// now calculate the RMSE for zeta
		RMSE_zeta.push_back(RMSE_mean(meas_zeta,
						interpolated_modeled_zeta));

	}

	modeled_zeta_in.close();

	// now for h
	// first set up the interpolation
	vector<int> h_ds_interp_node_num;
	vector<int> h_us_interp_node_num;
	vector<double> h_interpolated_distance_fraction;
	linear_interpolation_setup(meas_h_s_loc,
						modeled_s_loc,
						h_ds_interp_node_num,
						h_us_interp_node_num,
						h_interpolated_distance_fraction);

	// next open the modeled_zeta file
	ifstream modeled_h_in;
	modeled_h_in.open(final_h_fname.c_str());

	// the final zeta and h files have the format:
	// run_number t_ime zeta_1 zeta_2 ... zeta_n_nodes
	// so each model run has a row
	// reading the files we ignore the first two data elements
	vector<double> modeled_h(n_nodes);
	vector<double> interpolated_modeled_h;
	while (modeled_h_in >> temp1 >> temp2)
	{
		//cout << "run number is: " << temp1 << endl;

		// loop through a row collecting zeta values
		for(int i = 0; i<n_nodes; i++)
		{
			modeled_h_in >> temp2;
			modeled_h[i] = temp2;
		}

		// now get the interpolated modeled vector
		interpolated_modeled_h =
				linear_interpolate_vector(
							h_ds_interp_node_num,
							h_us_interp_node_num,
							h_interpolated_distance_fraction,
							modeled_h);

		// now calculate the RMSE for zeta
		RMSE_h.push_back(RMSE(meas_h,
						interpolated_modeled_h));

	}

	modeled_h_in.close();

	// now print to file
	int n_runs;
	// now prepare the RMSE file
	ofstream RMSE_out;
	RMSE_out.open(RMSE_fname.c_str());
	n_runs = RMSE_h.size();
	for(int i = 0; i<n_runs; i++)
	{
		RMSE_out << i+1 << " " << RMSE_zeta[i] << " " << RMSE_h[i] << endl;
	}
	RMSE_out.close();



}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function extracts a measured zeta and h record from
// final zeta and eta files
// false zeta and h are teh virtual hillslope final values
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void create_meas_zeta_h_from_modeled(int run_number,
					string steady_modeled_fname,
					string zeta_final_fname,
					string h_final_fname,
					string false_zeta_meas_fname,
					string false_h_meas_fname)
{
	int rn_counter;			// a counter used for examining
							// the final zeta and h files

	// get the s locations
	double temp_s_loc;
	double temp;
	int n_nodes;
	vector<double> modeled_s_loc;
	ifstream sm_in;
	sm_in.open(steady_modeled_fname.c_str());
	while (sm_in >> temp_s_loc >> temp >> temp >> temp >> temp)
	{
		modeled_s_loc.push_back(temp_s_loc);
	}
	sm_in.close();
	n_nodes = modeled_s_loc.size();

	// now load zeta
	vector<double> false_zeta_meas(n_nodes);
	ifstream zeta_in;
	zeta_in.open(zeta_final_fname.c_str());
	cout << "zeta_fname is: " << zeta_final_fname << endl;
	rn_counter = 0;
	while(rn_counter < run_number)
	{

		rn_counter++;
		//cout << "run is: " << rn_counter << endl;
		zeta_in >> temp >> temp;
		for (int i = 0; i<n_nodes; i++)
		{
			zeta_in >> temp;
			//cout << " zeta["<<i<<"] is: " << temp << endl;
			false_zeta_meas[i] = temp;
		}
	}
	zeta_in.close();

	// now h
	vector<double> false_h_meas(n_nodes);
	ifstream h_in;
	h_in.open(h_final_fname.c_str());
	rn_counter = 0;
	while(rn_counter < run_number)
	{
		rn_counter++;
		h_in >> temp >> temp;
		for (int i = 0; i<n_nodes; i++)
		{
			h_in >> temp;
			false_h_meas[i] = temp;
		}
	}
	h_in.close();

	// now print the two files
	ofstream zeta_out;
	zeta_out.open(false_zeta_meas_fname.c_str());
	for (int i = 0; i<n_nodes; i++)
	{
		zeta_out << modeled_s_loc[i] << " " << false_zeta_meas[i] << endl;
	}
	zeta_out.close();

	ofstream h_out;
	h_out.open(false_h_meas_fname.c_str());
	for (int i = 0; i<n_nodes; i++)
	{
		h_out << modeled_s_loc[i] << " " << false_h_meas[i] << endl;
	}
	h_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function sets up the interpolation by finding the locations
// of the interpolated nodes vs the predicted nodes.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void linear_interpolation_setup(vector<double> interp_x_loc,
					vector<double> modeled_x_loc,
					vector<int>& ds_interp_node_num,
					vector<int>& us_interp_node_num,
					vector<double>& interpolated_distance_fraction)
{
	int n_nodes = modeled_x_loc.size();
	int n_interp_nodes = interp_x_loc.size();		// number of nodes in the vector to be
													// interpolated
	int h_node_number = 0;				// the node number of the downslope distance
	int interp_node_number = 0;			// node number of the interpolated vector
	vector<int> temp_ds_node_num;		// a temporary vector for storing distance_h nodes
										// adjacent to interpolated nodes (downslope).
	vector<int> temp_us_node_num;		// a temporary vector for storing distance_h nodes
										// adjacent to interpolated nodes (upslope).
	vector<double> temp_idf;			// a vector that holds the farction of distance
										// that an interpolated value is between nodes
	double idf;
	// loop through interpolating vector. The loop stops when it runs out of either
	// interpolating nodes or distance_h_nodes
	while (h_node_number < n_nodes-1 && interp_node_number <n_interp_nodes)
	{
		// if the interpolated node lies between distance_h nodes...
		// the 0.00001 is due to roundoff error
		if (modeled_x_loc[h_node_number+1]+0.00001 > interp_x_loc[interp_node_number] &&
		    modeled_x_loc[h_node_number] <= interp_x_loc[interp_node_number])
		{
			// store the index of the adjacent nodes
			temp_ds_node_num.push_back(h_node_number+1);
			temp_us_node_num.push_back(h_node_number);

			// get the fractional distance of teh interpolated value
			idf = (interp_x_loc[interp_node_number]-modeled_x_loc[h_node_number])/
				  (modeled_x_loc[h_node_number+1]-modeled_x_loc[h_node_number]);

			temp_idf.push_back(idf);
			interp_node_number++;		// and increment the interpolating node number
		}
		else
		{
			h_node_number++;				// now increment the distance_h node number
											// the else statement is necessary for the case
											// when two interpolated node lie between the
											// same two downslope distance nodes.
		}
	}

	// once finished, update the intterpolating reverernce vectors
	ds_interp_node_num = temp_ds_node_num;
	us_interp_node_num = temp_us_node_num;
	interpolated_distance_fraction = temp_idf;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// interpolating function:
// takes information generated using the initialize_interpolation_nodes
// function and returns the interpolaed values of h
// using linear interpolation
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<double> linear_interpolate_vector(vector<int> ds_interp_node_num,
					vector<int> us_interp_node_num,
					vector<double> interpolated_fractional_distance,
					vector<double> modeled_data)
{
	vector<double> temp_interp_h;
	int n_interp_nodes = ds_interp_node_num.size();
	for(int i = 0; i<n_interp_nodes; i++)

	{
		temp_interp_h.push_back( (modeled_data[ ds_interp_node_num[i] ]
		                    -modeled_data[ us_interp_node_num[i] ])
		                    *interpolated_fractional_distance[i]
		                    +modeled_data[ us_interp_node_num[i] ] );
	}
	return temp_interp_h;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//	virtual_erate_history_analysis
// this function caltucates statistics that are used to examing
// if a virtual hillslope can be recreated
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void virtual_erate_history_analysis(int start_datafile,
					int end_datafile, int n_runs_per_datafile,
					int virtual_hs_datfile,
					int virtual_hs_runnumber,
					double t_spacing,
					vector<double> erate_RMSE_efolding_values,
					vector<double> window_times,
					string steady_in_fname,
					string virtual_analysis_out_fname)
{
	int n_efolding_values = erate_RMSE_efolding_values.size();
	int n_window_times = window_times.size();

	int start_run = 1;
	int end_run = n_runs_per_datafile;

	// the strings holding a number of filenames
	string begin_name;
	string end_name;
	string steady_modeled_fname;
	string zeta_final_fname;
	string h_final_fname;
	string file_number;
	string RMSE_fname;
	string erf_datfile_fname;
	string window_2p_out_fname;
	string window_5p_out_fname;
	string window_10p_out_fname;
	string RMSE_erate_fname;

	double temp;

	string virtual_zeta_meas_fname = "zeta_virtual.zdat";
	string virtual_h_meas_fname = "h_virtual.zdat";

	// first get the zeta and h values for the virtual hillslope
	steady_modeled_fname = "steady.sm";

	// get the final zeta fname
	begin_name = "zeta.";
	end_name = ".zdat";
	file_number = itoa(virtual_hs_datfile);
	zeta_final_fname = begin_name+file_number+end_name;

	// get the final zeta fname
	begin_name = "h.";
	end_name = ".hdat";
	h_final_fname = begin_name+file_number+end_name;

	cout << "zeta_fname: " << zeta_final_fname
		 << " and h final fname: " << h_final_fname << endl;

	// now create zeta and eta datafiles that will be read by
	// other functions
	create_meas_zeta_h_from_modeled(virtual_hs_runnumber,
					steady_modeled_fname,
					zeta_final_fname,
					h_final_fname,
					virtual_zeta_meas_fname,
					virtual_h_meas_fname);

	// get the rates at regular time intervasl of the virtual erosion
	// history
	// start by opening the erf datfile
	begin_name = "erf_df.";
	end_name = ".datfile";
	erf_datfile_fname = begin_name + file_number+end_name;
	cout << "erf datfile fname: " << erf_datfile_fname;
	ifstream erf_datfile;
	erf_datfile.open(erf_datfile_fname.c_str());

	vector<double> compare_reg_times;
	vector<double> compare_reg_rates;
	generate_regularly_spaced_erate_hist(
					t_spacing,
					virtual_hs_runnumber,
					erf_datfile,
					compare_reg_times,
					compare_reg_rates);
	erf_datfile.close();


	// open the outfiles and print the intitial line of data
	ofstream virtual_analysis_out;
	virtual_analysis_out.open(virtual_analysis_out_fname.c_str());

	// the columns are:
	// 1 data file number
	// 2 run number
	// 3 RMSE zeta
	// 4 RMSE h
	// 5 weighted RMSE
	// 6-(6+3*n_windows) the window times, repeated
	//       for 2%, 5%, and 10% windows
	// then n efolding times for RMSE of history
	virtual_analysis_out << "-99 -99 -99 -99 -99";
	for (int win =0; win<3; win++)
	{
		for (int wn = 0; wn<n_window_times;wn++)
		{
			virtual_analysis_out << " " <<window_times[wn];
		}
	}
	for (int i = 0; i< n_efolding_values; i++)
	{
		virtual_analysis_out << " " <<erate_RMSE_efolding_values[i];
	}
	virtual_analysis_out << endl;

	// now loop through the data files
	for (int i = start_datafile; i<= end_datafile; i++)
	{

		// get the final zeta fname
		begin_name = "zeta.";
		end_name = ".zdat";
		file_number = itoa(i);
		zeta_final_fname = begin_name+file_number+end_name;

		// get the final zeta fname
		begin_name = "h.";
		end_name = ".hdat";
		h_final_fname = begin_name+file_number+end_name;

		// get the erate history datafile
		begin_name = "erf_df.";
		end_name = ".datfile";
		erf_datfile_fname = begin_name + file_number+end_name;

		// now create temp datafiles
		// NOTE: the code is structured so that each subroutine
		// creates a datafile and then the results from this
		// datafile are read from these files
		// these files are thus temporary holding areas
		// for data
		begin_name = "RMSEtemp.";
		end_name = ".RMSE";
		RMSE_fname = begin_name+file_number+end_name;

		begin_name = "w2ptemp.";
		end_name = ".wdat";
		window_2p_out_fname = begin_name+file_number+end_name;

		begin_name = "w5ptemp.";
		end_name = ".wdat";
		window_5p_out_fname = begin_name+file_number+end_name;

		begin_name = "w10ptemp.";
		end_name = ".wdat";
		window_10p_out_fname = begin_name+file_number+end_name;

		begin_name = "erateRMSEtemp.";
		end_name = ".RMSE";
		RMSE_erate_fname = begin_name+file_number+end_name;

		compare_model_with_measured(virtual_zeta_meas_fname,
					virtual_h_meas_fname,
					zeta_final_fname,
					h_final_fname,
					steady_modeled_fname,
					RMSE_fname);

		compare_erate_windows(t_spacing,
					compare_reg_times,
					compare_reg_rates,
					start_run,
					end_run,
					window_times,
					window_2p_out_fname,
					window_5p_out_fname,
					window_10p_out_fname,
					erf_datfile_fname);

		compare_single_erate_hist_with_others(
					t_spacing,
					compare_reg_times,
					compare_reg_rates,
					start_run,
					end_run,
					erate_RMSE_efolding_values,
					erf_datfile_fname,
					RMSE_erate_fname);

		// now read the temporary data files, printing to
		// the analysis data file
		ifstream RMSE_in,window_2p_in,window_5p_in, window_10p_in,
				 RMSE_erate_in;
		RMSE_in.open(RMSE_fname.c_str());
		window_2p_in.open(window_2p_out_fname.c_str());
		window_5p_in.open(window_5p_out_fname.c_str());
		window_10p_in.open(window_10p_out_fname.c_str());
		RMSE_erate_in.open(RMSE_erate_fname.c_str());
		double RMSE_zeta, RMSE_h;
		//for (int rn = 1; rn<= n_runs_per_datafile; rn++)
		while (RMSE_in >> temp)
		{
			virtual_analysis_out << i << " " << temp;
			RMSE_in >> RMSE_zeta >> RMSE_h;
			virtual_analysis_out << " " << RMSE_zeta << " " << RMSE_h
								<< " " << (1.0/5.0)*(4.0*RMSE_h+RMSE_zeta);

			window_2p_in >> temp;
			for(int wn = 0; wn<n_window_times; wn++)
			{
				window_2p_in >> temp;
				virtual_analysis_out << " " << temp;
			}


			window_5p_in >> temp;
			for(int wn = 0; wn<n_window_times; wn++)
			{
				window_5p_in >> temp;
				virtual_analysis_out << " " << temp;
			}

			window_10p_in >> temp;
			for(int wn = 0; wn<n_window_times; wn++)
			{
				window_10p_in >> temp;
				virtual_analysis_out << " " << temp;
			}

			RMSE_erate_in >> temp;
			for(int wn = 0; wn<n_efolding_values; wn++)
			{
				RMSE_erate_in >> temp;
				virtual_analysis_out << " " << temp;
			}

			virtual_analysis_out << endl;

		}
		RMSE_in.close();
		window_2p_in.close();
		window_5p_in.close();
		window_10p_in.close();
		RMSE_erate_in.close();
	}
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//	virtual_erate_history_analysis
// this function caltucates statistics that are used to examing
// if a virtual hillslope can be recreated
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void virtual_erate_history_analysis(int start_datafile,
					int end_datafile, int n_runs_per_datafile,
					int virtual_hs_datfile,
					int virtual_hs_runnumber,
					double t_spacing,
					vector<double> erate_RMSE_efolding_values,
					vector<double> window_times,
					string steady_in_fname,
					string virtual_analysis_out_fname,
					string pathname)
{
	int n_efolding_values = erate_RMSE_efolding_values.size();
	int n_window_times = window_times.size();

	int start_run = 1;
	int end_run = n_runs_per_datafile;

	// the strings holding a number of filenames
	string begin_name;
	string end_name;
	string steady_modeled_fname;
	string zeta_final_fname;
	string h_final_fname;
	string file_number;
	string RMSE_fname;
	string erf_datfile_fname;
	string window_2p_out_fname;
	string window_5p_out_fname;
	string window_10p_out_fname;
	string RMSE_erate_fname;

	double temp;

	string virtual_zeta_meas_fname = pathname;
	virtual_zeta_meas_fname+="zeta_virtual.zdat";
	string virtual_h_meas_fname = pathname;
	virtual_h_meas_fname+="h_virtual.zdat";

	// first get the zeta and h values for the virtual hillslope
	steady_modeled_fname = "steady.sm";

	// get the final zeta fname
	begin_name = pathname+"zeta.";
	end_name = ".zdat";
	file_number = itoa(virtual_hs_datfile);
	zeta_final_fname = begin_name+file_number+end_name;

	// get the final zeta fname
	begin_name = pathname+"h.";
	end_name = ".hdat";
	h_final_fname = begin_name+file_number+end_name;

	//cout << "zeta_fname: " << zeta_final_fname
	//	 << " and h final fname: " << h_final_fname << endl
	//	 << " virt runnumber: " <<  virtual_hs_runnumber << endl;

	// now create zeta and eta datafiles that will be read by
	// other functions
	create_meas_zeta_h_from_modeled(virtual_hs_runnumber,
					steady_modeled_fname,
					zeta_final_fname,
					h_final_fname,
					virtual_zeta_meas_fname,
					virtual_h_meas_fname);

	// get the rates at regular time intervasl of the virtual erosion
	// history
	// start by opening the erf datfile
	begin_name = pathname;
	begin_name += "erf_df.";
	end_name = ".datfile";
	erf_datfile_fname = begin_name + file_number+end_name;
	cout << "erf datfile fname: " << erf_datfile_fname;
	ifstream erf_datfile;
	erf_datfile.open(erf_datfile_fname.c_str());

	vector<double> compare_reg_times;
	vector<double> compare_reg_rates;
	generate_regularly_spaced_erate_hist(
					t_spacing,
					virtual_hs_runnumber,
					erf_datfile,
					compare_reg_times,
					compare_reg_rates);
	erf_datfile.close();

	//int crt = compare_reg_times.size();
	//for (int ii = 0; ii<crt; ii++)
	//{
	//	cout << "time: " << compare_reg_times[ii] << " " << compare_reg_rates[ii] << endl;
	//}

	// open the outfiles and print the intitial line of data
	ofstream virtual_analysis_out;
	virtual_analysis_out.open(virtual_analysis_out_fname.c_str());

	// the columns are:
	// 1 data file number
	// 2 run number
	// 3 RMSE zeta
	// 4 RMSE h
	// 5 weighted RMSE
	// 6-(6+3*n_windows) the window times, repeated
	//       for 2%, 5%, and 10% windows
	// then n efolding times for RMSE of history
	virtual_analysis_out << "-99 -99 -99 -99 -99";
	for (int win =0; win<3; win++)
	{
		for (int wn = 0; wn<n_window_times;wn++)
		{
			virtual_analysis_out << " " <<window_times[wn];
		}
	}
	for (int i = 0; i< n_efolding_values; i++)
	{
		virtual_analysis_out << " " <<erate_RMSE_efolding_values[i];
	}
	virtual_analysis_out << endl;

	// now loop through the data files
	for (int i = start_datafile; i<= end_datafile; i++)
	{

		// get the final zeta fname
		begin_name = pathname+"zeta.";
		end_name = ".zdat";
		file_number = itoa(i);
		zeta_final_fname = begin_name+file_number+end_name;
		cout << endl << "data filename is: " << zeta_final_fname << endl;

		// get the final zeta fname
		begin_name = pathname+"h.";
		end_name = ".hdat";
		h_final_fname = begin_name+file_number+end_name;

		// get the erate history datafile
		begin_name = pathname+"erf_df.";
		end_name = ".datfile";
		erf_datfile_fname = begin_name + file_number+end_name;

		// now create temp datafiles
		// NOTE: the code is structured so that each subroutine
		// creates a datafile and then the results from this
		// datafile are read from these files
		// these files are thus temporary holding areas
		// for data
		begin_name = pathname+"RMSEtemp.";
		end_name = ".RMSE";
		RMSE_fname = begin_name+file_number+end_name;

		begin_name = pathname+"w2ptemp.";
		end_name = ".wdat";
		window_2p_out_fname = begin_name+file_number+end_name;

		begin_name = pathname+"w5ptemp.";
		end_name = ".wdat";
		window_5p_out_fname = begin_name+file_number+end_name;

		begin_name = pathname+"w10ptemp.";
		end_name = ".wdat";
		window_10p_out_fname = begin_name+file_number+end_name;

		begin_name = pathname+"erateRMSEtemp.";
		end_name = ".RMSE";
		RMSE_erate_fname = begin_name+file_number+end_name;

		//cout << "LINE 2956 getting RMSE" << endl;
		compare_model_with_measured(virtual_zeta_meas_fname,
					virtual_h_meas_fname,
					zeta_final_fname,
					h_final_fname,
					steady_modeled_fname,
					RMSE_fname);

		//cout << "LINE 2964 compare erate windows" << endl;
		//cout << "start run: " << start_run << " and end run: " << end_run << endl;
		compare_erate_windows(t_spacing,
					compare_reg_times,
					compare_reg_rates,
					start_run,
					end_run,
					window_times,
					window_2p_out_fname,
					window_5p_out_fname,
					window_10p_out_fname,
					erf_datfile_fname);

		//cout << "LINE 2976 compare RMSE" << endl;
		compare_single_erate_hist_with_others(
					t_spacing,
					compare_reg_times,
					compare_reg_rates,
					start_run,
					end_run,
					erate_RMSE_efolding_values,
					erf_datfile_fname,
					RMSE_erate_fname);
		cout << "LINE 2986 ending RMSE" << endl;

		// now read the temporary data files, printing to
		// the analysis data file
		ifstream RMSE_in,window_2p_in,window_5p_in, window_10p_in,
				 RMSE_erate_in;
		RMSE_in.open(RMSE_fname.c_str());
		window_2p_in.open(window_2p_out_fname.c_str());
		window_5p_in.open(window_5p_out_fname.c_str());
		window_10p_in.open(window_10p_out_fname.c_str());
		RMSE_erate_in.open(RMSE_erate_fname.c_str());
		double RMSE_zeta, RMSE_h;
		//for (int rn = 1; rn<= n_runs_per_datafile; rn++)
		window_2p_in >> temp;
		for(int wn = 0; wn<n_window_times; wn++)
		{
			window_2p_in >> temp;
		}
		window_5p_in >> temp;
		for(int wn = 0; wn<n_window_times; wn++)
		{
			window_5p_in >> temp;
		}
		window_10p_in >> temp;
		for(int wn = 0; wn<n_window_times; wn++)
		{
			window_10p_in >> temp;
		}
		RMSE_erate_in >> temp;
		for(int wn = 0; wn<n_efolding_values; wn++)
		{
			RMSE_erate_in  >> temp;
		}
		while (RMSE_in >> temp)
		{
			virtual_analysis_out << i << " " << temp;
			RMSE_in >> RMSE_zeta >> RMSE_h;
			virtual_analysis_out << " " << RMSE_zeta << " " << RMSE_h
								<< " " << (1.0/5.0)*(4.0*RMSE_h+RMSE_zeta);

			window_2p_in >> temp;
			for(int wn = 0; wn<n_window_times; wn++)
			{
				window_2p_in >> temp;
				virtual_analysis_out << " " << temp;
			}


			window_5p_in >> temp;
			for(int wn = 0; wn<n_window_times; wn++)
			{
				window_5p_in >> temp;
				virtual_analysis_out << " " << temp;
			}

			window_10p_in >> temp;
			for(int wn = 0; wn<n_window_times; wn++)
			{
				window_10p_in >> temp;
				virtual_analysis_out << " " << temp;
			}

			RMSE_erate_in >> temp;
			for(int wn = 0; wn<n_efolding_values; wn++)
			{
				RMSE_erate_in  >> temp;
				virtual_analysis_out << " " << temp;
			}

			virtual_analysis_out << endl;

		}
		RMSE_in.close();
		window_2p_in.close();
		window_5p_in.close();
		window_10p_in.close();
		RMSE_erate_in.close();
	}

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//	virtual_erate_history_analysis
// this function caltucates statistics that are used to examing
// if a virtual hillslope can be recreated
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void virtual_erate_history_analysis(int start_datafile,
					int end_datafile, int n_runs_per_datafile,
					int virtual_hs_datfile,
					int virtual_hs_runnumber,
					double t_spacing,
					vector<double> erate_RMSE_efolding_values,
					string steady_in_fname,
					string virtual_analysis_out_fname,
					string pathname)
{
	int n_efolding_values = erate_RMSE_efolding_values.size();

	int start_run = 1;
	int end_run = n_runs_per_datafile;

	// the strings holding a number of filenames
	string begin_name;
	string end_name;
	string steady_modeled_fname;
	string zeta_final_fname;
	string h_final_fname;
	string file_number;
	string RMSE_fname;
	string erf_datfile_fname;
	string RMSE_erate_fname;

	double temp;

	string virtual_zeta_meas_fname = pathname;
	virtual_zeta_meas_fname+="zeta_virtual.zdat";
	string virtual_h_meas_fname = pathname;
	virtual_h_meas_fname+="h_virtual.zdat";

	// first get the zeta and h values for the virtual hillslope
	steady_modeled_fname = "steady.sm";

	// get the final zeta fname
	begin_name = pathname+"zeta.";
	end_name = ".zdat";
	file_number = itoa(virtual_hs_datfile);
	zeta_final_fname = begin_name+file_number+end_name;

	// get the final zeta fname
	begin_name = pathname+"h.";
	end_name = ".hdat";
	h_final_fname = begin_name+file_number+end_name;

	//cout << "zeta_fname: " << zeta_final_fname
	//	 << " and h final fname: " << h_final_fname << endl
	//	 << " virt runnumber: " <<  virtual_hs_runnumber << endl;

	// now create zeta and eta datafiles that will be read by
	// other functions
	create_meas_zeta_h_from_modeled(virtual_hs_runnumber,
					steady_modeled_fname,
					zeta_final_fname,
					h_final_fname,
					virtual_zeta_meas_fname,
					virtual_h_meas_fname);

	// get the rates at regular time intervasl of the virtual erosion
	// history
	// start by opening the erf datfile
	begin_name = pathname;
	begin_name += "erf_df.";
	end_name = ".datfile";
	erf_datfile_fname = begin_name + file_number+end_name;
	cout << "erf datfile fname: " << erf_datfile_fname;
	ifstream erf_datfile;
	erf_datfile.open(erf_datfile_fname.c_str());

	vector<double> compare_reg_times;
	vector<double> compare_reg_rates;
	generate_regularly_spaced_erate_hist(
					t_spacing,
					virtual_hs_runnumber,
					erf_datfile,
					compare_reg_times,
					compare_reg_rates);
	erf_datfile.close();

	//int crt = compare_reg_times.size();
	//for (int ii = 0; ii<crt; ii++)
	//{
	//	cout << "time: " << compare_reg_times[ii] << " " << compare_reg_rates[ii] << endl;
	//}

	// open the outfiles and print the intitial line of data
	ofstream virtual_analysis_out;
	virtual_analysis_out.open(virtual_analysis_out_fname.c_str());

	// the columns are:
	// 1 data file number
	// 2 run number
	// 3 RMSE zeta
	// 4 RMSE h
	// 5 weighted RMSE
	// 6-(6+3*n_windows) the window times, repeated
	//       for 2%, 5%, and 10% windows
	// then n efolding times for RMSE of history
	virtual_analysis_out << "-99 -99 -99 -99 -99";

	for (int i = 0; i< n_efolding_values; i++)
	{
		virtual_analysis_out << " " <<erate_RMSE_efolding_values[i];
	}

	virtual_analysis_out << " -99" << endl;

	// now loop through the data files
	for (int i = start_datafile; i<= end_datafile; i++)
	{

		// get the final zeta fname
		begin_name = pathname+"zeta.";
		end_name = ".zdat";
		file_number = itoa(i);
		zeta_final_fname = begin_name+file_number+end_name;
		cout << endl << "data filename is: " << zeta_final_fname << endl;

		// get the final zeta fname
		begin_name = pathname+"h.";
		end_name = ".hdat";
		h_final_fname = begin_name+file_number+end_name;

		// get the erate history datafile
		begin_name = pathname+"erf_df.";
		end_name = ".datfile";
		erf_datfile_fname = begin_name + file_number+end_name;

		// now create temp datafiles
		// NOTE: the code is structured so that each subroutine
		// creates a datafile and then the results from this
		// datafile are read from these files
		// these files are thus temporary holding areas
		// for data
		begin_name = pathname+"RMSEtemp.";
		end_name = ".RMSE";
		RMSE_fname = begin_name+file_number+end_name;

		begin_name = pathname+"erateRMSEtemp.";
		end_name = ".RMSE";
		RMSE_erate_fname = begin_name+file_number+end_name;

		//cout << "LINE 2956 getting RMSE" << endl;
		compare_model_with_measured(virtual_zeta_meas_fname,
					virtual_h_meas_fname,
					zeta_final_fname,
					h_final_fname,
					steady_modeled_fname,
					RMSE_fname);


		//cout << "LINE 2976 compare RMSE" << endl;
		compare_single_erate_hist_with_others(
					t_spacing,
					compare_reg_times,
					compare_reg_rates,
					start_run,
					end_run,
					erate_RMSE_efolding_values,
					erf_datfile_fname,
					RMSE_erate_fname);
		cout << "LINE 2986 ending RMSE" << endl;

		// now read the temporary data files, printing to
		// the analysis data file
		ifstream RMSE_in,window_2p_in,window_5p_in, window_10p_in,
				 RMSE_erate_in;
		RMSE_in.open(RMSE_fname.c_str());
		RMSE_erate_in.open(RMSE_erate_fname.c_str());
		double RMSE_zeta, RMSE_h;
		//for (int rn = 1; rn<= n_runs_per_datafile; rn++)

		RMSE_erate_in >> temp;
		for(int wn = 0; wn<n_efolding_values; wn++)
		{
			RMSE_erate_in  >> temp;
		}
		while (RMSE_in >> temp)
		{
			virtual_analysis_out << i << " " << temp;
			RMSE_in >> RMSE_zeta >> RMSE_h;
			virtual_analysis_out << " " << RMSE_zeta << " " << RMSE_h
								<< " " << (1.0/5.0)*(4.0*RMSE_h+RMSE_zeta);



			RMSE_erate_in >> temp;
			vector<double> RMSE_erate_temp_vec(n_efolding_values);
			double avg_RMSE = 0;
			for(int wn = 0; wn<n_efolding_values; wn++)
			{
				RMSE_erate_in  >> temp;
				RMSE_erate_temp_vec[wn] = temp;
				virtual_analysis_out << " " << temp;
				avg_RMSE+= temp;
			}
			avg_RMSE = avg_RMSE/double(n_efolding_values);


			virtual_analysis_out << " " << avg_RMSE << endl;

		}
		RMSE_in.close();
		RMSE_erate_in.close();
	}

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// "I'm gonna take a lackadaisical ride on my back in the daycycle"
// -MC Paul Barman

#endif
