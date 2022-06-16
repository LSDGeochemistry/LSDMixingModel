// Part_util

#ifndef Part_UTIL_CPP
#define Part_UTIL_CPP

#include "Part_util.hpp"
#include "flowtube.hpp"
#include "FT_util.hpp"
#include "chronos_particle_info.hpp"
#include "CRN_tParticle_bins.hpp"
#include "mathutil.h"
#include <math.h>
#include <time.h>
#include <vector>
#include <list>
#include <iostream>
using namespace std;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function runs a flowtube at a steady erosion rate, and
// tests the particle inserter
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void part_ft_steady_erate(double SS_erate, double end_time,
					double ds_elev,
					vector<double> surf_erate,
					double insert_interval_time,
					double particle_print_interval,
					double zhe_print_interval,
					double age_print_interval,				//*****
					double eroded_catch_window,
					double max_age,							//*********
					int n_spacings,							//*********
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
					string age_pdf_out_fname,					//*******
					string age_cdf_out_fname,					// ******
					string zeta_out_fname, string eta_out_fname, string h_out_fname)
{
	// initialize some variables
	//double erate;			// the erosion rate

	double t_ime;
	double start_time;
	int tt;					// a counter for the time

	// the parameters for the in situ cosmogenics
	LSDCRNParameters CRNp;

	int flux_switch;				// see flowtube.h: determines flux law
	int prod_switch;				// see flowtube.h: determines soil production law
	double flux_us;					// see flowtube.h: flux from upslope
	double dt;						// time interval

	double start_depth;				// the starting depth in m
									// (LATER) it is set to 60 because
									// this is roughly 3x the e-folding
									// depth of the deepest muonogenic
									// production mechanism at rock density
	double vert_mix_vel;
	double horiz_mix_vel;
	double Omega;					// the activity of particles (that is the
									// proportion of particles that are moving
									// at any given time)

	//double d;						// depth (m)
	//double eff_d;					// effective depth (in g/cm^2)

	//double z_p;
									// the elevation of the introduced particle
	double part_conc;				// particle concentration in particles per kg

	vector<double> old_eta;			// the elevation of the soil-saprolite boudnary from the
									// last timestep
	vector<double> new_eta;			// the elevation of the soil-saprolite boundary from this
									// timestep
	vector< list<LSDCRNParticle> > eroded_bins;
									// a vector of lists of particles eroded from each hillslope
									// node
	vector< list<LSDCRNParticle> > temp_part_bins;
	vector< list<LSDCRNParticle> > particle_bins;
	// get the particle types
	Particle_info pi(particle_list_fname.c_str());
	vector<int> starting_pID = pi.get_type_index();
	int n_ptypes = starting_pID.size();

	ifstream pmfrac_in;
	pmfrac_in.open(particle_mfracs_fname.c_str());
	vector<double> starting_p_mfrac(n_ptypes);
	for (int i = 0; i<n_ptypes; i++)
	{
		pmfrac_in >> starting_p_mfrac[i];
	}
	pmfrac_in.close();

	// now open the zeta, eta, and h outfiles
	ofstream zeta_out,h_out,eta_out;
	zeta_out.open(zeta_out_fname.c_str());
	h_out.open(h_out_fname.c_str());
	eta_out.open(eta_out_fname.c_str());


	ofstream age_cdf_out,age_pdf_out;
	age_cdf_out.open(age_cdf_out_fname.c_str());
	age_pdf_out.open(age_pdf_out_fname.c_str());


	// get the model parameters
	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt;
	model_run_params_in.close();

	// an integer used for printing
	const int part_p_i = int (double(particle_print_interval/dt+0.5) );
	//const int zhe_p_i = int (double(zhe_print_interval/dt+0.5) );
	//const int age_p_i = int (double(age_print_interval/dt+0.5) );
	double next_catch = insert_interval_time;

	// get the parameters for the CRN particles
	ifstream CRN_parameter_in;
	CRN_parameter_in.open(CRN_parameter_fname.c_str());
	CRN_parameter_in >> temp >> start_depth >> temp >> vert_mix_vel
				     >> temp >> horiz_mix_vel >> temp >> Omega
				     >> temp >> part_conc;
	CRN_parameter_in.close();

	//cout << "LINE 121 got LSDCRNParameters " << endl;

	// open the datafile for the particle information
	ofstream particle_out;
	particle_out.open(particle_out_fname.c_str());
	ofstream eroded_particle_out;
	eroded_particle_out.open(eroded_particle_out_fname.c_str());



	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname,
						  ft_initial_profile_fname);


	//cout << "LINE 135 got flowtube " << endl;

	// initialize a CRN_particle_list
	LSDCRNParticle_bins CRN_tpb(ft_test);
	int n_bins = CRN_tpb.get_n_bins();
    vector< list<LSDCRNParticle> > eroded_catcher(n_bins+1);
    vector< list<LSDCRNParticle> > empty_eroded_catcher(n_bins+1);
    list<LSDCRNParticle>::iterator part_iter;	// list iterator

    //cout << "LINE 144 created CRN_tpart_bins " << endl;

	// raise the flowtube so the downstream boundary is at elevation
	// 100
	vector<double> zeta_old = ft_test.get_zeta();
	double zeta_bdry = zeta_old[ ft_test.get_n_nodes()-1 ];
	double increment = zeta_bdry-ds_elev;
	double zeta_zero = 100.0;
	ft_test.raise_zeta_eta_ds_bound(zeta_zero);
	ds_elev = zeta_zero-increment;
	//cout << "ds_elev: " << ds_elev << endl;

	//cout << "LINE 156 raised elevation " << endl;

	double insert_time_clock = 0;
	old_eta = ft_test.get_eta();
	vector<double> Delta_eta = old_eta;
	vector<double> old_bottom_depth = old_eta;
	vector<double> old_zeta = ft_test.get_zeta();
	vector<double> old_h = ft_test.get_h();
	int eta_sz = Delta_eta.size();
	for (int i = 0; i<eta_sz; i++)
	{
		//cout << "LINE 371 h["<<i<<"]: " << old_h[i] << " and zeta: " << old_eta[i] << endl;
		//old_bottom_depth[i] = old_eta[i]-start_depth;
		old_bottom_depth[i] = old_eta[i]-start_depth;
		Delta_eta[i] = start_depth;
	}
	int part_ID_start = 1;

	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;

	// insert initial particles
	part_ID_start =CRN_tpb.insert_particles(ft_test, Delta_eta, old_bottom_depth,
								 part_conc, starting_pID, starting_p_mfrac);


	// print initial conditions
	CRN_tpb.print_particle_stats_soil(t_ime, ft_test, particle_out);
	ft_test.print_h_s(zeta_out);
	ft_test.print_h_s(eta_out);
	ft_test.print_h_s(h_out);
	ft_test.print_zeta(t_ime, zeta_out);
	ft_test.print_eta(t_ime, eta_out);
	ft_test.print_h(t_ime, h_out);

	// variables for the sediment transport loop
	//vector<double> surface_change_rate(ft_test.get_n_nodes(),0.0);
	//									// the erosion rate of deposition rate

	// now loop through time
	int particle_trigger = 1;
	while(t_ime < end_time)
	{
		tt++;
		t_ime += dt;		// increment the time
		insert_time_clock+=dt;	// increment the insert time clock
								// particles are not inserted every timestep
								// rather they are inserted at intervals.
								// when the insert time clock exceeds the insert
								// interval particles are inserted
								// and the insert time clock is reset

		cout << "LINE 211 Part_util.cpp time is: " << t_ime << endl;

		ds_elev -= dt*SS_erate;

		// run flux
		ft_test.flux_timestep_elev_bc(dt, flux_us, ds_elev,flux_switch, prod_switch,
							surf_erate);

		// run particle motion if particles have been inserted.
		// particle_trigger == 1 when particles have been inserted
		if (particle_trigger == 1)
		{
			eroded_bins = CRN_tpb.particle_motion(dt, ft_test,
										Omega, vert_mix_vel,
										horiz_mix_vel,CRN_switch, CRNp);

			// if the time is within the eroded catch window, catch
			// all the particles in the eroded bin
			if (t_ime > next_catch-eroded_catch_window)
			{
				//cout << "LINE 423, time is: " << t_ime << " and next_catch: " << next_catch << endl;
				for(int bn = 0; bn<=n_bins; bn++)
				{
					if (eroded_bins[bn].size()>0)
					{
						part_iter = eroded_bins[bn].begin();
						while(part_iter != eroded_bins[bn].end())
						{
							eroded_catcher[bn].push_back(*part_iter);
							part_iter++;
						}
					}
				}
			}
		}

		// if the time elapsed sinse the last insertion of particles is gereater
		// than or equal to the insertion interval, then get the amount of
		// soil-saprolite boundary lowering that has occurred and insert
		// particles into the insertion zone
		if (insert_time_clock >= insert_interval_time - dt/2)
		{
			// if this is the first time we are inserting particles,
			// set the particle trigger to 1
			if (particle_trigger == 0)
			{
				particle_trigger = 1;
			}
			new_eta = ft_test.get_eta(); 	// get the updated eta

			//cout << "Time: " << t_ime << " Delta eta: " << endl;
			// calculate delta eta
			for(int ii = 0; ii< eta_sz; ii++)
			{
				Delta_eta[ii] = old_eta[ii]-new_eta[ii];
				//cout << "i = " << ii << " Delta_eta: " << Delta_eta[ii] << endl;
			}

			// insert particles into the insertion zone

			part_ID_start =CRN_tpb.insert_particles(ft_test, Delta_eta, old_bottom_depth,
								 part_conc, starting_pID, starting_p_mfrac);

			old_eta = new_eta;				// reset old eta
			insert_time_clock = 0;			// reset the insert time clock

		}

		// print the particle data to file
		int n_depthintervals_soil = 2;
		int n_depthintervals_parent = 3;
		double bottom_depth = 2.0;
		if (particle_trigger == 1 && tt%part_p_i== 0)
		{
			cout << "LINE 289, printing particles, time: " << t_ime << endl;
			//CRN_tpb.print_particle_stats(t_ime, ft_test, particle_out);
			//CRN_tpb.print_particle_stats_soil(t_ime, ft_test, particle_out);
			//CRN_tpb.print_particle_stats_vtk(t_ime, ft_test, vtk_fname);
			//CRN_tpb.cell_printing_vtk(t_ime, ft_test, vtk_cell_fname,
			//							n_depthintervals_soil, n_depthintervals_parent,
			//							bottom_depth);
			//CRN_tpb.cell_and_particle_printing_vtk(t_ime, ft_test, vtk_particle_fname, vtk_cell_fname,
			//										n_depthintervals_soil, n_depthintervals_parent,
			//								bottom_depth);
			int ref_frame_switch = 0;
			CRN_tpb.cell_and_particle_chemistry_printing_vtk(t_ime, ft_test,
											 pi, vtk_particle_fname, vtk_cell_fname,
											 n_depthintervals_soil, n_depthintervals_parent,
											 bottom_depth, ref_frame_switch);

			//CRN_tpb.print_eroded_stats(t_ime,eroded_catcher,ft_test,eroded_particle_out);
			//next_catch+=particle_print_interval;	// reset the next catch time
			//eroded_catcher = empty_eroded_catcher;
												// reset the eroded catcher
		}
/*
		if (particle_trigger == 1 && tt%zhe_p_i == 0)
		{
			cout << "LINE 207, printing zeta, eta and h, time: " << t_ime << endl;
			ft_test.print_zeta(t_ime, zeta_out);
			ft_test.print_eta(t_ime, eta_out);
			ft_test.print_h(t_ime, h_out);
		}


		if (particle_trigger == 1 && tt%age_p_i == 0)
		{
			cout << "LINE 302, printing age data, time: " << t_ime << endl;
			particle_bins = CRN_tpb.get_particle_bins();
			//cout << "LINE 304, got particle bins" << endl;
			double K_times_D = 1e-7;
			double D = 0.0001;
			double sigma = -0.42;

			CRN_tpb.print_age_cdfpdf_bins(t_ime, max_age, n_spacings,
					 					  K_times_D, D, sigma, ft_test,
				     					  age_cdf_out, age_pdf_out);

			CRN_tpb.print_eroded_stats(t_ime,eroded_catcher,ft_test,eroded_particle_out);
			//cout << "LINE 327, printed eroded stats" << endl;
			next_catch+=particle_print_interval;	// reset the next catch time
			//cout << "LINE 329, reset catch interval" << endl;
			eroded_catcher = empty_eroded_catcher;
			particle_bins = temp_part_bins;
			//cout << "LINE 332, reset particle bin vecs" << endl;
			//cout << "LINE 333, done printing age data, time: " << t_ime << endl;
		}

*/

		if (tt%1000 == 0)
		{
			cout << "time is: " << t_ime << endl;

		}
	}


	// now print the result to a profile file
	ofstream profile_out;
	profile_out.open(profile_out_fname.c_str());
	ft_test.export_input_profile(profile_out);
	profile_out.close();
	particle_out.close();
	eroded_particle_out.close();
	zeta_out.close();
	eta_out.close();
	h_out.close();
	age_cdf_out.close();
	age_pdf_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function runs a flowtube at a steady erosion rate, and
// tests the particle inserter
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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
					string zeta_out_fname, string eta_out_fname, string h_out_fname)
{
	// initialize some variables
	//double erate;			// the erosion rate

	double t_ime;
	double start_time;
	int tt;					// a counter for the time

	// the parameters for the in situ cosmogenics
	LSDCRNParameters CRNp;

	int flux_switch;				// see flowtube.h: determines flux law
	int prod_switch;				// see flowtube.h: determines soil production law
	double flux_us;					// see flowtube.h: flux from upslope
	double dt;						// time interval

	double start_depth;				// the starting depth in m
									// (LATER) it is set to 60 because
									// this is roughly 3x the e-folding
									// depth of the deepest muonogenic
									// production mechanism at rock density
	double vert_mix_vel;
	double horiz_mix_vel;
	double Omega;					// the activity of particles (that is the
									// proportion of particles that are moving
									// at any given time)

	//double d;						// depth (m)
	//double eff_d;					// effective depth (in g/cm^2)

	//double z_p;
									// the elevation of the introduced particle
	double part_conc;				// particle concentration in particles per kg

	vector<double> old_eta;			// the elevation of the soil-saprolite boudnary from the
									// last timestep
	vector<double> new_eta;			// the elevation of the soil-saprolite boundary from this
									// timestep
	vector< list<LSDCRNParticle> > eroded_bins;
	vector< list<LSDCRNParticle> > temp_part_bins;
	vector< list<LSDCRNParticle> > particle_bins;
									// a vector of lists of particles eroded from each hillslope

	// get the particle types
	Particle_info pi(particle_list_fname.c_str());
	vector<int> starting_pID = pi.get_type_index();
	int n_ptypes = starting_pID.size();

	ifstream pmfrac_in;
	pmfrac_in.open(particle_mfracs_fname.c_str());
	vector<double> starting_p_mfrac(n_ptypes);
	for (int i = 0; i<n_ptypes; i++)
	{
		pmfrac_in >> starting_p_mfrac[i];
	}
	pmfrac_in.close();


	// now open the zeta, eta, and h outfiles
	ofstream zeta_out,h_out,eta_out;
	zeta_out.open(zeta_out_fname.c_str());
	h_out.open(h_out_fname.c_str());
	eta_out.open(eta_out_fname.c_str());


	ofstream age_cdf_out,age_pdf_out;
	age_cdf_out.open(age_cdf_out_fname.c_str());
	age_pdf_out.open(age_pdf_out_fname.c_str());


	// get the model parameters
	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt;
	model_run_params_in.close();

	// an integer used for printing
	const int part_p_i = int (double(particle_print_interval/dt+0.5) );
	//const int zhe_p_i = int (double(zhe_print_interval/dt+0.5) );
	//const int age_p_i = int (double(age_print_interval/dt+0.5) );
	double next_catch = insert_interval_time;

	// get the parameters for the CRN particles
	ifstream CRN_parameter_in;
	CRN_parameter_in.open(CRN_parameter_fname.c_str());
	CRN_parameter_in >> temp >> start_depth >> temp >> vert_mix_vel
				     >> temp >> horiz_mix_vel >> temp >> Omega
				     >> temp >> part_conc;
	CRN_parameter_in.close();

	// open the datafile for the particle information
	ofstream particle_out;
	particle_out.open(particle_out_fname.c_str());
	ofstream eroded_particle_out;
	eroded_particle_out.open(eroded_particle_out_fname.c_str());

	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname,
						  ft_initial_profile_fname);

	// initialize a CRN_particle_list
	CRN_tParticle_bins CRN_tpb(ft_test);
	int n_bins = CRN_tpb.get_n_bins();
    vector< list<LSDCRNParticle> > eroded_catcher(n_bins+1);
    vector< list<LSDCRNParticle> > empty_eroded_catcher(n_bins+1);
    list<LSDCRNParticle>::iterator part_iter;	// list iterator

	// raise the flowtube so the downstream boundary is at elevation
	// 100
	vector<double> zeta_old = ft_test.get_zeta();
	//double zeta_bdry = zeta_old[ ft_test.get_n_nodes()-1 ];
	//double increment = zeta_bdry-ds_elev;
	double zeta_zero = 100.0;
	ft_test.raise_zeta_eta_ds_bound(zeta_zero);
	//ds_elev = zeta_zero-increment;
	//cout << "ds_elev: " << ds_elev << endl;

	double insert_time_clock = 0;
	old_eta = ft_test.get_eta();
	vector<double> Delta_eta = old_eta;
	vector<double> old_bottom_depth = old_eta;
	vector<double> old_zeta = ft_test.get_zeta();
	vector<double> old_h = ft_test.get_h();
	int eta_sz = Delta_eta.size();
	for (int i = 0; i<eta_sz; i++)
	{
		cout << "LINE 371 h["<<i<<"]: " << old_h[i] << " and zeta: " << old_eta[i] << endl;
		old_bottom_depth[i] = old_eta[i];
		Delta_eta[i] = start_depth;
	}
	int part_ID_start = 1;

	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;

	// insert initial particles
	part_ID_start =CRN_tpb.insert_particles(ft_test, Delta_eta, old_bottom_depth,
								 part_conc, starting_pID, starting_p_mfrac);

	// print initial conditions
	//CRN_tpb.print_particle_stats(t_ime, ft_test, particle_out);
	ft_test.print_h_s(zeta_out);
	ft_test.print_h_s(eta_out);
	ft_test.print_h_s(h_out);
	//ft_test.print_zeta(t_ime, zeta_out);
	//ft_test.print_eta(t_ime, eta_out);
	//ft_test.print_h(t_ime, h_out);

	cout << "LINE 379, n_nodes: " << ft_test.get_n_nodes() << endl;

	// variables for the sediment transport loop
	//vector<double> surface_change_rate(ft_test.get_n_nodes(),0.0);
										// the erosion rate of deposition rate

	// now loop through time
	int particle_trigger = 1;
	while(t_ime < end_time)
	{
		tt++;
		t_ime += dt;		// increment the time
		insert_time_clock+=dt;	// increment the insert time clock
								// particles are not inserted every timestep
								// rather they are inserted at intervals.
								// when the insert time clock exceeds the insert
								// interval particles are inserted
								// and the insert time clock is reset

		//ds_elev -= dt*SS_erate;

		// run flux
		ft_test.flux_timestep_flux_bc(dt, flux_us, SS_flux,flux_switch, prod_switch,
							surf_erate);

		// run particle motion if particles have been inserted.
		// particle_trigger == 1 when particles have been inserted
		if (particle_trigger == 1)
		{
			eroded_bins = CRN_tpb.particle_motion(dt, ft_test,
										Omega, vert_mix_vel,
										horiz_mix_vel, CRN_switch, CRNp);

			// if the time is within the eroded catch window, catch
			// all the particles in the eroded bin
			if (t_ime > next_catch-eroded_catch_window)
			{
				//cout << "LINE 423, time is: " << t_ime << " and next_catch: " << next_catch << endl;
				for(int bn = 0; bn<=n_bins; bn++)
				{
					if (eroded_bins[bn].size()>0)
					{
						part_iter = eroded_bins[bn].begin();
						while(part_iter != eroded_bins[bn].end())
						{
							eroded_catcher[bn].push_back(*part_iter);
							part_iter++;
						}
					}
				}
			}
		}

		// if the time elapsed sinse the last insertion of particles is gereater
		// than or equal to the insertion interval, then get the amount of
		// soil-saprolite boundary lowering that has occurred and insert
		// particles into the insertion zone
		if (insert_time_clock >= insert_interval_time - dt/2)
		{
			// if this is the first time we are inserting particles,
			// set the particle trigger to 1
			if (particle_trigger == 0)
			{
				particle_trigger = 1;
			}
			new_eta = ft_test.get_eta(); 	// get the updated eta

			//cout << "Time: " << t_ime << " Delta eta: " << endl;
			// calculate delta eta
			for(int ii = 0; ii< eta_sz; ii++)
			{
				Delta_eta[ii] = old_eta[ii]-new_eta[ii];
				//cout << "i = " << ii << " Delta_eta: " << Delta_eta[ii] << endl;
			}

			// insert particles into the insertion zone

			part_ID_start =CRN_tpb.insert_particles(ft_test, Delta_eta, old_bottom_depth,
								 part_conc, starting_pID, starting_p_mfrac);

			old_eta = new_eta;				// reset old eta
			insert_time_clock = 0;			// reset the insert time clock

		}

		// print the particle data to file
		int n_depthintervals_soil = 2;
		int n_depthintervals_parent = 3;
		double bottom_depth = 2.0;
		if (particle_trigger == 1 && tt%part_p_i== 0)
		{
			cout << "LINE 289, printing particles, time: " << t_ime << endl;
			CRN_tpb.print_particle_stats(t_ime, ft_test, particle_out);
			CRN_tpb.print_particle_stats_soil(t_ime, ft_test, particle_out);
			CRN_tpb.print_particle_stats_vtk(t_ime, ft_test, vtk_fname);
			//CRN_tpb.cell_printing_vtk(t_ime, ft_test, vtk_cell_fname,
			//							n_depthintervals_soil, n_depthintervals_parent,
			//							bottom_depth);
			//CRN_tpb.cell_and_particle_printing_vtk(t_ime, ft_test, vtk_particle_fname, vtk_cell_fname,
			//										n_depthintervals_soil, n_depthintervals_parent,
			//								bottom_depth);
			int ref_frame_switch = 0;
			CRN_tpb.cell_and_particle_chemistry_printing_vtk(t_ime, ft_test,
											 pi, vtk_particle_fname, vtk_cell_fname,
											 n_depthintervals_soil, n_depthintervals_parent,
											 bottom_depth, ref_frame_switch);

			//CRN_tpb.print_eroded_stats(t_ime,eroded_catcher,ft_test,eroded_particle_out);
			//next_catch+=particle_print_interval;	// reset the next catch time
			//eroded_catcher = empty_eroded_catcher;
												// reset the eroded catcher
		}
/*
		if (particle_trigger == 1 && tt%zhe_p_i == 0)
		{
			cout << "LINE 207, printing zeta, eta and h, time: " << t_ime << endl;
			ft_test.print_zeta(t_ime, zeta_out);
			ft_test.print_eta(t_ime, eta_out);
			ft_test.print_h(t_ime, h_out);
		}


		if (particle_trigger == 1 && tt%age_p_i == 0)
		{
			cout << "LINE 568, printing age data, time: " << t_ime << endl;
			particle_bins = CRN_tpb.get_particle_bins();
			double K_times_D = 1e-7;
			double D = 0.0001;
			double sigma = -0.42;
			CRN_tpb.print_age_cdfpdf_bulk(t_ime, max_age, n_spacings,
					 					  K_times_D, D, sigma, ft_test,
				     					  age_cdf_out, age_pdf_out);
			CRN_tpb.print_eroded_stats(t_ime,eroded_catcher,ft_test,eroded_particle_out);
			next_catch+=particle_print_interval;	// reset the next catch time
			eroded_catcher = empty_eroded_catcher;
			particle_bins = temp_part_bins;
			cout << "LINE 573, done printing age data, time: " << t_ime << endl;
		}

*/

		if (tt%100000 == 0)
		{
			cout << "time is: " << t_ime << endl;

		}
	}


	// now print the result to a profile file
	ofstream profile_out;
	profile_out.open(profile_out_fname.c_str());
	ft_test.export_input_profile(profile_out);
	profile_out.close();
	particle_out.close();
	eroded_particle_out.close();
	zeta_out.close();
	eta_out.close();
	h_out.close();
	age_cdf_out.close();
	age_pdf_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function runs a flowtube using an erate history file
// it also extracts data from a sampling file and checks the data against the
// sampled data
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void part_ft_erate_from_erh(string run_name, vector<double>& sample_s_locs, vector<double>& d_top, vector<double>& d_bottom,
							vector<double>& meas, vector<double>& unc, vector<double>& modelled, double& MLE,
							vector<double>& c_frac, vector<double>& pf1, vector<double>& pf2,
							vector<double>& pf3, vector<double>& pf4, vector<double>& pf5, vector<double>& C10Be)
{
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=
	// set up the infiles
	string sed_trans_param_ext = ".sed_trans_param.stparam";
	string sed_trans_param_fname = run_name+sed_trans_param_ext;
	string model_run_params_ext = ".model_run.param";
	string model_run_params_fname = run_name+model_run_params_ext;
	string CRN_parameter_ext = ".CRN_trans_param.CRNparam";
	string CRN_parameter_fname = run_name+CRN_parameter_ext;
	string ft_parameter_ext =  ".ft_details.param";
	string ft_parameter_fname = run_name+ft_parameter_ext;
	string profile_in_ext = ".profile.sm";
	string profile_in_fname = run_name+profile_in_ext;
	string particle_types_ext = ".five_comp.clist";
	string particle_types_fname = run_name+particle_types_ext;
	string part_mfracs_ext = ".mfrac.mfrac";
	string part_mfracs_fname = run_name+part_mfracs_ext;
	string profile_enrich_ext = ".profile_enrich.data";
	string profile_enrich_fname = run_name+profile_enrich_ext;
	string erh_ext = ".erh.erhdat";
	string erh_fname = run_name+erh_ext;

	// set up the outfiles
	string profile_out_ext = ".column_out.sm";
	string profile_out_fname = run_name+profile_out_ext;
	string particle_out_ext = ".p_trans_out.pout";
	string particle_out_fname = run_name+particle_out_ext;
	string eroded_pout_ext = ".ep_trans_out.pout";
	string eroded_pout_fname = run_name+eroded_pout_ext;
	string zeta_out_ext = ".zeta_trans.zdat";
	string zeta_out_fname = run_name+zeta_out_ext;
	string h_out_ext = ".h_trans.hdat";
	string h_out_fname = run_name+h_out_ext;
	string eta_out_ext = ".eta_trans.edat";
	string eta_out_fname = run_name+eta_out_ext;
	string age_cdf_out_ext = ".age_cdf.acdf";
	string age_cdf_out_fname = run_name+age_cdf_out_ext;
	string age_pdf_out_ext = ".age_pdf.apdf";
	string age_pdf_out_fname = run_name+age_pdf_out_ext;
	string vtk_cell_ext = ".Cell_data";
	string vtk_cell_fname = run_name+vtk_cell_ext;
	string vtk_particle_ext = ".Particle_data";
	string vtk_particle_fname = run_name+vtk_particle_ext;
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=

	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=
	// Initialize the variables
	double SS_flux;										// flux
	double constant_surface_change_rate;
	vector<double> surf_erate;							// erosion from the surface
														// in m/yr: this is negative for
														// erosion
	int CRN_switch;				// sets the cosmogenics key:
										// 0 == no CRN
										// 1 == all CRN
										// 2 == all, neutron only
										// 3 == 10Be full
										// default == all, neutron only

	double t_ime;
	double start_time;
	int tt;					// a counter for the time

	int flux_switch;				// see flowtube.h: determines flux law
	int prod_switch;				// see flowtube.h: determines soil production law
	double flux_us;					// see flowtube.h: flux from upslope
	double dt;						// time interval

	double start_depth;				// the starting depth in m
									// (LATER) it is set to 60 because
									// this is roughly 3x the e-folding
									// depth of the deepest muonogenic
									// production mechanism at rock density
	double vert_mix_vel;
	double horiz_mix_vel;
	double Omega;					// the activity of particles (that is the
									// proportion of particles that are moving
									// at any given time)

									// the elevation of the introduced particle
	double part_conc;				// particle concentration in particles per kg

	vector<double> old_eta;			// the elevation of the soil-saprolite boudnary from the
									// last timestep
	vector<double> new_eta;			// the elevation of the soil-saprolite boundary from this
									// timestep
	vector< list<LSDCRNParticle> > eroded_bins;
	vector< list<LSDCRNParticle> > temp_part_bins;
	vector< list<LSDCRNParticle> > particle_bins;
									// a vector of lists of particles eroded from each hillslope

	// initial in situ concentrations in atoms per gram
	double C_10Be,C_26Al,C_36Cl,C_14C,C_21Ne,C_3He;

	// set the parameters for production
	// 0 == granger
	// 1 == schaller
	int CRN_muon_param_switch;

	// the parameters for the in situ cosmogenics
	LSDCRNParameters CRNp;

	// note: scaling determined using cosmocalc:
	// Schaller reports pinedale at 42 53 26 N
	// this converts to an inclination of 61.7 degrees
	// Schaller reports elevation of 2298 masl
	// this converts to an atmospheric depth of 781 g/cm2
	// this converts to a dunai scaling of 5.99
	// multiply this by schaller's snow shielding (for nucleonic production
	// of 0.925 one gets a scaling of 5.54
	double single_scaling;

	// paramters for meteoric 10Be
	double M_supply_surface;	// in atoms/(cm^2*yr)
	double k_f10Be;				// in cm^2/g this is
									// an efolding depth of 20 g/cm^2
	double k2_f10Be;
	double chi_f10Be;			// fraction of supply that goes into the
								// shallow meteoric supply
	double deltad;				// in m (this gets converted
								// to cm in LSDParticle.cpp)

	// parameters for dealing with time
	double end_time;
	double insert_interval;
	double particle_printing_interval;
	double zhe_printing_interval;
	double eroded_catch_window;
	double age_printing_interval;
	double max_age;
	int n_spacings;
	int reference_frame_printing_switch; 		// a switch for the vtk files:
												// 0 == printing relative to base level
												// 1 == printing as a function of depth

	// parameters for dealing with the erate history
	// vectors for holding individual erate histories
	int n_ts;
	int curr_erate_index;
	vector<double> erate_times;
	vector<double> erate_rates;
	double erate;
	double ds_elev;
	double insert_interval_time;
	double old_erate;
	double old_erate_time;
	double new_erate;
	double new_erate_time;
	double erate_slope;
	double time_diff;
	double starting_time;

	// some chemistry
	vector<double> enrich_s_loc;
	vector<double> enrich_d_top;
	vector<double> enrich_d_bottom;
	vector<double> enrich_samp;
	vector<double> enrich_1sigma;

	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	//=-=-=-=-=-=-=-=
	// load files
	//=-=-=-=-=-=-=-=
	// load the profile data
	ifstream enrich_in;
	enrich_in.open(profile_enrich_fname.c_str());
	int n_enrich_samples;
	enrich_in >> n_enrich_samples;
	cout << "LINE 900 number_of enrich samples; " << n_enrich_samples << endl;
	double temp_s,temp_enrich,temp_top,temp_bottom, temp_1sigma;
	for (int i = 0; i< n_enrich_samples; i++)
	{
		enrich_in >> temp_s >> temp_top >> temp_bottom >> temp_enrich >> temp_1sigma;
		enrich_s_loc.push_back(temp_s);
		enrich_d_top.push_back(temp_top);
		enrich_d_bottom.push_back(temp_bottom);
		enrich_samp.push_back(temp_enrich);
		enrich_1sigma.push_back(temp_1sigma);
	}
	enrich_in.close();

	// get the model parameters
	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt >> temp >> CRN_switch >> temp >> end_time
					>> temp >> constant_surface_change_rate >> temp >> particle_printing_interval
					>> temp >> eroded_catch_window >> temp >> max_age
					>> temp >> n_spacings >>  temp >> insert_interval_time
					>> temp >> starting_time >> temp >> reference_frame_printing_switch;
	model_run_params_in.close();

	cout << "Flux switch: " << flux_switch << " prod_switch: " << prod_switch
		 << " Flux_us: " << flux_us << " CRN switch: " << CRN_switch << endl
	     << "dt: " << dt << " end_time: " << end_time << " surface erate: " << constant_surface_change_rate << endl
	     << "particle print interval: " << particle_printing_interval
	     << " ecatch: " << eroded_catch_window << " max age: " << max_age << " n_space: " << n_spacings << endl
	     << "insert interval: " << insert_interval_time << " start erate index: " << starting_time
	     << " ref frame switch: " << reference_frame_printing_switch << endl;

	// get the parameters for the CRN particles
	ifstream CRN_parameter_in;
	CRN_parameter_in.open(CRN_parameter_fname.c_str());
	CRN_parameter_in >> temp >> start_depth >> temp >> vert_mix_vel
				     >> temp >> horiz_mix_vel >> temp >> Omega
				     >> temp >> part_conc >> temp >> CRN_muon_param_switch
				     >> temp >> single_scaling >> temp >> C_10Be >> temp >> C_26Al
				     >> temp >> C_36Cl >> temp >> C_14C >> temp >> C_21Ne >> temp
				     >> C_3He >> temp >> M_supply_surface >> temp >> k_f10Be >> temp
				     >> deltad >> temp >> k2_f10Be >> temp >> chi_f10Be;
	CRN_parameter_in.close();

	cout << "start depth: " << start_depth << " vert mix vel: " << vert_mix_vel << endl
	     << "horiz mix vel: " << horiz_mix_vel << " Omega " << Omega << " part_conc: " << part_conc << endl
	     << "Muon switch: " << CRN_muon_param_switch << " scaling: " << single_scaling << endl
	     << "Init conc, 10Be: " << C_10Be << " 26Al: " << C_26Al << " C_36Cl: " << C_36Cl << " C_14C: " << C_14C << endl
	     << "C_21Ne: " << C_21Ne << " C_3He: " << C_3He << " M_supp_surface: " << M_supply_surface << endl
	     << "K_f10Be: " << k_f10Be << " deltad: " << deltad << endl;

	// load in the erate history
	ifstream erh_in;
	erh_in.open(erh_fname.c_str());
	end_time = load_erate_history(erh_in,erate_times,
						erate_rates);
	erh_in.close();

	// get the particle types
	Particle_info pi(particle_types_fname.c_str());
	vector<int> starting_pID = pi.get_type_index();
	int n_ptypes = starting_pID.size();

	cout << endl << endl << endl << "n particle types: " << n_ptypes << endl;

	ifstream pmfrac_in;
	pmfrac_in.open(part_mfracs_fname.c_str());
	vector<double> starting_p_mfrac(n_ptypes);
	for (int i = 0; i<n_ptypes; i++)
	{
		pmfrac_in >> starting_p_mfrac[i];
		cout << "type " << i << " frac: " << starting_p_mfrac[i] << endl;
	}
	pmfrac_in.close();

	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	// work with the outfiles
	// now open the zeta, eta, and h outfiles
	ofstream zeta_out,h_out,eta_out;
	zeta_out.open(zeta_out_fname.c_str());
	h_out.open(h_out_fname.c_str());
	eta_out.open(eta_out_fname.c_str());

	// open some file for printing
	ofstream age_cdf_out,age_pdf_out;
	age_cdf_out.open(age_cdf_out_fname.c_str());
	age_pdf_out.open(age_pdf_out_fname.c_str());

	// an integer used for printing
	const int part_p_i = int (double(particle_printing_interval/dt+0.5) );
	double next_catch = insert_interval;

	// open the datafile for the particle information
	ofstream particle_out;
	particle_out.open(particle_out_fname.c_str());
	ofstream eroded_particle_out;
	eroded_particle_out.open(eroded_pout_fname.c_str());

	// initialize the CRN parameters
	if (CRN_muon_param_switch == 1)
	{
		CRNp.set_Schaller_parameters();
	}
	CRNp.scale_F_values(vector<bool> nuclides_for_scaling);

	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname,
						  profile_in_fname);

	cout << "Line 1367 flowtube initialized!" << endl;

	// initialize the surface erosion rate
	int n_ft_nodes = ft_test.get_n_nodes();
	for (int node = 0; node < n_ft_nodes; node++)
	{
		surf_erate.push_back(constant_surface_change_rate);
	}

	// initialize a CRN_particle_list
	CRN_tParticle_bins CRN_tpb(ft_test);
	int n_bins = CRN_tpb.get_n_bins();
    vector< list<LSDCRNParticle> > eroded_catcher(n_bins+1);
    vector< list<LSDCRNParticle> > empty_eroded_catcher(n_bins+1);
    list<LSDCRNParticle>::iterator part_iter;	// list iterator

	// raise the flowtube so the downstream boundary is at elevation
	// 100
	vector<double> zeta_old = ft_test.get_zeta();
	double zeta_zero = 100.0;
	ft_test.raise_zeta_eta_ds_bound(zeta_zero);

	double rho_s = ft_test.get_rho_s();

	double insert_time_clock = 0;
	old_eta = ft_test.get_eta();
	vector<double> Delta_eta = old_eta;
	vector<double> old_bottom_depth = old_eta;
	vector<double> old_zeta = ft_test.get_zeta();
	vector<double> old_h = ft_test.get_h();
	int eta_sz = Delta_eta.size();
	int part_ID_start = 1;

	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;


	// insert initial particles
	// to insert particles throughout the dominan we need to se old botom depth as the
	// zeta elevations.
	for (int i = 0; i<eta_sz; i++)
	{
			//cout << "LINE 1411 h["<<i<<"]: " << old_h[i] << " and zeta: " << old_eta[i] << endl;
			old_bottom_depth[i] = old_zeta[i];
			Delta_eta[i] = start_depth;
	}
	part_ID_start =CRN_tpb.insert_particles(ft_test, Delta_eta, old_bottom_depth,
								 part_conc, starting_pID, starting_p_mfrac,
								 C_10Be, C_26Al, C_36Cl, C_14C, C_21Ne, C_3He);


	// print initial conditions
	//CRN_tpb.print_particle_stats(t_ime, ft_test, particle_out);
	ft_test.print_h_s(zeta_out);
	ft_test.print_h_s(eta_out);
	ft_test.print_h_s(h_out);
	ft_test.print_zeta(t_ime, zeta_out);
	ft_test.print_eta(t_ime, eta_out);
	ft_test.print_h(t_ime, h_out);

	cout << "LINE 379, n_nodes: " << ft_test.get_n_nodes() << endl;

	// now get information about erate
	n_ts = erate_times.size();
	end_time = erate_times[n_ts-1];

	curr_erate_index = 0;
	old_erate = erate_rates[curr_erate_index];
	old_erate_time = erate_times[curr_erate_index];
	new_erate = erate_rates[curr_erate_index+1];
	new_erate_time = erate_times[curr_erate_index+1];
	erate_slope = (new_erate-old_erate)/(new_erate_time-old_erate_time);

	// ramp up to the starting time
	while (t_ime+0.0001 < starting_time)
	{
		t_ime+= dt;
		tt++;

		// ramp the erate history
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
	}
	cout << "LINE 1123, ramped; time is: " << t_ime << endl;
	ds_elev = ft_test.get_zeta_ds();
	cout << "elev above ds_elev is: " << ft_test.get_zeta_ds() - ds_elev << " erate is: " << erate << endl;

	// now loop through time
	int particle_trigger = 1;
	while(t_ime < end_time)
	{
		tt++;
		t_ime += dt;		// increment the time
		insert_time_clock+=dt;	// increment the insert time clock
								// particles are not inserted every timestep
								// rather they are inserted at intervals.
								// when the insert time clock exceeds the insert
								// interval particles are inserted
								// and the insert time clock is reset
		//cout << "LINE 1441 time is: " << t_ime << endl;

		// get erate from erh
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
			cout << "updated erate, new erate, time: " << t_ime << " erate: " << new_erate << endl;
		}
		ds_elev -= erate*dt;

		//if (tt%1000 == 0)
		//{
		//	cout << "LINE 431 time is: " << t_ime << " and rate is: " << erate << " and ds_elev: " << ds_elev << endl;
		//}

		// run flux
		ft_test.flux_timestep_elev_bc(dt, flux_us, ds_elev,flux_switch, prod_switch,
							surf_erate);

		// run particle motion if particles have been inserted.
		// particle_trigger == 1 when particles have been inserted
		if (particle_trigger == 1)
		{
			eroded_bins = CRN_tpb.particle_motion(dt, ft_test,
										Omega, vert_mix_vel,
										horiz_mix_vel,CRN_switch, CRNp);

			// if the time is within the eroded catch window, catch
			// all the particles in the eroded bin
			if (t_ime > next_catch-eroded_catch_window)
			{
				//cout << "LINE 423, time is: " << t_ime << " and next_catch: " << next_catch << endl;
				for(int bn = 0; bn<=n_bins; bn++)
				{
					if (eroded_bins[bn].size()>0)
					{
						part_iter = eroded_bins[bn].begin();
						while(part_iter != eroded_bins[bn].end())
						{
							eroded_catcher[bn].push_back(*part_iter);
							part_iter++;
						}
					}
				}
			}
		}


		// if the time elapsed sinse the last insertion of particles is gereater
		// than or equal to the insertion interval, then get the amount of
		// soil-saprolite boundary lowering that has occurred and insert
		// particles into the insertion zone
		if (insert_time_clock >= insert_interval_time - dt/2)
		{
			// if this is the first time we are inserting particles,
			// set the particle trigger to 1
			if (particle_trigger == 0)
			{
				particle_trigger = 1;
			}
			new_eta = ft_test.get_eta(); 	// get the updated eta

			//cout << "Time: " << t_ime << " Delta eta: " << endl;
			// calculate delta eta
			for(int ii = 0; ii< eta_sz; ii++)
			{
				Delta_eta[ii] = old_eta[ii]-new_eta[ii];
				//cout << "i = " << ii << " Delta_eta: " << Delta_eta[ii] << endl;
			}

			// insert particles into the insertion zone

			part_ID_start =CRN_tpb.insert_particles(ft_test, Delta_eta, old_bottom_depth,
								 part_conc, starting_pID, starting_p_mfrac);

			old_eta = new_eta;				// reset old eta
			insert_time_clock = 0;			// reset the insert time clock

		}

		// print the particle data to file
		int n_depthintervals_soil = 5;
		int n_depthintervals_parent = 5;
		double bottom_depth = 2.0;
		tt = int(t_ime/dt);
		if (particle_trigger == 1 && tt%part_p_i== 0)
		{
			cout << "LINE 289, printing particles, time: " << t_ime << endl;
			ft_test.print_zeta(t_ime, zeta_out);
			ft_test.print_eta(t_ime, eta_out);
			ft_test.print_h(t_ime, h_out);
			CRN_tpb.cell_and_particle_chemistry_printing_vtk(t_ime, ft_test,
								 pi, vtk_particle_fname, vtk_cell_fname,
								 n_depthintervals_soil, n_depthintervals_parent,
								 bottom_depth, reference_frame_printing_switch);

			//CRN_tpb.print_eroded_stats(t_ime,eroded_catcher,ft_test,eroded_particle_out);
			//next_catch+=particle_print_interval;	// reset the next catch time
			//eroded_catcher = empty_eroded_catcher;
												// reset the eroded catcher
		}


		if (tt%int(dt*1000) == 0)
		{
			cout << "time is: " << t_ime << " and erate is: " << erate << endl;

		}
	}

	// calculate the MLE for both the meteoric and the in situ 10Be
	vector<double> sample_mean_age;
	vector<double> sample_mean_C10Be;
	vector<double> sample_enrich;
	vector<double> sample_frac_clay;
	vector<double> sample_frac_t1;
	vector<double> sample_frac_t2;
	vector<double> sample_frac_t3;
	vector<double> sample_frac_t4;
	vector<double> sample_frac_t5;
	CRN_tpb.calculate_sample_averages(enrich_s_loc, enrich_d_top, enrich_d_bottom,
										pi,sample_mean_age, sample_mean_C10Be, sample_enrich,
										sample_frac_clay,sample_frac_t1, sample_frac_t2,
										sample_frac_t3, sample_frac_t4,sample_frac_t5);

	int n_samples = enrich_d_top.size();

	sample_s_locs = enrich_s_loc;
	d_top = enrich_d_top;
	d_bottom = enrich_d_bottom;
	meas = enrich_samp;
	unc = enrich_1sigma;
	modelled = sample_enrich;
	c_frac = sample_frac_clay;
	pf1 = sample_frac_t1;
	pf2 = sample_frac_t2;
	pf3 = sample_frac_t3;
	pf4 = sample_frac_t4;
	pf5 = sample_frac_t5;
	C10Be = sample_mean_C10Be;

	double enrich_MLE;
	enrich_MLE = calculate_MLE(enrich_samp, sample_enrich, enrich_1sigma);
	cout << "MLE is: " << enrich_MLE << endl;


	// now print the result to a profile file
	ofstream profile_out;
	profile_out.open(profile_out_fname.c_str());
	ft_test.export_input_profile(profile_out);
	profile_out.close();
	particle_out.close();
	eroded_particle_out.close();
	zeta_out.close();
	eta_out.close();
	h_out.close();
	age_cdf_out.close();
	age_pdf_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function runs a flowtube at a steady erosion rate, and
// tests the particle inserter
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<double> part_ft_moraine_fit(string run_name, vector<double>& d_top_insitu, vector<double>& d_bottom_insitu,
							vector<double>& meas_insitu, vector<double>& modelled_insitu, double& MLE_insitu,
							vector<double>& d_top_meteor, vector<double>& d_bottom_meteor,
							vector<double>& meas_meteor, vector<double>& modelled_meteor, double& MLE_meteor)
{
	//long seed = time(NULL);               // seed for random number generator

	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=
	// set up the infiles
	string sed_trans_param_ext = ".sed_trans_param.stparam";
	string sed_trans_param_fname = run_name+sed_trans_param_ext;
	string model_run_params_ext = ".model_run.param";
	string model_run_params_fname = run_name+model_run_params_ext;
	string CRN_parameter_ext = ".CRN_trans_param.CRNparam";
	string CRN_parameter_fname = run_name+CRN_parameter_ext;
	string ft_parameter_ext =  ".ft_details.param";
	string ft_parameter_fname = run_name+ft_parameter_ext;
	string profile_in_ext = ".profile.sm";
	string profile_in_fname = run_name+profile_in_ext;
	string particle_types_ext = ".four_comp.clist";
	string particle_types_fname = run_name+particle_types_ext;
	string part_mfracs_ext = ".mfrac.mfrac";
	string part_mfracs_fname = run_name+part_mfracs_ext;
	string profile_10Be_ext = ".profile_10Be.data";
	string profile_10Be_fname = run_name+profile_10Be_ext;
	string profile_f10Be_ext = ".profile_f10Be.data";
	string profile_f10Be_fname = run_name+profile_f10Be_ext;

	// set up the outfiles
	string profile_out_ext = ".column_out.sm";
	string profile_out_fname = run_name+profile_out_ext;
	string particle_out_ext = ".p_trans_out.pout";
	string particle_out_fname = run_name+particle_out_ext;
	string eroded_pout_ext = ".ep_trans_out.pout";
	string eroded_pout_fname = run_name+eroded_pout_ext;
	string zeta_out_ext = ".zeta_trans.zdat";
	string zeta_out_fname = run_name+zeta_out_ext;
	string h_out_ext = ".h_trans.hdat";
	string h_out_fname = run_name+h_out_ext;
	string eta_out_ext = ".eta_trans.edat";
	string eta_out_fname = run_name+eta_out_ext;
	string age_cdf_out_ext = ".age_cdf.acdf";
	string age_cdf_out_fname = run_name+age_cdf_out_ext;
	string age_pdf_out_ext = ".age_pdf.apdf";
	string age_pdf_out_fname = run_name+age_pdf_out_ext;
	string vtk_cell_ext = ".Cell_data";
	string vtk_cell_fname = run_name+vtk_cell_ext;
	string vtk_particle_ext = ".Particle_data";
	string vtk_particle_fname = run_name+vtk_particle_ext;
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=


	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=
	// Initialize the variables
	double SS_flux;										// flux
	double constant_surface_change_rate;
	vector<double> surf_erate;							// erosion from the surface
														// in m/yr: this is negative for
														// erosion

	int CRN_switch;				// sets the cosmogenics key:
										// 0 == no CRN
										// 1 == all CRN
										// 2 == all, neutron only
										// 3 == 10Be full
										// default == all, neutron only

	double t_ime;
	double start_time;
	int tt;					// a counter for the time

	int flux_switch;				// see flowtube.h: determines flux law
	int prod_switch;				// see flowtube.h: determines soil production law
	double flux_us;					// see flowtube.h: flux from upslope
	double dt;						// time interval

	double start_depth;				// the starting depth in m
									// (LATER) it is set to 60 because
									// this is roughly 3x the e-folding
									// depth of the deepest muonogenic
									// production mechanism at rock density
	double vert_mix_vel;
	double horiz_mix_vel;
	double Omega;					// the activity of particles (that is the
									// proportion of particles that are moving
									// at any given time)

									// the elevation of the introduced particle
	double part_conc;				// particle concentration in particles per kg

	vector<double> old_eta;			// the elevation of the soil-saprolite boudnary from the
									// last timestep
	vector<double> new_eta;			// the elevation of the soil-saprolite boundary from this
									// timestep
	vector< list<LSDCRNParticle> > eroded_bins;
	vector< list<LSDCRNParticle> > temp_part_bins;
	vector< list<LSDCRNParticle> > particle_bins;
									// a vector of lists of particles eroded from each hillslope

	// initial in situ concentrations in atoms per gram
	double C_10Be,C_26Al,C_36Cl,C_14C,C_21Ne,C_3He;

	// set the parameters for production
	// 0 == granger
	// 1 == schaller
	int CRN_muon_param_switch;

	// the parameters for the in situ cosmogenics
	LSDCRNParameters CRNp;

	// note: scaling determined using cosmocalc:
	// Schaller reports pinedale at 42 53 26 N
	// this converts to an inclination of 61.7 degrees
	// Schaller reports elevation of 2298 masl
	// this converts to an atmospheric depth of 781 g/cm2
	// this converts to a dunai scaling of 5.99
	// multiply this by schaller's snow shielding (for nucleonic production
	// of 0.925 one gets a scaling of 5.54
	double single_scaling;

	// paramters for meteoric 10Be
	double M_supply_surface;	// in atoms/(cm^2*yr)
	double k_f10Be;				// in cm^2/g this is
									// an efolding depth of 20 g/cm^2
	double k2_f10Be;
	double chi_f10Be;			// fraction of supply that goes into the
								// shallow meteoric supply
	double deltad;				// in m (this gets converted
								// to cm in LSDParticle.cpp)

	// parameters for dealing with time
	double end_time;
	double insert_interval;
	double particle_printing_interval;
	double zhe_printing_interval;
	double eroded_catch_window;
	double age_printing_interval;
	double max_age;
	int n_spacings;

	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	//=-=-=-=-=-=-=-=
	// load files
	//=-=-=-=-=-=-=-=
	// load the profile of in situ 10Be
	ifstream InSitu_10Be_in;
	InSitu_10Be_in.open(profile_10Be_fname.c_str());
	int n_10Be_samples;
	InSitu_10Be_in >> n_10Be_samples;
	vector<double> InSitu_d_top;
	vector<double> InSitu_d_bottom;
	vector<double> InSitu_C10Be;
	vector<double> InSitu_1sigma;
	double temp_top,temp_bottom,temp_C10Be,temp_1sigma;
	for (int i = 0; i< n_10Be_samples; i++)
	{
		InSitu_10Be_in >> temp_top >> temp_bottom >> temp_C10Be >> temp_1sigma;
		InSitu_d_top.push_back(temp_top);
		InSitu_d_bottom.push_back(temp_bottom);
		InSitu_C10Be.push_back(temp_C10Be);
		InSitu_1sigma.push_back(temp_1sigma);
	}
	InSitu_10Be_in.close();

	// load the profile of meteoric 10Be
	ifstream Meteor_10Be_in;
	Meteor_10Be_in.open(profile_f10Be_fname.c_str());
	int n_f10Be_samples;
	Meteor_10Be_in >> n_f10Be_samples;
	cout << "LINE 1288 number_of meteoric samples; " << n_f10Be_samples << endl;
	vector<double> Meteor_d_top;
	vector<double> Meteor_d_bottom;
	vector<double> Meteor_C10Be;
	vector<double> Meteor_1sigma;
	double temp_fC10Be;
	for (int i = 0; i< n_f10Be_samples; i++)
	{
		Meteor_10Be_in >> temp_top >> temp_bottom >> temp_fC10Be >> temp_1sigma;
		Meteor_d_top.push_back(temp_top);
		Meteor_d_bottom.push_back(temp_bottom);
		Meteor_C10Be.push_back(temp_fC10Be);
		Meteor_1sigma.push_back(temp_1sigma);
		//cout << "d top: " << temp_top << " d bot: " << temp_bottom
		//       << " conc: " << temp_fC10Be << " temp unc: " << temp_1sigma << endl;
	}
	Meteor_10Be_in.close();

	// get the model parameters
	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt >> temp >> CRN_switch >> temp >> end_time
					>> temp >> constant_surface_change_rate >> temp >> particle_printing_interval
					>> temp >> eroded_catch_window >> temp >> max_age
					>> temp >> n_spacings;
	model_run_params_in.close();

	//cout << "Flux switch: " << flux_switch << " prod_switch: " << prod_switch
	//	 << " Flux_us: " << flux_us << " CRN switch: " << CRN_switch << endl
	//     << "dt: " << dt << " end_time: " << end_time << " surface erate: " << constant_surface_change_rate << endl
	//     << "particle print interval: " << particle_printing_interval
	//     << " ecatch: " << eroded_catch_window << " max age: " << max_age << " n_space: " << n_spacings << endl;

	// get the parameters for the CRN particles
	ifstream CRN_parameter_in;
	CRN_parameter_in.open(CRN_parameter_fname.c_str());
	CRN_parameter_in >> temp >> start_depth >> temp >> vert_mix_vel
				     >> temp >> horiz_mix_vel >> temp >> Omega
				     >> temp >> part_conc >> temp >> CRN_muon_param_switch
				     >> temp >> single_scaling >> temp >> C_10Be >> temp >> C_26Al
				     >> temp >> C_36Cl >> temp >> C_14C >> temp >> C_21Ne >> temp
				     >> C_3He >> temp >> M_supply_surface >> temp >> k_f10Be >> temp
				     >> deltad >> temp >> k2_f10Be >> temp >> chi_f10Be;
	CRN_parameter_in.close();

	//cout << "start depth: " << start_depth << " vert mix vel: " << vert_mix_vel << endl
	//     << "horiz mix vel: " << horiz_mix_vel << " Omega " << Omega << " part_conc: " << part_conc << endl
	//     << "Muon switch: " << CRN_muon_param_switch << " scaling: " << single_scaling << endl
	//     << "Init conc, 10Be: " << C_10Be << " 26Al: " << C_26Al << " C_36Cl: " << C_36Cl << " C_14C: " << C_14C << endl
	//     << "C_21Ne: " << C_21Ne << " C_3He: " << C_3He << " M_supp_surface: " << M_supply_surface << endl
	//     << "K_f10Be: " << k_f10Be << " deltad: " << deltad << endl;

	// get the particle types
	Particle_info pi(particle_types_fname.c_str());
	vector<int> starting_pID = pi.get_type_index();
	int n_ptypes = starting_pID.size();

	ifstream pmfrac_in;
	pmfrac_in.open(part_mfracs_fname.c_str());
	vector<double> starting_p_mfrac(n_ptypes);
	for (int i = 0; i<n_ptypes; i++)
	{
		pmfrac_in >> starting_p_mfrac[i];
	}
	pmfrac_in.close();

	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	// work with the outfiles
	// now open the zeta, eta, and h outfiles
	ofstream zeta_out,h_out,eta_out;
	zeta_out.open(zeta_out_fname.c_str());
	h_out.open(h_out_fname.c_str());
	eta_out.open(eta_out_fname.c_str());

	// open some file for printing
	ofstream age_cdf_out,age_pdf_out;
	age_cdf_out.open(age_cdf_out_fname.c_str());
	age_pdf_out.open(age_pdf_out_fname.c_str());

	// an integer used for printing
	const int part_p_i = int (double(particle_printing_interval/dt+0.5) );
	double next_catch = insert_interval;

	// open the datafile for the particle information
	ofstream particle_out;
	particle_out.open(particle_out_fname.c_str());
	ofstream eroded_particle_out;
	eroded_particle_out.open(eroded_pout_fname.c_str());


	// initialize the CRN parameters
	if (CRN_muon_param_switch == 1)
	{
		CRNp.set_Schaller_parameters();
	}
	CRNp.scale_F_values(vector<bool> nuclides_for_scaling);

	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname,
						  profile_in_fname);

	// initialize the surface erosion rate
	int n_ft_nodes = ft_test.get_n_nodes();
	for (int node = 0; node < n_ft_nodes; node++)
	{
		surf_erate.push_back(constant_surface_change_rate);
	}

	// initialize a CRN_particle_list
	CRN_tParticle_bins CRN_tpb(ft_test);
	int n_bins = CRN_tpb.get_n_bins();
    vector< list<LSDCRNParticle> > eroded_catcher(n_bins+1);
    vector< list<LSDCRNParticle> > empty_eroded_catcher(n_bins+1);
    list<LSDCRNParticle>::iterator part_iter;	// list iterator

	// raise the flowtube so the downstream boundary is at elevation
	// 100
	vector<double> zeta_old = ft_test.get_zeta();
	double zeta_zero = 100.0;
	ft_test.raise_zeta_eta_ds_bound(zeta_zero);

	double rho_s = ft_test.get_rho_s();

	double insert_time_clock = 0;
	old_eta = ft_test.get_eta();
	vector<double> Delta_eta = old_eta;
	vector<double> old_bottom_depth = old_eta;
	vector<double> old_zeta = ft_test.get_zeta();
	vector<double> old_h = ft_test.get_h();
	int eta_sz = Delta_eta.size();
	int part_ID_start = 1;

	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;

	// insert initial particles
	// to insert particles throughout the dominan we need to se old botom depth as the
	// zeta elevations.
	for (int i = 0; i<eta_sz; i++)
	{
			cout << "LINE 371 h["<<i<<"]: " << old_h[i] << " and zeta: " << old_eta[i] << endl;
			old_bottom_depth[i] = old_zeta[i];
			Delta_eta[i] = start_depth;
	}
	part_ID_start =CRN_tpb.insert_particles(ft_test, Delta_eta, old_bottom_depth,
								 part_conc, starting_pID, starting_p_mfrac,
								 C_10Be, C_26Al, C_36Cl, C_14C, C_21Ne, C_3He);

	// print initial conditions
	//CRN_tpb.print_particle_stats(t_ime, ft_test, particle_out);
	ft_test.print_h_s(zeta_out);
	ft_test.print_h_s(eta_out);
	ft_test.print_h_s(h_out);
	ft_test.print_zeta(t_ime, zeta_out);
	ft_test.print_eta(t_ime, eta_out);
	ft_test.print_h(t_ime, h_out);

	cout << "LINE 379, n_nodes: " << ft_test.get_n_nodes() << endl;

	// now loop through time
	int particle_trigger = 1;
	while(t_ime < end_time)
	{
		tt++;
		t_ime += dt;		// increment the time
		insert_time_clock+=dt;	// increment the insert time clock
								// particles are not inserted every timestep
								// rather they are inserted at intervals.
								// when the insert time clock exceeds the insert
								// interval particles are inserted
								// and the insert time clock is reset

		//ds_elev -= dt*SS_erate;

		// run flux
		//cout << "running flux ";
		ft_test.flux_timestep_flux_bc(dt, flux_us, SS_flux,flux_switch, prod_switch,
							surf_erate);
		//cout << "...ran flux" << endl;

		//cout << "running motion";
		// run particle motion if particles have been inserted.
		// particle_trigger == 1 when particles have been inserted
		if (particle_trigger == 1)
		{
			eroded_bins = CRN_tpb.particle_motion(dt, ft_test,
										Omega, vert_mix_vel,
										horiz_mix_vel, CRN_switch, CRNp);

			if (chi_f10Be == 0)
			{
				CRN_tpb.update_fallout_10Be_bins(dt, M_supply_surface,
						rho_s, k_f10Be, deltad, CRNp);
			}
			else
			{
				CRN_tpb.update_fallout_10Be_bins(dt, M_supply_surface,
						rho_s, k_f10Be, k2_f10Be, chi_f10Be, deltad, CRNp);
			}

			// if the time is within the eroded catch window, catch
			// all the particles in the eroded bin
			if (t_ime > next_catch-eroded_catch_window)
			{
				//cout << "LINE 423, time is: " << t_ime << " and next_catch: " << next_catch << endl;
				for(int bn = 0; bn<=n_bins; bn++)
				{
					if (eroded_bins[bn].size()>0)
					{
						part_iter = eroded_bins[bn].begin();
						while(part_iter != eroded_bins[bn].end())
						{
							eroded_catcher[bn].push_back(*part_iter);
							part_iter++;
						}
					}
				}
			}
		}
		//cout << "...ran motion" << endl;


		//cout << "running insertion";
		// if the time elapsed sinse the last insertion of particles is gereater
		// than or equal to the insertion interval, then get the amount of
		// soil-saprolite boundary lowering that has occurred and insert
		// particles into the insertion zone
		if (insert_time_clock >= insert_interval - dt/2)
		{
			// if this is the first time we are inserting particles,
			// set the particle trigger to 1
			if (particle_trigger == 0)
			{
				particle_trigger = 1;
			}
			new_eta = ft_test.get_eta(); 	// get the updated eta

			//cout << "Time: " << t_ime << " Delta eta: " << endl;
			// calculate delta eta
			for(int ii = 0; ii< eta_sz; ii++)
			{
				Delta_eta[ii] = old_eta[ii]-new_eta[ii];
				//cout << "i = " << ii << " Delta_eta: " << Delta_eta[ii] << endl;
			}

			// insert particles into the insertion zone
			part_ID_start =CRN_tpb.insert_particles(ft_test, Delta_eta, old_bottom_depth,
								 part_conc, starting_pID, starting_p_mfrac);

			old_eta = new_eta;				// reset old eta
			insert_time_clock = 0;			// reset the insert time clock

		}
		//cout << "...ran insertion" << endl;

		//cout << "running printing" << endl;
		// print the particle data to file
		int n_depthintervals_soil = 2;
		int n_depthintervals_parent = 3;
		double bottom_depth = 2.0;
		if (particle_trigger == 1 && tt%part_p_i== 0)
		{
			cout << "LINE 289, printing particles, time: " << t_ime << endl;
			ft_test.print_zeta(t_ime, zeta_out);
			ft_test.print_eta(t_ime, eta_out);
			ft_test.print_h(t_ime, h_out);
			int ref_frame_switch = 0;
			CRN_tpb.cell_and_particle_chemistry_printing_vtk(t_ime, ft_test,
											 pi, vtk_particle_fname, vtk_cell_fname,
											 n_depthintervals_soil, n_depthintervals_parent,
											 bottom_depth, ref_frame_switch);

		}
		//cout << "...ran printing" << endl;

		if (tt%2500 == 0)
		{
			cout << "time is: " << t_ime << endl;

		}
	}

	// now print the result to a profile file
	ofstream profile_out;
	profile_out.open(profile_out_fname.c_str());
	ft_test.export_input_profile(profile_out);

	// close files
	profile_out.close();
	particle_out.close();
	eroded_particle_out.close();
	zeta_out.close();
	eta_out.close();
	h_out.close();
	age_cdf_out.close();
	age_pdf_out.close();

	// calculate the MLE for both the meteoric and the in situ 10Be
	int bn = 0;
	vector<double> sample_mean_age;
	vector<double> sample_mean_C10Be;
	vector<double> sample_mean_Cf10Be;
	CRN_tpb.calculate_sample_averages(bn, InSitu_d_top, InSitu_d_bottom,
										sample_mean_age, sample_mean_C10Be,
										sample_mean_Cf10Be);

	int n_samples = InSitu_d_top.size();
	//for (int samp = 0; samp < n_samples; samp++)
	//{
	//	cout << "sample " << samp << " meas 10Be: " << InSitu_C10Be[samp]
	//		 << " modelled: " << sample_mean_C10Be[samp] << endl;
	//}

	double InSitu_MLE;
	//double sigma = 2e8;
	vector<double> meas_sm;
	vector<double> mod_sm;
	vector<double> sig_sm;
	for (int samp = 0; samp < n_samples; samp++)
	{
		meas_sm.push_back(InSitu_C10Be[samp]/1e5);
		mod_sm.push_back(sample_mean_C10Be[samp]/1e5);
		sig_sm.push_back(InSitu_1sigma[samp]/1e5);
	}
	InSitu_MLE = calculate_MLE(meas_sm, mod_sm, sig_sm);
	//cout << "MLE is: " << InSitu_MLE << endl;

	// now do it for the fallout samples
	CRN_tpb.calculate_sample_averages(bn, Meteor_d_top, Meteor_d_bottom,
										sample_mean_age, sample_mean_C10Be,
										sample_mean_Cf10Be);

	//
	//for (int samp = 0; samp < n_samples; samp++)
	//{
	//	cout << "sample " << samp << " meas meteoric 10Be: " << Meteor_C10Be[samp]
	//		 << " modelled: " << sample_mean_Cf10Be[samp] << endl;
	//}

	double Meteoric_MLE;
	n_samples = Meteor_d_top.size();
	//double sigma = 2e8;
	vector<double> meas_sm_m;
	vector<double> mod_sm_m;
	vector<double> sig_sm_m;
	for (int samp = 0; samp < n_samples; samp++)
	{
		meas_sm_m.push_back(Meteor_C10Be[samp]/1e6);
		mod_sm_m.push_back(sample_mean_Cf10Be[samp]/1e6);
		sig_sm_m.push_back(Meteor_1sigma[samp]/1e6);
	}
	Meteoric_MLE = calculate_MLE(meas_sm_m, mod_sm_m, sig_sm_m);
	//cout << "meteoric MLE is: " << Meteoric_MLE << endl;

	d_top_insitu = InSitu_d_top;

	//cout << "size d_top_insitu: " << d_top_insitu.size();
	//for (int s = 0;  s < d_top_insitu.size(); s++)
	//{
	//	cout << "LINE 1662: " << d_top_insitu[s] << " " << InSitu_d_top[s] << endl;
	//}

	d_bottom_insitu = InSitu_d_bottom;
	meas_insitu = meas_sm;
	modelled_insitu = mod_sm;
	MLE_insitu = InSitu_MLE;
	d_top_meteor = Meteor_d_top;
	d_bottom_meteor = Meteor_d_bottom;
	meas_meteor = meas_sm_m;
	modelled_meteor = mod_sm_m;
	MLE_meteor = Meteoric_MLE;


	vector<double> MLE_vec(2);
	MLE_vec[0] = InSitu_MLE;
	MLE_vec[1] = Meteoric_MLE;




	return MLE_vec;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Automation and file manipulation
//
void modify_parameter_files(string run_name, double start_depth, double vert_mix_vel, double part_conc,
					double single_scaling, double M_supply_surface,
					double k_f10Be, double deltad, double k2_f10Be, double chi_f10Be,
					double W_0, double gamma)
{


	ifstream master_mrun;
	master_mrun.open("./master_params/CRN_trans_param.CRNparam");

	char line[500];
	vector<string> file_lines;

	string temp;
	while(master_mrun.getline(line,500))
	{
		temp = line;
		file_lines.push_back(temp);
	}
	master_mrun.close();

	//int n_files = file_lines.size();
	//for(int i = 0; i<n_files; i++)
	//{
	//	cout << file_lines[i] << endl;
	//}

	ofstream CRN_out;
	string CRN_fname = ".CRN_trans_param.CRNparam";
	CRN_fname = run_name+CRN_fname;
	//cout << "CRN_fname: " << CRN_fname << endl;
	CRN_out.open(CRN_fname.c_str());
	//cout << "opened file" << endl;

	CRN_out << "start_depth: " << start_depth << endl;
	CRN_out << "vert_mix_vel: " << vert_mix_vel << endl;
	CRN_out << file_lines[2] << endl << file_lines[3] << endl;
	CRN_out << "part_conc: " << part_conc << endl;
	CRN_out << file_lines[5] << endl;
	CRN_out << "single_scaling: " << single_scaling << endl;
	CRN_out << file_lines[7] << endl << file_lines[8] << endl << file_lines[9] << endl
			<< file_lines[10] << endl << file_lines[11] << endl << file_lines[12] << endl;
	CRN_out << "M_supply_surface: " << M_supply_surface << endl;
	CRN_out << "k_f10Be: " << k_f10Be << endl;
	CRN_out << "deltad: " << deltad << endl;
	CRN_out << "k2_f10Be: " << k2_f10Be << endl;
	CRN_out << "chi_f10Be: " << chi_f10Be << endl;

	CRN_out.close();

	ifstream st_in;
	st_in.open("./master_params/sed_trans_param.stparam");
	vector<string> st_lines;
	//cout << "opened file" << endl;
	while(st_in.getline(line,500))
	{
		temp = line;
		st_lines.push_back(temp);
	}
	st_in.close();

	//n_files = st_lines.size();
	//for(int i = 0; i<n_files; i++)
	//{
	//	cout << st_lines[i] << endl;
	//}

	ofstream ST_out;
	string ST_fname = ".sed_trans_param.stparam";
	ST_fname = run_name+ST_fname;
	ST_out.open(ST_fname.c_str());

	ST_out << st_lines[0] << endl << st_lines[1] << endl << st_lines[2]
	       << endl << st_lines[3] << endl;
	ST_out << "W_0: " << W_0 << endl;
	ST_out << "gamma: " << gamma << endl;

	ST_out.close();
}

#endif
