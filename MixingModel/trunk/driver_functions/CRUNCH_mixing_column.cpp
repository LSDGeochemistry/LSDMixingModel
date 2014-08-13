#include <iostream>
#include <fstream>
#include <vector>
#include "../tParticle.hpp"
#include "../CRN_tParticle_bins.hpp"
#include "../mathutil.hpp"
#include "../CRUNCH_engine.hpp"
#include "../flowtube.hpp"
#include "../FT_util.hpp"
#include "../VolumeParticleInfo.hpp"
#include "../CRUNCH_bins.hpp"
using namespace std;

int main ()
{
	//long seed = time(NULL);               // seed for random number generator
	string run_name = "c:/code/devel_projects/MixingModel/Runs/Run1/run1";
	//string run_name = "M:/papers/mixing_model_2014/source/runs/run1/run1";

	//string crunch_pname = "M:/papers/mixing_model_2014/source/CRUNCH_binary/";
	//string run_pname = "M:/papers/mixing_model_2014/source/runs/run1/"; 
	string crunch_pname = "c:/code/devel_projects/MixingModel/CRUNCH_binary/";
	string run_pname = "c:/code/devel_projects/MixingModel/Runs/Run1/"; 
	
	string vtk_particle_fname = run_pname+"/basic_particles";
	string vtk_cell_fname = run_pname+"/CRUNCH_cells";


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
	string VolumeParticleInfo_ext = ".VolumeParticleData.in";
	string VolumeParticleInfo_fname = run_name+VolumeParticleInfo_ext;
	//string profile_10Be_ext = ".profile_10Be.data";
	//string profile_10Be_fname = run_name+profile_10Be_ext;
	//string profile_f10Be_ext = ".profile_f10Be.data";
	//string profile_f10Be_fname = run_name+profile_f10Be_ext;

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
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=


  ofstream zeta_out,h_out,eta_out;
	zeta_out.open(zeta_out_fname.c_str());
	h_out.open(h_out_fname.c_str());
	eta_out.open(eta_out_fname.c_str());


	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=
	// Initialize the variables
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	
	// the number of partitions in the mixed and unmixed zone
	int n_PDZ_intervals;
	int n_CAZ_intervals;	
	
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
	vector< list<CRN_tParticle> > eroded_bins;
	vector< list<CRN_tParticle> > temp_part_bins;
	vector< list<CRN_tParticle> > particle_bins;
									// a vector of lists of particles eroded from each hillslope

	// load a volume particle info object
	cout << "LINE 115 loading particle info" << endl;
	cout << "line 115, particle info name: " << VolumeParticleInfo_fname << endl;
	VolumeParticleInfo vpi(VolumeParticleInfo_fname.c_str());
	cout << "LINE 115 loaded particle info" << endl;
	cout << "LINE 117, density type 0: " << vpi.get_type_density(0) << endl;
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	



  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  //
  // SET UP COSMO
  // 
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	// initial in situ concentrations in atoms per gram
	double C_10Be,C_26Al,C_36Cl,C_14C,C_21Ne,C_3He;

	// set the parameters for production
	// 0 == granger
	// 1 == schaller
	int CRN_muon_param_switch;

	// the parameters for the in situ cosmogenics
	CRN_parameters CRNp;

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
								// to cm in tParticle.cpp)
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  
  
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // 
  // IMPORT PARAMETERS
  //
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	// parameters for dealing with time
	double end_time;
	double insert_interval;	
	double weathering_time_interval;
	double particle_printing_interval;
	double zhe_printing_interval;
	double eroded_catch_window;
	double age_printing_interval;
	double max_age;
	int n_spacings;

	// get the model parameters
	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt >> temp >> CRN_switch >> temp >> end_time
					>> temp >> constant_surface_change_rate >> temp >> particle_printing_interval
					>> temp >> eroded_catch_window >> temp >> max_age
					>> temp >> n_spacings >> temp >> insert_interval 
          >> temp >> weathering_time_interval;
	model_run_params_in.close();
	cout << "LINE 209, got model_parameters" << endl;

	//cout << "Flux switch: " << flux_switch << " prod_switch: " << prod_switch
	//	 << " Flux_us: " << flux_us << " CRN switch: " << CRN_switch << endl
	//     << "dt: " << dt << " end_time: " << end_time << " surface erate: " << constant_surface_change_rate << endl
	//     << "particle print interval: " << particle_printing_interval
	//     << " catch: " << eroded_catch_window << " max age: " << max_age << " n_space: " << n_spacings 
  //     << " particle insert interval: " << insert_interval 
  //     << " weathering time interval: " << weathering_time_interval << endl;

	// get the parameters for the CRN particles
	ifstream CRN_parameter_in;
	CRN_parameter_in.open(CRN_parameter_fname.c_str());
	CRN_parameter_in >> temp >> start_depth >> temp >> vert_mix_vel
				     >> temp >> horiz_mix_vel >> temp >> Omega
				     >> temp >> part_conc >> temp >> CRN_muon_param_switch
				     >> temp >> single_scaling >> temp >> C_10Be >> temp >> C_26Al
				     >> temp >> C_36Cl >> temp >> C_14C >> temp >> C_21Ne >> temp
				     >> C_3He >> temp >> M_supply_surface >> temp >> k_f10Be >> temp
				     >> deltad >> temp >> k2_f10Be >> temp >> chi_f10Be
             >> temp >> n_PDZ_intervals >> temp >> n_CAZ_intervals;
	CRN_parameter_in.close();
	cout << "LINE 228, got CRN_parameters" << endl;

	//cout << "start depth: " << start_depth << " vert mix vel: " << vert_mix_vel << endl
	//     << "horiz mix vel: " << horiz_mix_vel << " Omega " << Omega << " part_conc: " << part_conc << endl
	//     << "Muon switch: " << CRN_muon_param_switch << " scaling: " << single_scaling << endl
	//     << "Init conc, 10Be: " << C_10Be << " 26Al: " << C_26Al << " C_36Cl: " << C_36Cl << " C_14C: " << C_14C << endl
	//     << "C_21Ne: " << C_21Ne << " C_3He: " << C_3He << " M_supp_surface: " << M_supply_surface << endl
	//     << "K_f10Be: " << k_f10Be << " deltad: " << deltad << endl;

	
	// an integer used for printing
	const int part_p_i = int (double(particle_printing_interval/dt+0.5) );
	cout << "printing interval is: " <<  particle_printing_interval 
       << " and part pi is: " << part_p_i << endl;
	double next_catch = insert_interval;	
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if(	weathering_time_interval > particle_printing_interval)
  {
    cout << "Fatal error, the weathering time interval must be <= particle_print interval!" << endl;
    exit(0);
  }




  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // 
  // CRUNCH INITIATION
  //
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	int tot_intervals = n_PDZ_intervals+n_CAZ_intervals;

  CRUNCH_bins Geochem_bins(n_PDZ_intervals , n_CAZ_intervals, start_depth, vpi);
  //cout << "Now getting "
  
  // create the crunch engine
	string master_fname = "master_crunch.in";
	CRUNCH_engine Ceng(crunch_pname, run_pname, master_fname);  
	
	string value_name = "spatial_profile";
	Ceng.modify_CRUNCH_value(value_name, double(weathering_time_interval));
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  
  
	// initialize the CRN parameters
	if (CRN_muon_param_switch == 1)
	{
		CRNp.set_Schaller_parameters();
	}
	CRNp.scale_F_values(single_scaling);
	cout << "scaled to schaller" << endl;

	// initialize a flowtube
	flowtube ft_test = flowtube_initializer(sed_trans_param_fname,
						  ft_parameter_fname,
						  profile_in_fname);
	cout << "LINE 292 initialized the flowtube" << endl;

	// initialize the surface erosion rate
	int n_ft_nodes = ft_test.get_n_nodes();
	for (int node = 0; node < n_ft_nodes; node++)
	{
		surf_erate.push_back(constant_surface_change_rate);
	}
	cout << "LINE 300 got surface erate" << endl;

	// initialize a CRN_particle_list
	CRN_tParticle_bins CRN_tpb(ft_test);
	int n_bins = CRN_tpb.get_n_bins();
  vector< list<CRN_tParticle> > eroded_catcher(n_bins+1);
  vector< list<CRN_tParticle> > empty_eroded_catcher(n_bins+1);
  list<CRN_tParticle>::iterator part_iter;	// list iterator

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
  double weathering_time_clock = 0;
  
	// insert initial particles
	// to insert particles throughout the domain we need to set old botom depth as the
	// zeta elevations.
	for (int i = 0; i<eta_sz; i++)
	{
			old_bottom_depth[i] = old_zeta[i];
			Delta_eta[i] = start_depth;
			//cout << "LINE 350 h["<<i<<"]: " << old_h[i] << " and zeta: " << old_zeta[i] << endl
			//     << " and start depth: " << start_depth << " and Delta_eta: " << Delta_eta[i] << endl;
	}

  part_ID_start = CRN_tpb.insert_particles_volumetric(ft_test, Delta_eta, old_bottom_depth,
										C_10Be, C_26Al, C_36Cl, C_14C, C_21Ne, C_3He,
										vpi);

	//cout << "LINE 379, n_nodes: " << ft_test.get_n_nodes() << endl;

	
	// now loop through time
	int particle_trigger = 1;
	while(t_ime < end_time)
	{
		tt++;
		t_ime += dt;		// increment the time
		
		cout << "Time is: " << t_ime << endl; 
		
		insert_time_clock+=dt;	// increment the insert time clock
								// particles are not inserted every timestep
								// rather they are inserted at intervals.
								// when the insert time clock exceeds the insert
								// interval particles are inserted
								// and the insert time clock is reset
    weathering_time_clock+=dt;  // similar to the insert clock, but this determines
                // weathering. Paricles move about without weathering, then
                // once a weathering interval has elapsed all particles are weatherd
                // for that interval and the mass updated. 

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
				cout << "i = " << ii << " Delta_eta: " << Delta_eta[ii] << endl;
			}
			
			// insert the particles
      part_ID_start = CRN_tpb.insert_particles_volumetric(ft_test, Delta_eta, old_bottom_depth,
										C_10Be, C_26Al, C_36Cl, C_14C, C_21Ne, C_3He,
										vpi);
      					 
			old_eta = new_eta;				// reset old eta
			insert_time_clock = 0;			// reset the insert time clock
		}               // !end particle insertion
		//cout << "...ran insertion" << endl;

    // now checking weathering clock
    if (weathering_time_clock >= weathering_time_interval - dt/2)
    {
      cout << "Time is: " << t_ime << " and weathering now" << endl;
    
      // run a crunch timestep
      Geochem_bins.run_CRUNCH_timestep(CRN_tpb, ft_test, Ceng);
      
      weathering_time_clock = 0;
    }

		//cout << "running printing" << endl;
		// print the particle data to file
		if (particle_trigger == 1 && tt%part_p_i== 0)
		{
			cout << "LINE 440, printing particles, time: " << t_ime << endl;
			ft_test.print_zeta(t_ime, zeta_out);
			ft_test.print_eta(t_ime, eta_out);
			ft_test.print_h(t_ime, h_out);
			int ref_frame_switch = 0;
			
			// print basic particle information
			CRN_tpb.vtk_print_basic_volume_particles(t_ime, vtk_particle_fname, 
                                               ref_frame_switch);

      // print cell data
      Geochem_bins.vtk_cell_bundler(t_ime, ref_frame_switch, 
                                    vtk_cell_fname, Ceng);

		}
		//cout << "...ran printing" << endl;

		if (tt%2500 == 0)
		{
			cout << "time is: " << t_ime << endl;

		}
  }
   
  eta_out.close();
  zeta_out.close();
  h_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

