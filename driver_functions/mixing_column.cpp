//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// mixing_column
//
// This initiates a mixing column. It can mix vertically or downslope depending on the
// parameter files
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
// either version 3 of the License, or (at your option) any later version.
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
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include "../tParticle.hpp"
#include "../CRN_tParticle_bins.hpp"
//#include "../CRUNCH_engine.hpp"
#include "../flowtube.hpp"
#include "../FT_util.hpp"
#include "../VolumeParticleInfo.hpp"

#include "../TNT/tnt.h"
#include "../LSDStatsTools.hpp"
using namespace std;

int main(int argc, char *argv[])
{

  if (argc != 2)
  {
    cout << "ERROR: you need to enter a path to the input files" << endl;
    cout << "Runs are generally stored in the folder ../Runs/run1/" << endl
         << "or ../Runs/run2/ and so on" << endl;
    cout << "Also, do you have all the driver files? You need: " << endl
         << "a sedimenent transport parameter file" << endl
         << "a model run paramters file (mostly for printing, end time, etc" << endl
         << "a particle parameter file" << endl
         << "a flowtube details file (for spacing of nodes" << endl
         << "a profile file (for measured soil thickness and elevation)" << endl
         << "a volume particle info file for mineral and size type information" << endl
         << "and a master_crunch file for default CrunchFlow parameters" << endl;
    exit(0);
  }

  string run_pname = argv[1];
  string lchar = run_pname.substr(run_pname.length()-2,1);
  string slash = "/";
  cout << "lchar is " << lchar << " and slash is " << slash << endl;

  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl;
    run_pname = run_pname+slash;
  }
  cout << "The pathname is: " << run_pname << endl;

  // This sets the path to the crunch files.
  //string crunch_pname = "M:/papers/mixing_model_2014/source/CRUNCH_binary/";
  //string run_pname = "M:/papers/mixing_model_2014/source/runs/run1/";
  string crunch_pname = "C:/Workspace/github/CRUNCH_binary/";
  //string run_pname = "../Runs/Run2/";

  // This sets up the paths for vtk files, which are used for visualisation
  string vtk_particle_fname = run_pname+"/basic_particles";
  string vtk_cell_fname = run_pname+"/CRUNCH_cells";
  string vtk_fname = run_pname+"/CRN";
  string path_to_data = "/";

  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=
  // set up the infiles
  // SMM July 2018: This is a bit of a stupid way of doing it: probably better
  // for these paths to be in an object, but no time to fix that now.
  string sed_trans_param_ext = "sed_trans_param.stparam";
  string sed_trans_param_fname = run_pname+sed_trans_param_ext;
  string model_run_params_ext = "model_run.param";
  string model_run_params_fname = run_pname+model_run_params_ext;
  string CRN_parameter_ext = "CRN_trans_param.CRNparam";
  string CRN_parameter_fname = run_pname+CRN_parameter_ext;
  string ft_parameter_ext =  "ft_details.param";
  string ft_parameter_fname = run_pname+ft_parameter_ext;
  string profile_in_ext = "profile.sm";
  string profile_in_fname = run_pname+profile_in_ext;
  string VolumeParticleInfo_ext = "VolumeParticleData.in";
  string VolumeParticleInfo_fname = run_pname+VolumeParticleInfo_ext;
  string erate_file_ext = "erate_hist.param";
  string erate_file_fname = run_pname+erate_file_ext;
  string base_level_file_ext = "base_level_hist.param";
  string base_level_hist_file_fname = run_pname+base_level_file_ext;

	// test to see if the first file exists
	ifstream test_in(model_run_params_fname.c_str());
	if (not test_in)
	{
    cout << "Can't find the model run params. Did you make sure to direct" << endl
         << "the code to the correct directory and do you have all the" << endl
         << "driver files? You need: " << endl
         << "a sedimenent transport parameter file" << endl
         << "a model run paramters file (mostly for printing, end time, etc" << endl
         << "a particle parameter file" << endl
         << "a flowtube details file (for spacing of nodes" << endl
         << "a profile file (for measured soil thickness and elevation)" << endl
         << "a volume particle info file for mineral and size type information" << endl
         << "and a master_crunch file for default CrunchFlow parameters" << endl;
         exit(0);
  }

	//string profile_10Be_ext = ".profile_10Be.data";
	//string profile_10Be_fname = run_name+profile_10Be_ext;
	//string profile_f10Be_ext = ".profile_f10Be.data";
	//string profile_f10Be_fname = run_name+profile_f10Be_ext;

	// set up the outfiles
	string profile_out_ext = "column_out.sm";
	string profile_out_fname = run_pname+profile_out_ext;
	string particle_out_ext = "p_trans_out.pout";
	string particle_out_fname = run_pname+particle_out_ext;
	string eroded_pout_ext = "ep_trans_out.pout";
	string eroded_pout_fname = run_pname+eroded_pout_ext;
	string surface_eroded_pout_ext = "surface_ep_trans_out.pout";
	string surface_eroded_pout_fname = run_pname+surface_eroded_pout_ext;
	string zeta_out_ext = "zeta_trans.zdat";
	string zeta_out_fname = run_pname+zeta_out_ext;
	string h_out_ext = "h_trans.hdat";
	string h_out_fname = run_pname+h_out_ext;
	string eta_out_ext = "eta_trans.edat";
	string eta_out_fname = run_pname+eta_out_ext;
	string hillslope_out_ext = "hillslope_out.pout";
	string hillslope_out_fname = run_pname+hillslope_out_ext;
    string ft_properties_out_ext = "ft_properties.out";
    string ft_properties_out_fname = run_pname+ft_properties_out_ext;
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=


  ofstream zeta_out,h_out,eta_out,particle_out,hillslope_out,ft_properties_out,eroded_particle_out,surface_eroded_particle_out;
	zeta_out.open(zeta_out_fname.c_str());
	h_out.open(h_out_fname.c_str());
	eta_out.open(eta_out_fname.c_str());
    particle_out.open(particle_out_fname.c_str());
    eroded_particle_out.open(eroded_pout_fname.c_str());
	surface_eroded_particle_out.open(surface_eroded_pout_fname.c_str());
    hillslope_out.open(hillslope_out_fname.c_str());
    ft_properties_out.open(ft_properties_out_fname.c_str());


	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=
	// Initialize the variables
	//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

	// the number of partitions in the mixed and unmixed zone
    // This sets the template for how many finite volumes you will have throughout
	int n_PDZ_intervals;
	int n_CAZ_intervals;
	int ref_frame_switch;

	double SS_flux;										// flux
	double constant_surface_change_rate;
	vector<double> surf_erate;							// erosion from the surface
														// in m/yr: this is negative for
														// erosion
// These are variable for the new boundary condition
  double SS_erate;
  double ds_elev;
	int CRN_switch;		// sets the cosmogenics key:
										// 0 == no CRN
										// 1 == all CRN
										// 2 == all, neutron only
										// 3 == 10Be full
										// default == all, neutron only

	double t_ime;
	double start_time;
	int tt;					// a counter for the time

	int flux_switch;				// see flowtube.h: determines flux law
	int lower_boundary_condition;    // sets the lower boundary condition as either a flux (1) or elevation (2) values
    int prod_switch;				// see flowtube.h: determines soil production law
	double flux_us;					// see flowtube.h: flux from upslope
	double dt;						// time interval

	double start_depth;				// the starting depth in m
									// (LATER) it is set to 60 because
									// this is roughly 3x the e-folding
									// depth of the deepest muonogenic
									// production mechanism at rock density
	double vert_mix_vel;            // the vertical mixing velocity, where used in code is called vert_vel_fluctuating
	double horiz_mix_vel;           // horizontal mixing does not seem to be used anywhere
	double Omega;					// the activity of particles (that is the
									// proportion of particles that are moving
									// at any given time)

									// the elevation of the introduced particle
	double part_conc;				// particle concentration in particles per kg
	int part_switch;				//Particle insert switch, 1 = volumetric insertion, 2 = representative insertion


	vector<double> old_eta;			// the elevation of the soil-saprolite boudnary from the
									// last timestep
	vector<double> new_eta;			// the elevation of the soil-saprolite boundary from this
									// timestep

    // These a "bins" that contain particles. The plotting functions aggregate these bins and also they are
    // used to create average values.
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

	vector<double> starting_p_mfrac = vpi.get_type_mfracs();
	vector<int> starting_pID = vpi.get_type_index();
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  //
  // SET UP COSMO
  //
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
	// initial in situ concentrations in atoms per gram
	double C_10Be,C_f10Be,C_26Al,C_36Cl,C_14C,C_21Ne,C_3He;

	// set the parameters for production
	// 0 == granger
	// 1 == schaller
    // 2 == CRONUS
	int CRN_muon_param_switch;

	// the parameters for the in situ cosmogenics
	CRN_parameters CRNp;

	////////////////////////This part is superseded following the update to the cosmo code?
    // note: scaling determined using cosmocalc:
	// Schaller reports pinedale at 42 53 26 N
	// this converts to an inclination of 61.7 degrees
	// Schaller reports elevation of 2298 masl
	// this converts to an atmospheric depth of 781 g/cm2
	// this converts to a dunai scaling of 5.99
	// multiply this by schaller's snow shielding (for nucleonic production
	// of 0.925 one gets a scaling of 5.54
	double single_scaling;
    //////////////////////////
    // paramters for meteoric 10Be
	double M_supply_surface;	// in atoms/(cm^2*yr)
	double k_f10Be;				// in cm^2/g this is
									// an efolding depth of 20 g/cm^2
	double k2_f10Be;
	double chi_f10Be;			// fraction of supply that goes into the
								// shallow meteoric supply
	double deltad;				// in m (this gets converted
								// to cm in tParticle.cpp)
    //Added in to update the cosmo code to calculate specific site scaling factors
    double lon;                 // longitude
    double lat;                 // latitude
    double site_elev;           // site elevation
    double Fsp;                 // fsp is the fraction (between 0 and 1) of production at sea level
                                // and high latitude due to spallation (as opposed to muons).
                                // This argument is optional and defaults to 0.978, which is the value
                                // used by Stone (2000) for Be-10.


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
    //ofstream particle_out;
  //For the erosion rate changing
  int erate_end_time;
  //For the base level changes
  int base_level_end_time;
  double base_level_change;
	// get the model parameters from a file.
    // The parameter file includes the names of the parameters, these are read as
    // "temp" and then ignored after
	string temp;
	ifstream model_run_params_in;
	model_run_params_in.open(model_run_params_fname.c_str());
	model_run_params_in >> temp >> flux_switch >> temp >> prod_switch
					>> temp >> flux_us >> temp >> dt >> temp >> CRN_switch >> temp >> end_time
					>> temp >> constant_surface_change_rate >> temp >> particle_printing_interval
					>> temp >> eroded_catch_window >> temp >> max_age
					>> temp >> n_spacings >> temp >> insert_interval
          >> temp >> weathering_time_interval >> temp >> ref_frame_switch
          >> temp >> SS_flux >> temp >> lower_boundary_condition >> temp >> SS_erate >> temp >> ds_elev;
	model_run_params_in.close();
	cout << "LINE 209, got model_parameters" << endl;
	cout << "WTI: " << weathering_time_interval << " RFS: " << ref_frame_switch
	     << " SS_f: " << SS_flux << endl;

    cout << "BC=" << lower_boundary_condition <<endl;

    if (lower_boundary_condition ==1)
    {cout << "Flux boundary condition used" << endl;}
    else if (lower_boundary_condition ==2)
    {cout << "Elevation boundary condition used" << endl;}

	//cout << "Flux switch: " << flux_switch << " prod_switch: " << prod_switch
	//	 << " Flux_us: " << flux_us << " CRN switch: " << CRN_switch << endl
	//     << "dt: " << dt << " end_time: " << end_time << " surface erate: " << constant_surface_change_rate << endl
	//     << "particle print interval: " << particle_printing_interval
	//     << " catch: " << eroded_catch_window << " max age: " << max_age << " n_space: " << n_spacings
  //     << " particle insert interval: " << insert_interval
  //     << " weathering time interval: " << weathering_time_interval << endl
  //     << "Refernce time switch" << ref_frame_switch << endl;

	// get the parameters for the CRN particles
	ifstream CRN_parameter_in;
	CRN_parameter_in.open(CRN_parameter_fname.c_str());
	CRN_parameter_in >> temp >> start_depth >> temp >> vert_mix_vel
				     >> temp >> horiz_mix_vel >> temp >> Omega
				     >> temp >> part_conc >> temp >> part_switch >> temp >> CRN_muon_param_switch
				     >> temp >> single_scaling >> temp >> C_10Be >> temp >> C_f10Be >> temp >> C_26Al
				     >> temp >> C_36Cl >> temp >> C_14C >> temp >> C_21Ne >> temp
				     >> C_3He >> temp >> M_supply_surface >> temp >> k_f10Be >> temp
				     >> deltad >> temp >> k2_f10Be >> temp >> chi_f10Be
             >> temp >> n_PDZ_intervals >> temp >> n_CAZ_intervals >> temp >> lat >> temp >> lon >> temp >> site_elev >> temp >> Fsp;
	CRN_parameter_in.close();
	cout << "LINE 228, got CRN_parameters" << endl;

	cout << "start depth: " << start_depth << " vert mix vel: " << vert_mix_vel << endl
	    << "horiz mix vel: " << horiz_mix_vel << " Omega " << Omega << " part_conc: " << part_conc << " part_switch: " << part_switch << endl
	    << "Muon switch: " << CRN_muon_param_switch << " scaling: " << single_scaling << endl
	    << "Init conc, 10Be: " << C_10Be << " 26Al: " << C_26Al << " C_36Cl: " << C_36Cl << " C_14C: " << C_14C << endl
	    << "C_21Ne: " << C_21Ne << " C_3He: " << C_3He << " M_supp_surface: " << M_supply_surface << endl
	   << "K_f10Be: " << k_f10Be << " deltad: " << deltad << endl;
  // get the parameters for the CRN particles
  cout << "Getting the erosion history" << endl;

  ifstream erate_hist_in;
  erate_hist_in.open(erate_file_fname.c_str());
  erate_hist_in >> erate_end_time >> constant_surface_change_rate;
  // erate_hist_in.close();
  cout << "Got the erosion history, erate is:" << constant_surface_change_rate << "and will end at:"<< erate_end_time << endl;

  cout << "Base level lowering is: " << SS_erate << endl;

  ifstream base_level_hist_in;
  base_level_hist_in.open(base_level_hist_file_fname.c_str());
  base_level_hist_in >> base_level_end_time >> base_level_change;
  cout << "Got the, base level change is:" << base_level_change << "and will end at:"<< base_level_end_time << endl;


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
//	int tot_intervals = n_PDZ_intervals+n_CAZ_intervals;
//
//  CRUNCH_bins Geochem_bins(n_PDZ_intervals , n_CAZ_intervals, start_depth, vpi);
//  //cout << "Now getting "
//
//  // create the crunch engine
//	string master_fname = "master_crunch.in";
//	CRUNCH_engine Ceng(crunch_pname, run_pname, master_fname);
//	cout << "LINE 306 Got crunch enginge" << endl;
//
//
//	string value_name = "spatial_profile";
//	Ceng.modify_CRUNCH_value(value_name, double(weathering_time_interval));
  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


	// initialize the CRN parameters
	if (CRN_muon_param_switch == 1)
	{
		CRNp.set_Schaller_parameters();
        cout << "scaled to schaller" << endl;
	}
    else if (CRN_muon_param_switch == 2)
    {
        CRNp.set_newCRONUS_parameters();
        cout << "scaled to CRONUS" << endl;
    }
    vector<bool> nuclides_for_scaling;
	CRNp.scale_F_values(nuclides_for_scaling, lat, lon, site_elev, Fsp);


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

	/// raise the flowtube so the downstream boundary is at elevation
	/// 100, this is done to avoid negative elevations due to lowering over time however this here is set to zeta_zero
	vector<double> zeta_old = ft_test.get_zeta();
	double zeta_zero = 1000.0;
	ft_test.raise_zeta_eta_ds_bound(zeta_zero);

	double rho_s = ft_test.get_rho_s();

	double insert_time_clock = 0;
	old_eta = ft_test.get_eta();
	vector<double> old_zeta = ft_test.get_zeta();
	vector<double> Delta_zeta = old_zeta;
	vector<double> old_bottom_depth = old_zeta;
	vector<double> old_h = ft_test.get_h();
	vector<double> last_insertion_zeta = ft_test.get_zeta();
	vector<double> this_insertion_zeta = ft_test.get_zeta();
	int zeta_sz = Delta_zeta.size();
	int part_ID_start = 1;

    /// Set the ds elevation for elevation boundary condition
    double ds_elevation = ft_test.get_zeta_ds();
    cout << "got ds_elevation: " << ds_elevation << endl;

	// reset the time
	t_ime = 0;
	start_time = 0;
	tt = 0;
	double old_t_time = 0;
  double weathering_time_clock = 0;


  // for the insertion of particles, we need to ensure that the insertion occurs
  // BELOW the zone in which CrunchFlow is evaluated. This means particles should
  // not be inserted within the bottom few crunch cells. We can calculate this
  // by compring the particle insert interval with the erosion rate.
  double depth_between_insertion = insert_interval*(-constant_surface_change_rate);

  // add a little fudge factor (10%) to this just in case the eta lowering
  // is faster than surface lowering with transient soil thickness
  depth_between_insertion= depth_between_insertion*1.1;
  cout << "Added depth below start depth: " << depth_between_insertion << endl;

	// insert initial particles
	// to insert particles throughout the domain we need to set old botom depth as the
	// zeta elevations.
	for (int i = 0; i<zeta_sz; i++)
	{
			old_bottom_depth[i] = old_zeta[i];
			last_insertion_zeta[i] = old_zeta[i];
			Delta_zeta[i] = start_depth*1.25;
			//cout << "LINE 350 h["<<i<<"]: " << old_h[i] << " and zeta: " << old_zeta[i] << endl
			//     << " and start depth: " << start_depth << " and Delta_eta: " << Delta_eta[i] << endl;
	}

	switch ( part_switch )
		{case 1:
      	part_ID_start = CRN_tpb.insert_particles_volumetric(ft_test, Delta_zeta, old_bottom_depth,
										C_10Be, C_f10Be, C_26Al, C_36Cl, C_14C, C_21Ne, C_3He, vpi);
		break;
		case 2 :								
	  	part_ID_start = CRN_tpb.insert_particles(ft_test, Delta_zeta, old_bottom_depth, part_conc, starting_pID, starting_p_mfrac,
 										C_10Be, C_f10Be, C_26Al, C_36Cl, C_14C, C_21Ne, C_3He);
		break;								 									
		}

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

    vector<double> bc_h_temp = ft_test.get_zeta();
      vector<double> bc_h = bc_h_temp;

		// ds_elev -= dt*SS_erate;

		/// Runs a flux timestep
		//cout << "running flux ";
        // cout << "SS_flux is: " << SS_flux << endl;
    if (t_ime == base_level_end_time)
          { cout << "Changing the base level rate" << endl;
          base_level_hist_in >> base_level_end_time >> base_level_change;

          cout << "Got the  history, new base level change is:" << base_level_change << "and will end at:"<< base_level_end_time << endl;
        }


    if (t_ime == erate_end_time)
      { cout << "Changing the erate" << endl;
      erate_hist_in >> erate_end_time >> constant_surface_change_rate;

      cout << "Got the erosion history, new erate is:" << constant_surface_change_rate << "and will end at:"<< erate_end_time << endl;
      surf_erate.clear();
      int n_ft_nodes = ft_test.get_n_nodes();
      	for (int node = 0; node < n_ft_nodes; node++)
      	{
      		surf_erate.push_back(constant_surface_change_rate);
      	}
      	cout << "Changed the erate successfully" << endl;
      }
        ///LK 21/8/18 Implementing a statement to allow both a flux and elevation boundary condition
        ///Flux = 1
        ///Elevation = 2
		if (lower_boundary_condition == 1)
        {    ft_test.flux_timestep_flux_bc(dt, flux_us, SS_flux,flux_switch, prod_switch,
							surf_erate);

		}
        else if (lower_boundary_condition == 2)
        {  ft_test.flux_timestep_elev_bc(dt,
							flux_us, ds_elev,
							 flux_switch,  prod_switch,
							 surf_erate);
            	ds_elev -= dt*base_level_change;
              cout << "New ds_elev is: "<< ds_elev << endl;
        }
        else if (lower_boundary_condition == 3)
        {  ft_test.flux_timestep_varying_elev_bc(dt,
              flux_us, ds_elev,
               flux_switch,  prod_switch,
               surf_erate,constant_surface_change_rate,base_level_change);

               double ds_elev_temp = ft_test.update_ds_elev(dt, base_level_change, prod_switch,
                      constant_surface_change_rate, ds_elev, bc_h);
                       ds_elev= ds_elev_temp;
                         cout << "New ds_elev is: "<< ds_elev << endl;


        }


                             // }
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
			// if (t_ime >= old_t_time+part_p_i-eroded_catch_window && t_ime <= old_t_time+part_p_i)
			// {
			// 	//cout << "LINE 423, time is: " << t_ime << " and next_catch: " << next_catch << endl;
			// 	for(int bn = 0; bn<=n_bins; bn++)
			// 	{
			// 		if (eroded_bins[bn].size()>0)
			// 		{
			// 			part_iter = eroded_bins[bn].begin();
			// 			while(part_iter != eroded_bins[bn].end())
			// 			{
			// 				eroded_catcher[bn].push_back(*part_iter);
			// 				part_iter++;
			// 			}
			// 		CRN_tpb.print_eroded_stats(t_ime,eroded_bins , ft_test, eroded_particle_out);	
			// 		}
			// 		if (t_ime == old_t_time+part_p_i)
			// 			old_t_time = t_ime;
			// 	}
			// }

			// Setup window to catch particles
			if (t_ime >= old_t_time+part_p_i-eroded_catch_window && t_ime <= old_t_time+part_p_i+eroded_catch_window)
			{
				double print_time =	old_t_time+part_p_i;
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
					
					CRN_tpb.print_surface_eroded_particles(print_time,eroded_bins , ft_test, surface_eroded_particle_out);
					}
					if (t_ime == old_t_time+part_p_i+eroded_catch_window)
						old_t_time = t_ime-eroded_catch_window;
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
			this_insertion_zeta = ft_test.get_zeta(); 	// get the updated eta

			// cout << "Time: " << t_ime << " Delta eta: " << endl;
			// calculate delta eta
			for(int ii = 0; ii< zeta_sz; ii++)
			{
				Delta_zeta[ii] = last_insertion_zeta[ii]-this_insertion_zeta[ii];
				// cout << "i = " << ii << " Delta_zeta: " << Delta_zeta[ii] << endl;
			}

			// insert the particles
		switch ( part_switch )
		{case 1:
      	part_ID_start = CRN_tpb.insert_particles_volumetric(ft_test, Delta_zeta, old_bottom_depth,
										C_10Be, C_f10Be, C_26Al, C_36Cl, C_14C, C_21Ne, C_3He, vpi);
		break;
		case 2 :								
	  	part_ID_start = CRN_tpb.insert_particles(ft_test, Delta_zeta, old_bottom_depth, part_conc, starting_pID, starting_p_mfrac,
 										C_10Be, C_f10Be, C_26Al, C_36Cl, C_14C, C_21Ne, C_3He);
		break;								 									
		}
			last_insertion_zeta = this_insertion_zeta;				// reset old eta
			insert_time_clock = 0;			// reset the insert time clock
		}               // !end particle insertion




		//cout << "running printing" << endl;
		// prints various data to file
		if (particle_trigger == 1 && tt%part_p_i== 0)
		{
			cout << "LINE 440, printing particles, time: " << t_ime << endl;
			ft_test.print_eta(t_ime, eta_out);
			ft_test.print_h(t_ime, h_out);
            ft_test.print_zeta(t_ime, zeta_out);
            ft_test.export_input_profile(hillslope_out);
			ft_test.print_ft_properties(ft_properties_out);
            CRN_tpb.print_particle_stats(t_ime, ft_test, particle_out);
            CRN_tpb.print_eroded_stats(t_ime,eroded_bins , ft_test, eroded_particle_out);
            //int ref_frame_switch = 1;

			// print basic particle information
//			CRN_tpb.vtk_print_basic_volume_particles(t_ime, vtk_particle_fname,
//                                               ref_frame_switch);

            // print cosmo data
//            CRN_tpb.print_particle_stats_vtk(t_ime, ft_test,
//								vtk_fname);
//


            // print cell data
//            Geochem_bins.vtk_cell_bundler(t_ime, ref_frame_switch,
//                                    vtk_cell_fname, Ceng, CRN_tpb);

		}
		//cout << "...ran printing" << endl;

		if (tt%2500 == 0)
		{
			cout << "time is: " << t_ime << endl;

		}
  }
  // Not sure this erate clse is needed
  erate_hist_in.close();
  base_level_hist_in.close();
  eta_out.close();
  h_out.close();
  zeta_out.close();
  particle_out.close();
  eroded_particle_out.close();
  surface_eroded_particle_out.close();
  hillslope_out.close();
  ft_properties_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

