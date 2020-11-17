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

#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include "../tParticle.hpp"
#include "../CRN_tParticle_bins.hpp"
#include "../CRUNCH_engine.hpp"
#include "../flowtube.hpp"
#include "../FT_util.hpp"
#include "../VolumeParticleInfo.hpp"
#include "../TNT/tnt.h"
#include "../LSDStatsTools.hpp"
#include "../parameterparser.hpp"

using namespace std;

int main (int nNumberofArgs,char *argv[])
{

  // Get the arguments
  vector<string> path_and_file = DriverIngestor(nNumberofArgs,argv);


  string path_name = path_and_file[0];
  string f_name = path_and_file[1];

  // load parameter parser object
  parameterparser LSDPP(path_name,f_name);

// maps for setting default parameters
  map<string,int> int_default_map;
  map<string,float> float_default_map;
  map<string,bool> bool_default_map;
  map<string,string> string_default_map;

//Set up the default parameters
///Model run parameters
int_default_map["flux_switch"] = ;
int_default_map["prod_switch"] = ;
float_default_map["flux_us"] = ;
int_default_map["dt"] = ;
int_default_map["CRN_switch"] = ;
int_default_map["end_time"] = ;
float_default_map["surf_erate"] = ;
int_default_map["particle_printing_interval"] = ;
int_default_map["eroded_catch_interval"] = ;
int_default_map["max_age"] = ;
int_default_map["n_spacings"] = ;
int_default_map["particle_insert_interval"] = ;
int_default_map["weathering_time_interval"] = ;
int_default_map["ref_frame_switch"] = ;
float_default_map["SS_flux"] = ;
int_default_map["Lower_boundary_condition"] = ;
///Crunch Switch - add in a bit to read in if a crunchflow input file is present then use it!
bool_default_map["Crunch_switch"] = false;
///CRN parameters
float_default_map["start_depth"] = ;
float_default_map["vert_mix_vel"] = 0;
float_default_map["horiz_mix_vel"] = 0;
float_default_map["Omega"] = 0.5;
float_default_map["part_conc"] = ;
int_default_map["CRN_muon_param_switch"] = 2;
float_default_map["single_scaling"] = ;
float_default_map["C_10Be_initial"] = 0;
float_default_map["C_26AL_initial"] = 0;
float_default_map["C_36Cl_initial"] = 0;
float_default_map["C_14C_initial"] = 0;
float_default_map["C_21Ne_initial"] = 0;
float_default_map["C_3He_initial"] = 0;
float_default_map["M_supply_surface"] = ;
float_default_map["k_f10Be"] = ;
float_default_map["k2_f10Be"] = ;
float_default_map["Chi_f10Be"] = ;
int_default_map["delta_d"] = ;
int_default_map["n_PDZ_intervals"] = ;
int_default_map["n_CAZ_intervals"] = ;
float_default_map["lat"] = ;
float_default_map["lon"] = ;
float_default_map["site_elev"] = ;
float_default_map["Fsp"] = ;


///Sediment transport parameters
float_default_map["rho_s"] = ;
float_default_map["rho_r"] = ;
float_default_map["K_h"] = ;
float_default_map["S_c"] = ;
float_default_map["W_0"] = ;
float_default_map["Gamma"] = ;
float_default_map["N"] = ;
float_default_map["N_0"] = ;
float_default_map["N_m"] = ;
float_default_map["Beta"] = ;
float_default_map["K_g"] = ;

