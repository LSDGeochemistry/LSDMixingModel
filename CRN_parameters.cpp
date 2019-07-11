//tParticle.cpp
// implementation of the dParticle object
#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include "mathutil.hpp"
#include "CRN_parameters.hpp"
#include "LSDStatsTools.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// the CRN_parameters object
//
// Sets constants based on the Vermeesh (2007) paper adjusted for both the Granger and Schaller
// formulations. Uses a single given scaling factor to produce scaling factors for each production method.
// Does this individually for every CRN.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// function to set CRN parameters
void CRN_parameters::create()
{
	S_t = 1;

	//decay constants from Vermeesh 2007
	lambda_10Be = 456e-9;		// in yr-1
	lambda_26Al = 980e-9;		// in yr-1
	lambda_14C = 121e-6;		// in yr-1
	lambda_36Cl = 230e-8;		// in yr-1

    //production rate at sea level, high latitude
	// from Vermeesh 2007
	P0_10Be = 5.11;					// in a/g/yr
	P0_26Al = 30.31;				// in a/g/yr
	P0_14C = 5.86;					// in a/g/yr
	P0_36Cl = 55.45;				// in a/g/yr
	P0_21Ne = 20.29;				// in a/g/yr
	P0_3He = 97.40;					// in a/g/yr

	//attenutation lengths in g/cm^2
	Gamma[0] = 160;
	Gamma[1] = 738.6;
	Gamma[2] = 2688;
	Gamma[3] = 4360;

	
    // Relative Productions by different sources for each nuclide
    // dimensionless
	F_10Be[0] = 0.9724;
	F_10Be[1] = 0.0186;
	F_10Be[2] = 0.004;
	F_10Be[3] = 0.005;

	// dimensionless
	F_26Al[0] = 0.9655;
	F_26Al[1] = 0.0233;
	F_26Al[2] = 0.005;
	F_26Al[3] = 0.0062;

	// dimensionless
	F_14C[0] = 0.83;
	F_14C[1] = 0.0691;
	F_14C[2] = 0.0809;
	F_14C[3] = 0.02;

	// dimensionless
	F_36Cl[0] = 0.903;
	F_36Cl[1] = 0.0447;
	F_36Cl[2] = 0.05023;
	F_36Cl[3] = 0.0;
}

void CRN_parameters::set_Granger_parameters()
{
	S_t = 1;

	// from Vermeesh 2007
	lambda_10Be = 456e-9;		// in yr-1
	lambda_26Al = 980e-9;		// in yr-1
	lambda_14C = 121e-6;		// in yr-1
	lambda_36Cl = 230e-8;		// in yr-1

	// from Vermeesh 2007
	P0_10Be = 5.11;					// in a/g/yr
	P0_26Al = 30.31;				// in a/g/yr
	P0_14C = 5.86;					// in a/g/yr
	P0_36Cl = 55.45;				// in a/g/yr
	P0_21Ne = 20.29;				// in a/g/yr
	P0_3He = 97.40;					// in a/g/yr

	// in g/cm^2
	Gamma[0] = 160;
	Gamma[1] = 738.6;
	Gamma[2] = 2688;
	Gamma[3] = 4360;

	// dimensionless
	F_10Be[0] = 0.9724;
	F_10Be[1] = 0.0186;
	F_10Be[2] = 0.004;
	F_10Be[3] = 0.005;

	// dimensionless
	F_26Al[0] = 0.9655;
	F_26Al[1] = 0.0233;
	F_26Al[2] = 0.005;
	F_26Al[3] = 0.0062;

	// dimensionless
	F_14C[0] = 0.83;
	F_14C[1] = 0.0691;
	F_14C[2] = 0.0809;
	F_14C[3] = 0.02;

	// dimensionless
	F_36Cl[0] = 0.903;
	F_36Cl[1] = 0.0447;
	F_36Cl[2] = 0.05023;
	F_36Cl[3] = 0.0;
}

// function to set CRN parameters
// based on the Vermeesch approximation of the Schaller et al (2000)
// formulation
void CRN_parameters::set_Schaller_parameters()
{
	S_t =1;

	// from Vermeesh 2007
	lambda_10Be = 456e-9;		// in yr-1
	lambda_26Al = 980e-9;		// in yr-1
	lambda_14C = 121e-6;		// in yr-1
	lambda_36Cl = 230e-8;		// in yr-1

	// from Vermeesh 2007
	P0_10Be = 5.11;					// in a/g/yr
	P0_26Al = 30.31;				// in a/g/yr
	P0_14C = 5.86;					// in a/g/yr
	P0_36Cl = 55.45;				// in a/g/yr
	P0_21Ne = 20.29;				// in a/g/yr
	P0_3He = 97.40;					// in a/g/yr

	// in g/cm^2
	Gamma[0] = 160;
	Gamma[1] = 738.6;
	Gamma[2] = 2688;
	Gamma[3] = 4360;

	// dimensionless
	F_10Be[0] = 0.964;
	F_10Be[1] = 0.0266;
	F_10Be[2] = -0.0074;
	F_10Be[3] = 0.0168;

	// dimensionless
	F_26Al[0] = 0.9575;
	F_26Al[1] = 0.0315;
	F_26Al[2] = -0.009;
	F_26Al[3] = 0.02;

	// dimensionless
	F_14C[0] = 0.83;
	F_14C[1] = 0.1363;
	F_14C[2] = 0.0137;
	F_14C[3] = 0.02;

	// dimensionless
	F_36Cl[0] = 0.903;
	F_36Cl[1] = 0.0793;
	F_36Cl[2] = 0.0177;
	F_36Cl[3] = 0.0;
}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// 10Be is set to a new production curve provided by Shasta Marrero
// All others: sets the parameters to those used by Braucher et al 2009
// as implemented in cosmocalc v2.0
// http://www.ucl.ac.uk/~ucfbpve/cosmocalc/updates.html
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRN_parameters::set_newCRONUS_parameters()
{
  S_t = 1;

  // from Vermeesh 2007
  // 10Be from Chmeleff/Korschinek 10Be decay constant;
  lambda_10Be = 500e-9;    // in yr-1
  lambda_26Al = 980e-9;    // in yr-1
  lambda_14C = 121e-6;     // in yr-1
  lambda_36Cl = 230e-8;    // in yr-1          

  // from Data provided by Shasta mararreo for 10Be, everyting
  // else is from the Braucher constants
  //
  // All but 10Be are calibrated to the Stone scaling
  // Also linke to the nishizumii standards
  // These come with Cosmocalc version 2.0
  // http://www.ucl.ac.uk/~ucfbpve/cosmocalc/updates.html
  P0_10Be = 4.075213;          // in a/g/yr
  P0_26Al = 31.10;         // in a/g/yr
  P0_14C = 15.21;          // in a/g/yr
  P0_36Cl = 58.95;         // in a/g/yr
  P0_21Ne = 18.23;         // in a/g/yr
  P0_3He = 121.59;         // in a/g/yr

  // in g/cm^2
  Gamma[0] = 160;
  Gamma[1] = 1459.76761923;
  Gamma[2] = 11039.2402217;
  Gamma[3] = 4320;

  // dimensionless
  F_10Be[0] = 0.98374;
  F_10Be[1] = 0.0137188126531;
  F_10Be[2] = 0.00252519164093;
  F_10Be[3] = 0.0;

  // dimensionless
  F_26Al[0] = 0.9699;
  F_26Al[1] = 0.0275;
  F_26Al[2] = 0.000;
  F_26Al[3] = 0.0026;

  // dimensionless
  F_14C[0] = 0.83;
  F_14C[1] = 0.15;
  F_14C[2] = 0.0;
  F_14C[3] = 0.02;

  // dimensionless
  F_36Cl[0] = 0.9456;
  F_36Cl[1] = 0.0324;
  F_36Cl[2] = 0.00;
  F_36Cl[3] = 0.022;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//Implementing new way to calculate scaling function, first find pressure as function of elevation,
//then find a scaling factor using stone 2000 scaling factor. Then update F values along with update 10Be
//with new depths.

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function converts elevation to atmospheric pressure
// It follows the cronus calculator NCEPatm_2 function
// written by Greg Balco
// SMM
// 4/12/2014
//
// Looks up surface pressure and 1000 mb temp from NCEP reanalysis
// and calculates site atmospheric pressures using these as inputs to the
// standard atmosphere equation. 
//
// Syntax: pressure = NCEPatm_2(lat,lon,site_elv);
// 
// Requires:
///       lat: latitude (DD). Southern hemisphere is negative.
//       lon: longitude (DD). Western hemisphere is negative.
//           Tries to deal with 0-360 longitudes gracefully.
//       site_elv: elevation (m).
//
// Returns site pressure in hPa.
//
// Vectorized. Send vectors of equal length.
//
// Note: this must load the data file NCEP2.mat whenever called. 
// Repeated calls to this function will be slow for this reason. 
//
// Also: This function is OK but not great for Antarctica.
// Use antatm.m instead. 
//
// Remember: it is always better to estimate the average pressure at your 
// site using a pressure-altitude relation obtained from nearby station
// data.
//
// Original m code Written by Greg Balco -- UW Cosmogenic Nuclide Lab
// balcs@u.washington.edu
// October, 2007
// Part of the CRONUS-Earth online calculators: 
//      http://hess.ess.washington.edu/math
//
// Copyright 2001-2007, University of Washington
// All rights reserved
// Developed in part with funding from the National Science Foundation.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double CRN_parameters::NCEPatm_2(double lat, double lon, double site_elev)
{
  // deal with negative longitudes
  if(lon < 0)
  {
    lon = lon+360.0;
  }
  
  //cout << "LSDCRNP, line 1618, Site lat: " << lat << " and long: " << lon << endl;
  
  // check to see if data is loaded:
  if (int(gm_hgt.size()) != 8)
  {
    //Assumes NCEP data is in directory above. maybe this should be included in main mixing_column.cpp  
    string path_to_NCEP_data = "./";
//      cout << "You didn't load the NCEP data. Doing that now. " << endl;
//      cout << "Enter path to data files: " << endl;
//      cin >> path_to_NCEP_data;
    load_parameters_for_atmospheric_scaling(path_to_NCEP_data);
  }
  
  // now, interpolate sea level pressure and temperature
  //cout << "interpolating pressure " << endl;
  double site_slp = interp2D_bilinear(NCEPlat, NCEPlon, meanslp, 
                        lat, lon);
  //cout << "Did pressure, now temp:" << endl;                      
  double site_T = interp2D_bilinear(NCEPlat, NCEPlon, meant1000, 
                        lat, lon);
  //cout << "Did temp" << endl;                      
  
  
  double site_T_degK = site_T + 273.15;

  // Some More parameters
  double gmr = -0.03417; // Assorted constants (this has come from Greg Balco's code)
  double dtdz = 0.0065;  // Lapse rate from standard atmosphere 
  
  // Calculate site pressure using the site-specific SLP and T1000 with the
  // standard atmosphere equation.

  //cout << site_T_degK << endl;
  //cout << "Log term: " << log(site_T_degK) - log(site_T_degK - (site_elev*dtdz)) << endl;
  //cout << "Exp term: " << exp( (gmr/dtdz)*( log(site_T_degK) - log(site_T_degK - (site_elev*dtdz)) ) ) << endl;
//    cout <<"latitude is: " << lat << endl;
//    cout <<"longitude is: " << lon << endl;
//cout << "site elevation is: " << site_elev << endl;
  double out = site_slp*exp( (gmr/dtdz)*( log(site_T_degK) - log(site_T_degK - (site_elev*dtdz)) ) );
  
  //cout << endl;
//  cout << "Site sea level pressure: " << site_slp << " and site Temp: "<< site_T 
//      << " and pressure: " << out << endl << endl <<endl;
  return out;
}

// this function gets the parameters used to convert elevation to 
// pressure
void CRN_parameters::load_parameters_for_atmospheric_scaling(string path_to_NCEP_data)
{
  cout.precision(8);
  
  // first load the levels
  levels.push_back(1000);
  levels.push_back(925);
  levels.push_back(850);
  levels.push_back(700);
  levels.push_back(600);
  levels.push_back(500);
  levels.push_back(400);
  levels.push_back(300);
  
  // the dimensions of the data
  int n_levels = 8;
  int NRows = 73;
  int NCols = 145;
  Array2D<double> new_slp(NRows,NCols,0.0);
  Array2D<double> new_meant(NRows,NCols,0.0);
    
  // now load the mean sea level pressure
  string filename = "NCEP2.bin";
  filename = path_to_NCEP_data+filename;
  //cout << "Loading mean sea level, file is: " << endl << filename << endl;

  ifstream ifs_data(filename.c_str(), ios::in | ios::binary);
  if( ifs_data.fail() )
  {
    cout << "\nFATAL ERROR: the data file \"" << filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }  
  

  double temp;
  //cout << "The size of a double is: " << sizeof(temp) << endl;
  for (int i=0; i<NCols; ++i)
  {
    for (int j=0; j<NRows; ++j)
    {
      ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
      new_slp[j][i] = temp;
      //cout << "new_slp["<<j+1<<"]["<<i+1<<"]: " << new_slp[j][i] << endl;
    }
  }
  
  for (int i=0; i<NCols; ++i)
  {
    for (int j=0; j<NRows; ++j)
    {
      ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
      new_meant[j][i] = temp;
      //cout << "new_meant100["<<j+1<<"]["<<i+1<<"]: " << new_meant[j][i] << endl;
    }
  }  
  
  // now get the indices
  vector<double> temp_lat(NRows,0.0);
  for (int i=0; i<NRows; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    temp_lat[i] = temp;
    //cout << "Lat["<<i+1<<"]: " << temp_lat[i] << endl;
  }
  vector<double> temp_long(NCols,0.0);
  for (int i=0; i<NCols; ++i)
  {
    ifs_data.read(reinterpret_cast<char*>(&temp), sizeof(temp));
    temp_long[i] = temp;
    //cout << "Long["<<i+1<<"]: " << temp_long[i] << endl;
  }  
  
  ifs_data.close();
  
  
  // now the data with levels
  filename = "NCEP_hgt.bin";
  filename = path_to_NCEP_data+filename;
  //cout << "Loading hgt, file is: " << endl << filename << endl;

  ifstream ifs_data2(filename.c_str(), ios::in | ios::binary);
  if( ifs_data2.fail() )
  {
    cout << "\nFATAL ERROR: the data file \"" << filename
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }  
  

  // get the gm heights
  vector< Array2D<double> > vec_hgt_gm_array;
  for (int lvl = 0; lvl < n_levels; lvl++)
  {
    Array2D<double> current_hgt_array(NRows,NCols,0.0);
    for (int i=0; i<NCols; ++i)
    {
      for (int j=0; j<NRows; ++j)
      {
        ifs_data2.read(reinterpret_cast<char*>(&temp), sizeof(temp));
        current_hgt_array[j][i] = temp;
        //cout << "new_slp["<<lvl<<"]["<<j+1<<"]["<<i+1<<"]: " << current_hgt_array[j][i] << endl;
      }
    }
    //cout << "new_slp["<<j+1<<"]["<<i+1<<"]: " << new_slp[j][i] << endl;
    vec_hgt_gm_array.push_back(current_hgt_array.copy());
  }

  // now the gp heights
  vector< Array2D<double> > vec_hgt_gp_array;
  for (int lvl = 0; lvl < n_levels; lvl++)
  {
    Array2D<double> current_hgt_array(NRows,NCols,0.0);
    for (int i=0; i<NCols; ++i)
    {
      for (int j=0; j<NRows; ++j)            
      {
        ifs_data2.read(reinterpret_cast<char*>(&temp), sizeof(temp));
        current_hgt_array[j][i] = temp;
        //cout << "new_slp["<<lvl<<"]["<<j+1<<"]["<<i+1<<"]: " << current_hgt_array[j][i] << endl;
      }
    }
    //cout << "new_slp["<<j+1<<"]["<<i+1<<"]: " << new_slp[j][i] << endl;
    vec_hgt_gp_array.push_back(current_hgt_array.copy());
  }
  ifs_data2.close();
  
  // now update the data elements
  NCEPlat = temp_lat;
  NCEPlon = temp_long;
  
  meanslp = new_slp.copy();
  meant1000 = new_meant.copy();
  gp_hgt = vec_hgt_gp_array;
  gm_hgt = vec_hgt_gp_array;  
  
  //cout << "Size lat: " << NCEPlat.size() << " size long: " << NCEPlon.size() << endl;
  //cout << "size slp:" << meanslp.dim1() << " " << meanslp.dim2() << endl;
  //cout << "size t1000: " << meant1000.dim1() << " " << meant1000.dim2() << endl;
  

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Scaling from the Stone 2000 paper
// Units:
// latitude in decimal degrees
// pressure in hPa
// fsp is the fraction (between 0 and 1) of production at sea level
// and high latitude due to spallation (as opposed to muons).
// This argument is optional and defaults to 0.978, which is the value
// used by Stone (2000) for Be-10. The corresponding value for Al-26
// is 0.974. Note that using 0.844 for Be-10 and 0.826 for Al-26 will
// closely reproduce the Lal, 1991 scaling factors as long as the standard
// atmosphere is used to convert sample elevation to atmospheric pressure.
// Also note that this function will yield the scaling factor for spallation
// only when fsp=1, and that for muons only when fsp=0.
//
// IMPORTANT: This (and the Rc version) is probably the best scaling method!
// See https://cosmognosis.wordpress.com/2014/01/07/high-altitude-low-latitude-calibration-sites-i/
//
// Elevation can be converted to pressure with the functions
// stdatm.m (general use) and antatm.m (Antarctica).
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
double CRN_parameters::stone2000sp(double lat, double lon, double site_elev, double Fsp)
{
    
    string path_to_NCEP_data = "/driver_functions/";
    //Find the pressure 
    double P = NCEPatm_2(lat, lon, site_elev);
    cout << "p is: " << P << endl;
    
  if (Fsp > 1)
  {
    Fsp = 0.978;
  }
  
  if (fabs(lat) > 90)
  {
    cout << "Your latitude is > 90! Defaulting to 45 degrees" << endl;
    lat = 45;
  }

  // Spallogenic production at index latitudes;
  // enter constants from Table 1
  vector<double> a;
  a.push_back(31.8518);
  a.push_back(34.3699);
  a.push_back(40.3153);
  a.push_back(42.0983);
  a.push_back(56.7733);
  a.push_back(69.0720);
  a.push_back(71.8733);

  vector<double> b;
  b.push_back(250.3193);
  b.push_back(258.4759);
  b.push_back(308.9894);
  b.push_back(512.6857);
  b.push_back(649.1343);
  b.push_back(832.4566);
  b.push_back(863.1927);

  vector<double> c;
  c.push_back(-0.083393);
  c.push_back(-0.089807);
  c.push_back(-0.106248);
  c.push_back(-0.120551);
  c.push_back(-0.160859);
  c.push_back(-0.199252);
  c.push_back(-0.207069);

  vector<double> d;
  d.push_back(7.4260e-5);
  d.push_back(7.9457e-5);
  d.push_back(9.4508e-5);
  d.push_back(1.1752e-4);
  d.push_back(1.5463e-4);
  d.push_back(1.9391e-4);
  d.push_back(2.0127e-4);

  vector<double> e;
  e.push_back(-2.2397e-8);
  e.push_back(-2.3697e-8);
  e.push_back(-2.8234e-8);
  e.push_back(-3.8809e-8);
  e.push_back(-5.0330e-8);
  e.push_back(-6.3653e-8);
  e.push_back(-6.6043e-8);

  vector<double> ilats;
  ilats.push_back(0);
  ilats.push_back(10);
  ilats.push_back(20);
  ilats.push_back(30);
  ilats.push_back(40);
  ilats.push_back(50);
  ilats.push_back(60);

  // calculate index latitudes at given P's
  double lat0  = a[0] + (b[0] * exp(P/(-150.0))) + (c[0]*P) + (d[0]*(P*P)) + (e[0]*(P*P*P));
  double lat10 = a[1] + (b[1] * exp(P/(-150.0))) + (c[1]*P) + (d[1]*(P*P)) + (e[1]*(P*P*P));
  double lat20 = a[2] + (b[2] * exp(P/(-150.0))) + (c[2]*P) + (d[2]*(P*P)) + (e[2]*(P*P*P));
  double lat30 = a[3] + (b[3] * exp(P/(-150.0))) + (c[3]*P) + (d[3]*(P*P)) + (e[3]*(P*P*P));
  double lat40 = a[4] + (b[4] * exp(P/(-150.0))) + (c[4]*P) + (d[4]*(P*P)) + (e[4]*(P*P*P));
  double lat50 = a[5] + (b[5] * exp(P/(-150.0))) + (c[5]*P) + (d[5]*(P*P)) + (e[5]*(P*P*P));
  double lat60 = a[6] + (b[6] * exp(P/(-150.0))) + (c[6]*P) + (d[6]*(P*P)) + (e[6]*(P*P*P));

  vector<double> lat_at_specifics(7,0.0);
  lat_at_specifics[0] = lat0;
  lat_at_specifics[1] = lat10;
  lat_at_specifics[2] = lat20;
  lat_at_specifics[3] = lat30;
  lat_at_specifics[4] = lat40;
  lat_at_specifics[5] = lat50;
  lat_at_specifics[6] = lat60;

  //northernize southern-hemisphere inputs
  lat = fabs(lat);

  //set high lats to 60;
  if(lat > 60)
  {
    lat = 60.0;
  }

  // interpoloate elevation
  double S = interp1D_ordered(ilats,lat_at_specifics, lat);

  // Production by muons

  //constants
  vector<double> mk;
  mk.push_back(0.587);
  mk.push_back(0.600);
  mk.push_back(0.678);
  mk.push_back(0.833);
  mk.push_back(0.933);
  mk.push_back(1.000);
  mk.push_back(1.000);

  // index latitudes at given P's
  vector<double> m_index_at_given_P;
  m_index_at_given_P.push_back(mk[0]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[1]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[2]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[3]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[4]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[5]*exp( (1013.25-P)/242.0));
  m_index_at_given_P.push_back(mk[6]*exp( (1013.25-P)/242.0));
   
  // interpolate for actual elevation
  double M = interp1D_ordered(ilats, m_index_at_given_P, lat);

  //cout << "S: " << S << " M: " << M << " Fsp: " << Fsp << endl;

  // Combine spallogenic and muogenic production; return
  double Fm = 1 - Fsp
      ;
  double out = ((S * Fsp) + (M * Fm));
  
  //cout << "Stone 2000 scaling is: "  << out << endl;

  return out;
}




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function is similar to the scaling function for total nuclides but
// it allows the user to pick the nuclides they want scaled, and
// also uses a faster newton raphson iterator to get the correct scaling
//
// The bool vector nuclides_for_scaling has four elements
// nuclides_for_scaling[0] = true: calculate 10Be
// nuclides_for_scaling[1] = true: calculate 26Al
// nuclides_for_scaling[2] = true: calculate 36Cl
// nuclides_for_scaling[3] = true: calculate 14C
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRN_parameters::scale_F_values(vector<bool> nuclides_for_scaling, double lat, double lon, double site_elev, double Fsp)
{
  // first check the boolean vector. If it is an incorrect size, print a warning
  // and then default to all true
  if(nuclides_for_scaling.size() != 4)
  {
    cout << "nulides for scaling vector is the wrong size! Defaulting to a " << endl;
    cout << "vector that calculates all nuclides" << endl;
    vector<bool> temp_vec(4,true);
    nuclides_for_scaling = temp_vec;
  }

  // set up the parameters for the newton-raphson iteration
  double single_scaling = stone2000sp(lat, lon, site_elev, Fsp);    
  double initial_guess = 0;         // an initial test depth
  double new_x;                 // the new test depth
  double displace_x = 1e-6;     // how far you diplace the test depth
  double displace_scaling;      // the scaling after displacement
  double scaling_this_step;     // the scaling at the current step
  double fx;                    // function for Newton Raphson to find root
  double fx_displace;           // function after displacement
  double fx_derivative;         // derivative of the root finding function
  double x_change;              // the change in the scaling
  vector<double> F(4,0.0);      // holds the F values, cycles between nuclides
  double tolerance = 1e-7;      // the tolerance over which the scaling can change

  // go through the nuclides, selecting each based on the booleans fed to 
  // the routine  
  // first 10Be  
  if(nuclides_for_scaling[0])
  {
    // replace the F values
    F[0] =  F_10Be[0];
    F[1] =  F_10Be[1];
    F[2] =  F_10Be[2];
    F[3] =  F_10Be[3];
    new_x = initial_guess;
  
    int iterations = 0;
  
    do
    {
      iterations++;
      // get the scaling this step
      scaling_this_step =  exp(-new_x/Gamma[0])*F[0]+
                           exp(-new_x/Gamma[1])*F[1]+
                           exp(-new_x/Gamma[2])*F[2]+
                           exp(-new_x/Gamma[3])*F[3];
                           
      // create the function for root finding
      fx = scaling_this_step-single_scaling;
      
      // now displace the test
      displace_scaling =  exp(-(new_x+displace_x)/Gamma[0])*F[0]+
                           exp(-(new_x+displace_x)/Gamma[1])*F[1]+
                           exp(-(new_x+displace_x)/Gamma[2])*F[2]+
                           exp(-(new_x+displace_x)/Gamma[3])*F[3];

      fx_displace =  displace_scaling-single_scaling;
    
      fx_derivative = (fx_displace-fx)/displace_x;
      
      if(fx_derivative != 0)
      {
        new_x = new_x-fx/fx_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        x_change = fx/fx_derivative;
        //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
      }
      else
      {
        x_change = 0;
      }
  
    } while(fabs(x_change) > tolerance);      
    
//    cout << "========================================================" << endl;
//    cout << "TESTING SCALING IN LSDCRNP" << endl;
//    cout << "LINE 1742, scaling is: " << new_x << endl;
//    cout << "Iterations were: " << iterations << endl;
//    cout << "========================================================" << endl;
//    
    // now reset the F_values
    F_10Be[0] = exp(-new_x/Gamma[0])*F_10Be[0];
    F_10Be[1] = exp(-new_x/Gamma[1])*F_10Be[1];
    F_10Be[2] = exp(-new_x/Gamma[2])*F_10Be[2];
    F_10Be[3] = exp(-new_x/Gamma[3])*F_10Be[3];  
//    cout << "new scaling is:" << F_10Be[0] << endl;
  }

  // now 26Al
  if(nuclides_for_scaling[1])
  {
    // replace the F values
    F[0] =  F_26Al[0];
    F[1] =  F_26Al[1];
    F[2] =  F_26Al[2];
    F[3] =  F_26Al[3];
    new_x = initial_guess;
  
    //int iterations = 0;
  
    do
    {
      //iterations++;
      // get the scaling this step
      scaling_this_step =  exp(-new_x/Gamma[0])*F[0]+
                           exp(-new_x/Gamma[1])*F[1]+
                           exp(-new_x/Gamma[2])*F[2]+
                           exp(-new_x/Gamma[3])*F[3];
                           
      // create the function for root finding
      fx = scaling_this_step-single_scaling;
      
      // now displace the test
      displace_scaling =  exp(-(new_x+displace_x)/Gamma[0])*F[0]+
                           exp(-(new_x+displace_x)/Gamma[1])*F[1]+
                           exp(-(new_x+displace_x)/Gamma[2])*F[2]+
                           exp(-(new_x+displace_x)/Gamma[3])*F[3];

      fx_displace =  displace_scaling-single_scaling;
    
      fx_derivative = (fx_displace-fx)/displace_x;
      
      if(fx_derivative != 0)
      {
        new_x = new_x-fx/fx_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        x_change = fx/fx_derivative;
        //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
      }
      else
      {
        x_change = 0;
      }
  
    } while(fabs(x_change) > tolerance);      
    
    //cout << "========================================================" << endl;
    //cout << "TESTING SCALING IN LSDCRNP" << endl;
    //cout << "LINE 1742, scaling is: " << new_x << endl;
    //cout << "Iterations were: " << iterations << endl;
    //cout << "========================================================" << endl;
    
    // now reset the F_values
    F_26Al[0] = exp(-new_x/Gamma[0])*F_26Al[0];
    F_26Al[1] = exp(-new_x/Gamma[1])*F_26Al[1];
    F_26Al[2] = exp(-new_x/Gamma[2])*F_26Al[2];
    F_26Al[3] = exp(-new_x/Gamma[3])*F_26Al[3];  
  }

  // now 36Cl
  if(nuclides_for_scaling[2])
  {
    // replace the F values
    F[0] =  F_36Cl[0];
    F[1] =  F_36Cl[1];
    F[2] =  F_36Cl[2];
    F[3] =  F_36Cl[3];
    new_x = initial_guess;
  
    //int iterations = 0;
  
    do
    {
      //iterations++;
      // get the scaling this step
      scaling_this_step =  exp(-new_x/Gamma[0])*F[0]+
                           exp(-new_x/Gamma[1])*F[1]+
                           exp(-new_x/Gamma[2])*F[2]+
                           exp(-new_x/Gamma[3])*F[3];
                           
      // create the function for root finding
      fx = scaling_this_step-single_scaling;
      
      // now displace the test
      displace_scaling =  exp(-(new_x+displace_x)/Gamma[0])*F[0]+
                           exp(-(new_x+displace_x)/Gamma[1])*F[1]+
                           exp(-(new_x+displace_x)/Gamma[2])*F[2]+
                           exp(-(new_x+displace_x)/Gamma[3])*F[3];

      fx_displace =  displace_scaling-single_scaling;
    
      fx_derivative = (fx_displace-fx)/displace_x;
      
      if(fx_derivative != 0)
      {
        new_x = new_x-fx/fx_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        x_change = fx/fx_derivative;
        //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
      }
      else
      {
        x_change = 0;
      }
  
    } while(fabs(x_change) > tolerance);      
    
    //cout << "========================================================" << endl;
    //cout << "TESTING SCALING IN LSDCRNP" << endl;
    //cout << "LINE 1742, scaling is: " << new_x << endl;
    //cout << "Iterations were: " << iterations << endl;
    //cout << "========================================================" << endl;
    
    // now reset the F_values
    F_36Cl[0] = exp(-new_x/Gamma[0])*F_36Cl[0];
    F_36Cl[1] = exp(-new_x/Gamma[1])*F_36Cl[1];
    F_36Cl[2] = exp(-new_x/Gamma[2])*F_36Cl[2];
    F_36Cl[3] = exp(-new_x/Gamma[3])*F_36Cl[3];  
  }

  // now 14C
  if(nuclides_for_scaling[3])
  {
    // replace the F values
    F[0] =  F_14C[0];
    F[1] =  F_14C[1];
    F[2] =  F_14C[2];
    F[3] =  F_14C[3];
    new_x = initial_guess;
  
    //int iterations = 0;
  
    do
    {
      //iterations++;
      // get the scaling this step
      scaling_this_step =  exp(-new_x/Gamma[0])*F[0]+
                           exp(-new_x/Gamma[1])*F[1]+
                           exp(-new_x/Gamma[2])*F[2]+
                           exp(-new_x/Gamma[3])*F[3];
                           
      // create the function for root finding
      fx = scaling_this_step-single_scaling;
      
      // now displace the test
      displace_scaling =  exp(-(new_x+displace_x)/Gamma[0])*F[0]+
                           exp(-(new_x+displace_x)/Gamma[1])*F[1]+
                           exp(-(new_x+displace_x)/Gamma[2])*F[2]+
                           exp(-(new_x+displace_x)/Gamma[3])*F[3];

      fx_displace =  displace_scaling-single_scaling;
    
      fx_derivative = (fx_displace-fx)/displace_x;
      
      if(fx_derivative != 0)
      {
        new_x = new_x-fx/fx_derivative;
      
        // check to see if the difference in erosion rates meet a tolerance
        x_change = fx/fx_derivative;
        //cout << "Change is: " << eff_e_change << " and erosion rate is: " << eff_e_new << endl;
      }
      else
      {
        x_change = 0;
      }
  
    } while(fabs(x_change) > tolerance);      
    
    //cout << "========================================================" << endl;
    //cout << "TESTING SCALING IN LSDCRNP" << endl;
    //cout << "LINE 1742, scaling is: " << new_x << endl;
    //cout << "Iterations were: " << iterations << endl;
    //cout << "========================================================" << endl;
    
    // now reset the F_values
    F_14C[0] = exp(-new_x/Gamma[0])*F_14C[0];
    F_14C[1] = exp(-new_x/Gamma[1])*F_14C[1];
    F_14C[2] = exp(-new_x/Gamma[2])*F_14C[2];
    F_14C[3] = exp(-new_x/Gamma[3])*F_14C[3];  
  }


}
