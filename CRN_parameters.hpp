//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// LSDCRNParameters.hpp
//
// Land Surface Dynamics Cosmogenic Radionuclide Parameters Object
//
// This keeps track of paramters used to calculate the evolution of 
// in situ cosmogenic nuclides. 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for calculating concntration of environmental tracers, CRNs, TCN, fallout
//  nuclides
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
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
// either version 2 of the License, or (at your option) any later version.
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
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include "TNT/tnt.h"
#include "LSDStatsTools.hpp"
using namespace std;
using namespace TNT;

#ifndef CRN_parameters_H
#define CRN_parameters_H

/// @brief This class contains parameters used in cosmogenic nuclide calculations
/// It sits seperately from the particle object since it applies to an
/// entire environment and not just an individual particle. 
/// Seperating the object in this way reduces memory redundancy 
class CRN_parameters
{
  public:
    /// @brief The default constructor. It is the only possible constructor
    LSDCRNParameters()         { create(); }
  
    /// This is a friend class so that it can be called from the particle 
    friend class LSDCRNParticle;

    /// @brief This resets the F, Gamma and P0 values so that they conform to 
    /// Granger and Smith 2000 scaling. Adopted from from Vermeesh 2007
    /// @author SMM
    /// @date 01/01/2010
    void set_Granger_parameters();
  
    /// @brief This resets the F, Gamma and P0 values so that they conform to 
    /// Schaller (2009) scaling. Adopted from from Vermeesh 2007
    /// @author SMM
    /// @date 01/01/2010
    void set_Schaller_parameters();
  
	//void scale_F_values(double single_scaling);	// parameter values
	void update_10Be_decay(double new_decay)	{ lambda_10Be = new_decay; }
	void update_10Be_P0(double new_P0)			{ P0_10Be = new_P0; }
    
    void scale_F_values(vector<bool> nuclides_for_scaling);
    
    double NCEPatm_2(double lat, double lon, double site_elev); 
    
    /// @brief function for loading parameters that allow pressure calculation
    /// from elevation
    /// @author SMM
    /// @date 02/12/2014
    void load_parameters_for_atmospheric_scaling(string path_to_params);
    
    double stone2000sp(double lat, double Fsp);

  private:
    /// @brief This is called by the default constructor. 
    /// It is the only possible constructor
    void create();

    /// The attenudation lengths in g/cm^2. Each element in the array
    /// refers to a different production mechanism
    double Gamma[4];			// attenuation legths in g/cm^2

    /// F values allocating production to spallation and muon production for 10Be
    double F_10Be[4];
  
    /// F values allocating production to spallation and muon production for 14C
    double F_14C[4];
  
    /// F values allocating production to spallation and muon production for 26Al
    double F_26Al[4];
  
    /// F values allocating production to spallation and muon production for 36Cl
    double F_36Cl[4];
  
    /// F values allocating production to spallation and muon production for 21Ne
    double F_21Ne[4];
  
    /// F values allocating production to spallation and muon production for 3He
    double F_3He[4];

    /// topographic shielding
    /// other shielding and scaling calucalted using the scale_F_values function
    double S_t;
    

    /// decay rate for 10Be in yr-1
    double lambda_10Be;
    
    /// decay rate for 26Al in yr-1	
    double lambda_26Al;
    
    /// decay rate for 14C in yr-1	
    double lambda_14C;	
    
    /// decay rate for 36Cl in yr-1	
    double lambda_36Cl;		
  
    /// production rate for 10Be in a/g/yr
    double P0_10Be;	
    
    /// production rate for 26Al in a/g/yr		
    double P0_26Al;		
    
    /// production rate for 14C in a/g/yr
    double P0_14C;			
    
    /// production rate for 36Cl in a/g/yr
    double P0_36Cl;		
    
    /// production rate for 21Ne in a/g/yr
    double P0_21Ne;			
    
    /// production rate for 3He in a/g/yr
    double P0_3He;			
    
    /// latitude
    double lat;
    
    /// longitude
    double lon;
    
    /// The elevation of the site
    double site_elev;
    double Fsp;
    
    /// levels: the levels for the atmospheric scaling of pressure
    vector<double> levels;
  
    /// This is an index for the latitudes for atmospheric scaling
    vector<double> NCEPlat;
  
    /// This is an index for the longitudes for atmospheric scaling
    vector<double> NCEPlon;
  
    /// This is an array holding sea level perssures
    Array2D<double> meanslp;
  
    /// This is an array healing mean temperatures
    Array2D<double> meant1000;
  
    /// This is a vector of arrays holding something called gp_hgt;
    vector< Array2D<double> > gp_hgt;
  
    /// This is a vector of arrays holding something called gp_hgt;
    vector< Array2D<double> > gm_hgt;  
  

    
};

#endif
