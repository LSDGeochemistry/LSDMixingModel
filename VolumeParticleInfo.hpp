//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// VolumeParticleInfo
// An object for keeping track of particle information for geochemical and sediment
// tracing experiments
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group mixing model
//  for exploring hillslope mixing and particle weathering
//
// Developed by:
//  Simon M. Mudd
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


#include<iostream>
#include<fstream>
#include<vector>
#include<string>
using namespace std;

#ifndef VolumeParticleInfo_HPP
#define VolumeParticleInfo_HPP

/// This class holds information about particles so that computations don't
/// need to be repeated for each particle
class VolumeParticleInfo
{
 public:
    
    /// @brief Default create function
    /// @author SMM
    /// @date 01/01/11
    VolumeParticleInfo()				{ create(); }
    
    /// @brief Create function that starts the particle from a particle file
    /// @param The name of the parameter file (with extension) 
    /// @author SMM
    /// @date 01/01/11    
    VolumeParticleInfo( const char* fname)           { create(fname); }
    
    /// @brief Returns the mass fraction of the particle
    /// @param Type_index the particle type (an integer. Could be, for example an integrer for a mineral type) 
    /// @param GSDindex The index for the grain size (instead of solving for an double grain size)
    /// @returns The mass fraction (as a double)
    /// @author SMM
    /// @date 01/01/11   
    double return_mass_fraction(int Type_index, int GSDindex);
    
    /// @brief Returns the surface area of the particle
    /// @param Type_index the particle type (an integer. Could be, for example an integrer for a mineral type) 
    /// @param GSDindex The index for the grain size (instead of solving for an double grain size)
    /// @param mass The mass of the particle
    /// @returns The surface area (as a double)
    /// @author SMM
    /// @date 01/01/11      
    double return_surface_area(int Type_index, int GSDindex, double mass);

    // getter functions
    int get_n_types()					{ return n_types; }
    int get_n_sizes()					{ return n_sizes; }
    double get_ParticleTargetMass()		{ return ParticleTargetMass; }
    string get_type_name(int index);
    
    double get_type_density(int index);
    vector<double> get_type_mfracs() { return type_mfracs; }
    vector<int> get_type_index()   { return type_index; }


 private:
    int n_types;
    int n_sizes;
    vector<int> type_index;
    vector<string> type_name;
    vector<double> size_index;			// you get particle size by taking 2^(size index)
    vector<double> type_mfracs;
    vector<double> size_mfracs;
    vector<double> densities;			// in kg/m^3
    vector<double> lambdas;				// roughness, used in surface area calucaltion
    double ParticleTargetMass;			// the target mass of individual particles when particles are
                       // created the inserter tries to get as close to this mass whilst still
                       // conserving mass

	double parent_fracs[20][20];		// this stores the parent fraction array. It
                       	// is designed to hold 20 types and 20 particle
                       	// size classes. This can be expanded later
	double surface_area_multiplier[20][20];
                       // this converts mass to surface area for a given
                       // mass needs to be given in kg and area will be
                       // returned in m^2

    void create();
    void create(const char*);

};

#endif
