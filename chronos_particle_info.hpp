//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// chronos_particle_info
// An object that holds properties of particles
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
#include<vector>
#include<string>
#include "LSDParticle.hpp"
using namespace std;

#ifndef Particle_info_H
#define Particle_info_H

/// Particle_info class. Holds information about particles in general so that 
///  one does not need to recalculate things for every particle individually. 
class Particle_info
{
 public:
    
    /// @brief Default constructor. Throws and error since you need io initiate object with a file
    /// @author SMM
    /// @date 01/01/2008
    Particle_info()				{ create(); }
    
    /// @brief The constructior that takes a particle information file which does things like defines the
    ///  types, enters mass fractions of minerals, etc
    /// @param fname The name of the partle info file
    /// @author SMM
    /// @date 01/01/2008
    Particle_info( const char* fname)           { create(fname); }

    /// @brief This prints the particle information to screen for bug checking
    /// @author SMM
    /// @date 01/01/2008
    void details_to_screen();
    
    /// @brief The creates time vectors: Used to create times at which weathering extents are calculated
    ///  The wethering extents are calculated along a vector of times, which then can be used for each particle
    ///  this is only used for the time based weathering (i.e., not using CRUNCH).
    ///  It ALSO calculates all the mass fractions through time. Calling this function can take some time
    ///   because it sets in motion calculation of a numer of power law equations over many timesteps for many
    ///   minerals
    /// @param dt The time increment
    /// @param n_timesteps The number of time increments in the time vector
    /// @author SMM
    /// @date 01/01/2008
    void generate_time_vecs(double dt,int n_timesteps);
    
    /// @brief Prints the time vectors to file. Used for debugging
    /// @param fname the name of the file
    /// @author SMM
    /// @date 01/01/2008
    void print_time_vecs(const char* fname);
    
    /// @brief This calculates particle fractions (the fraction of the volume of particles) 
    ///  based on the mass fraction of particles
    /// @param chi_frac These are the mass fractions
    /// @return The volme fractions of the particles
    /// @author SMM
    /// @date 01/01/2008
    vector<double> get_particle_fractions(vector<double> chi_frac);
    
    /// @brief Returns the fraction of the height. Not quite sure what this is (SMM, 23/07/2018)
    /// @author SMM
    /// @date 01/01/2008
    double retrieve_h_frac(int Timestep, int type);
    
    /// @brief This gets the mass rafction remaining of type designated for a given timestep
    ///  based on weathering calculations done in the generate_time_vecs function
    /// @param Timestep The integer timestep. Needs to have the same dt as in generate_time_vecs
    /// @param type The type of the particle
    /// @author SMM
    /// @date 01/01/2008
    vector<double> retrieve_mass_fracs(int Timestep, int type);
    
    vector<double> calculate_m_0_vec(double dx, double h_0, double phi);

    vector<double> calculate_h_congruent(LSDParticle& tpart, double deltat, double phi);
    
    /// @brief This calculates the weathering extent of a particle at a given time on demand
    ///  Useful if you just want to simulate the evolution of a single particle, but not useful for
    ///  multi-paticle simulations since computationally expensive. Use generate_time_vecs instead
    /// @author SMM
    /// @date 01/01/2008
    vector<double> calculate_weathering_on_demand(double t_ime,int type);
    
    /// @brief Get the density of the particle types
    /// @return A vector of the densities
    /// @author SMM
    /// @date 01/01/08
    vector<double> get_rho_p_vec()			{ return rho_p_vec; }
    
    /// @brief Gets the vector of the type indices
    ///  This is a slightly stupid way of doing things but it was before I started using dicts for such purposes
    /// @return The type indices in a vector
    /// @author SMM
    /// @date 01/01/2008
    vector<int> get_type_index()			{ return type_index; }

 private:
    int n_types;
    
    /// An index into the type names.
    vector<int> type_index;
    
    /// The names of the types
    vector<string> type_name;   
    
    /// Vector of the a coefficient of weathering for each type
    vector<double> a_vec;
    
    /// Vector of the alpha coefficient for weathering for each type
    vector<double> alpha_vec;
    
    /// Vector of the b coefficient for weathering of each type
    vector<double> b_vec;
    
    /// Vector of the beta coefficient for weathering of each type
    vector<double> beta_vec;
    
    /// Vector of the density of the primary mineral of each type
    vector<double> rho_p_vec;
    
    /// Vector of the clay mineral product density wor each type
    vector<double> rho_c_vec;
    
    /// The molar wieght of the primary mineral of each type
    vector<double> w_p_vec;
    
    /// The molar weight of the clay mineral product of each type
    vector<double> w_c_vec;
    
    /// The stoichiometric coefficient of the primary mineral of each type (for the weathering reaction)
    vector<double> stoic_p_vec;
    
    /// The stoichiometric coefficient of the clay mineral of each type (for the weathering reaction)
    vector<double> stoic_c_vec;
    
    
    vector<double> D_vec;
    vector<double> chi_vec;
    
    /// The stoichiometric coefficient of SiO2 realeased by the weathering reaction
    vector<double> SiO2_stoic;
    
    /// The stoichiometric coefficient of SiO2 realeased by the weathering reaction
    vector<double> Na_stoic;
    
    /// The stoichiometric coefficient of SiO2 realeased by the weathering reaction
    vector<double> Ca_stoic;
    
    /// The stoichiometric coefficient of SiO2 realeased by the weathering reaction
    vector<double> K_stoic;
    
    /// The stoichiometric coefficient of SiO2 realeased by the weathering reaction
    vector<double> Mg_stoic;

    int n_timesteps;
    double dt;
    vector< vector<double> > m_p_frac_vecs;		// vector of vectors describing
    							// the fraction of mass remaining in
    							// primary mineral
    vector< vector<double> > m_c_frac_vecs;		// vector of vectors describing
    							// the fraction of mass produced in
    							// clay mineral
    vector< vector<double> > h_frac_vecs;		// vector of vectors describing
    							// the proportion of original height
    							// of the particle
    vector< vector<double> > m_rate_frac_vecs;		// vector of vectors describing
    							// the rate as a fraction of the original
    							// mass of the particle in units yr^-1

    void create();
    void create(const char*);

};

#endif

