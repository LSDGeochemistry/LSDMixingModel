//dParticle.h
// header file for the discrete particle object
// this version of the dicrete particle keeps track
// of two integers, the time and the position

#include<iostream>
#include<vector>
#include<string>
#include "tParticle.hpp"
using namespace std;

#ifndef Particle_info_H
#define Particle_info_H

class Particle_info
{
 public:
    Particle_info()				{ create(); }
    Particle_info( const char* fname)           { create(fname); }

    void details_to_screen();
    void generate_time_vecs(double dt,int n_timesteps);
    void print_time_vecs(const char* fname);
    vector<double> get_particle_fractions(vector<double> chi_frac);
    double retrieve_h_frac(int Timestep, int type);
    vector<double> retrieve_mass_fracs(int Timestep, int type);
    vector<double> calculate_m_0_vec(double dx, double h_0, double phi);

    vector<double> calculate_h_congruent(tParticle& tpart, double deltat, double phi);
    vector<double> calculate_weathering_on_demand(double t_ime,int type);

    vector<double> get_rho_p_vec()			{ return rho_p_vec; }
    vector<int> get_type_index()			{ return type_index; }

 private:
    int n_types;
    vector<int> type_index;
    vector<string> type_name;
    vector<double> a_vec;
    vector<double> alpha_vec;
    vector<double> b_vec;
    vector<double> beta_vec;
    vector<double> rho_p_vec;
    vector<double> rho_c_vec;
    vector<double> w_p_vec;
    vector<double> w_c_vec;
    vector<double> stoic_p_vec;
    vector<double> stoic_c_vec;
    vector<double> D_vec;
    vector<double> chi_vec;
    vector<double> SiO2_stoic;
    vector<double> Na_stoic;
    vector<double> Ca_stoic;
    vector<double> K_stoic;
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

