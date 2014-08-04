#include<iostream>
#include<vector>
#include<string>
using namespace std;

#ifndef VolumeParticleInfo_HPP
#define VolumeParticleInfo_HPP

class VolumeParticleInfo
{
 public:
    VolumeParticleInfo()				{ create(); }
    VolumeParticleInfo( const char* fname)           { create(fname); }

    double return_mass_fraction(int Type_index, int GSDindex);
    double return_surface_area(int Type_index, int GSDindex, double mass);

    int get_n_types()					{ return n_types; }
    int get_n_sizes()					{ return n_sizes; }
    double get_ParticleTargetMass()		{ return ParticleTargetMass; }
    string get_type_name(int index);
    double get_type_density(int index);

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
