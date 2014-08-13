#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<math.h>
#include "VolumeParticleInfo.hpp"
using namespace std;

#ifndef VolumeParticleInfo_CPP
#define VolumeParticleInfo_CPP

void VolumeParticleInfo::create()
{
	cout << "VolumeParticleInfo line 14 error, you need to have a filename\n";
}

// this function loads the parameters from the file with the name fname
// file format
// n_types
// type_index type_name mass_fraction density lambda
// ...
//
// n_sizes
// type_index log_2_size mfrac
void VolumeParticleInfo::create( const char* fname )
{
	ifstream in;
	in.open(fname);

	vector<int> t_type_index;
	int tti;
	vector<string> t_type_name;
	string ttn;
	vector<double> t_size_index;
	double tsi;

	vector<double> t_type_mfrac;
	double ttmf;
	vector<double> t_size_mfrac;
	double tsmf;

	vector<double> t_densities;
	double trho;
	vector<double> t_lambdas;
	double tl;

	int t_ntypes;
	int t_nsizes;
	double t_ParticleTargetMass;

	in >> t_ntypes >> t_ParticleTargetMass;
	ParticleTargetMass = t_ParticleTargetMass;

	for (int i = 0; i< t_ntypes; i++)
	{
		in >> tti >> ttn >> ttmf >> trho >> tl;
		t_type_index.push_back(tti);
		t_type_name.push_back(ttn);
		t_type_mfrac.push_back(ttmf);
		t_densities.push_back(trho);
		t_lambdas.push_back(tl);
	}

	// input the size classes
	// the file format is
	// 2 2 0.1
	// 3 2.33 0.2
	// 4 2.66 0.1
	// where the first column is the index
	// where the particle size is 2^(second column) and the second column is
	// the mass fraction of that particle size.

	// get the number of particle sizes
	in >> t_nsizes;
	for (int i = 0; i<t_nsizes; i++)
	{
		t_size_index.push_back(0.0);
		t_size_mfrac.push_back(0.0);
	}

	int tGSDi;
	for (int i = 0; i< t_nsizes; i++)
	{
		in >> tGSDi >> tsi >> tsmf;
		t_size_index[tGSDi] = tsi;
		t_size_mfrac[tGSDi] = tsmf;
	}

	// set the private data
  n_types = t_ntypes;
  n_sizes = t_nsizes;
  type_index = t_type_index;
  type_name = t_type_name;
	size_index = t_size_index;
	type_mfracs = t_type_mfrac;
	size_mfracs = t_size_mfrac;
	densities = t_densities;
	lambdas = t_lambdas;

	// set the parent fraction and multiplier array to zero
	for (int i = 0; i<20; i++)
	{
		for (int j = 0; j<20; j++)
		{
			parent_fracs[i][j]= 0;
			surface_area_multiplier[i][j] = 0;
		}
	}

	// now get the fractions and the surface area multipliers
	for (int type = 0; type<n_types; type++)
	{
		for (int sizet = 0; sizet<n_sizes; sizet++)
		{
			parent_fracs[type][sizet] = type_mfracs[type]*size_mfracs[sizet];
			surface_area_multiplier[type][sizet] = 6*lambdas[type]/(0.001*pow(2.0,t_size_index[sizet])*densities[type]);
						// the factor of 0.001 is because the sizes are reported in log base 2 of particle diamter in mm
						// but the surface area is reported in metres^2
		}
	}

	in.close();

}

string VolumeParticleInfo::get_type_name(int index)
{
	if (index > n_types)
	{
		cout << "error, VolumeParticleInfo line 126; index exceeds number of types";
	}
	else
	{
		return type_name[index];
	}
}

double VolumeParticleInfo::get_type_density(int index)
{
	if (index > n_types)
	{
		cout << "error, VolumeParticleInfo line 139; index exceeds number of types";
	}
	else
	{
		return densities[index];
	}
}

double VolumeParticleInfo::return_mass_fraction(int Type_index, int GSDindex)
{
	return parent_fracs[Type_index][GSDindex];
}

double VolumeParticleInfo::return_surface_area(int Type_index, int GSDindex, double mass)
{
	return surface_area_multiplier[Type_index][GSDindex]*mass;
}

#endif
