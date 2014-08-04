#include <iostream>
#include <string>
#include "VolumeParticleInfo.h"
using namespace std;

int main()
{
	VolumeParticleInfo vpi("VolumeParticleData.in");

	int n_types = vpi.get_n_types();
	int n_sizes = vpi.get_n_sizes();

	cout << "n_types: " << n_types << " and n_sizes: " << n_sizes << endl;

	for(int i = 0; i<n_types; i++)
	{
		cout << "type_index: " << i << " and type: " << vpi.get_type_name(i) << endl;
	}

	double sum =0;
	// get the mass fraction matrix
	for (int type = 0; type<n_types; type++)
	{
		for (int sizet =0; sizet<n_sizes; sizet++)
		{
			cout << "type: " << type << " size: " << sizet
			     << " mfrac: " << vpi.return_mass_fraction(type, sizet)
			     << endl;
			sum +=vpi.return_mass_fraction(type, sizet);
		}
	}


	// get some sample masses
	cout << endl << endl;
	double mass = 1;
	for (int type = 0; type<n_types; type++)
	{
		for (int sizet =0; sizet<n_sizes; sizet++)
		{
			cout << "type: " << type << " size: " << sizet
			     << " surface area: " << vpi.return_surface_area(type, sizet,mass)
			     << endl;
		}
	}

	cout << "total mass fraction is: " << sum << endl;
}
