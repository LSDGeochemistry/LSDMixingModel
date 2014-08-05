#include <iostream>
#include "../CRUNCH_engine.hpp"
using namespace std;

int main ()
{
	string crunch_pname = "/home/smudd/SMMDataStore/CRUNCH_binary/";
	string run_pname = "/home/smudd/SMMDataStore/devel_projects/MixingModel/trunk/test_run";
  string master_fname = "master_crunch.in";

	CRUNCH_engine Ceng(crunch_pname, run_pname, master_fname);
	//Ceng.print_master();

	double new_tspacing = 20;
	string value_name = "spatial_profile";
	Ceng.modify_CRUNCH_value(value_name, new_tspacing);

	int n_ts = 1;
	int n_conditions;
	vector<double> pH_values;
	vector<double> spacings;
	vector<double> top_depths;
	vector<double> bottom_depths;
 	list< vector<double> > concentrations;
 	list< vector<double> > mineral_vfracs;
 	list< vector<double> > mineral_bsa;
 	list< vector<double> > reaction_rates;
	Ceng.parse_CRUNCH_files(n_ts, n_conditions, pH_values, spacings,
							top_depths, bottom_depths,
						    concentrations, mineral_vfracs,
							mineral_bsa, reaction_rates);

	cout << "the number of conditions is: " << n_conditions << endl;

	Ceng.set_spacing_line_from_depths(top_depths, bottom_depths);

	Ceng.create_CRUNCH_in_file(n_conditions, pH_values,
							top_depths, bottom_depths,
						    concentrations, mineral_vfracs,
							mineral_bsa);

	Ceng.parse_parent_material_file();
	//Ceng.call_CRUNCH();

	//Ceng.print_master();
}
