#include <iostream>
#include "../CRUNCH_engine.hpp"
using namespace std;

int main ()
{
	string crunch_pname = "/cygdrive/m/papers/mixing_model_2014/source/CRUNCH_binary/";
	string run_pname = "/cygdrive/m/papers/mixing_model_2014/source/runs/run1/";
  string master_fname = "master_crunch.in";

  //cout << "Crunch: " << crunch_pname << endl;
  //cout << "Run: " << run_pname << endl;
  //cout << "Master: " << master_fname << endl;

	CRUNCH_engine Ceng(crunch_pname, run_pname, master_fname);

	//Ceng.print_master();

	double new_tspacing = 20;
	string value_name = "spatial_profile";
	Ceng.modify_CRUNCH_value(value_name, new_tspacing);

	int n_ts = 1;
	vector<double> spacings;
	vector<double> top_depths;
	vector<double> bottom_depths;
	
	// set up some top and bottom depths
	double ds = 0.1;       // depth spacing
	double td = 0;         // top depth
	double bd = td+0.1;    // bottom depth
	
	// first section of depths
	top_depths.push_back(td);
	bottom_depths.push_back(bd);
	for(int i = 0; i<3; i++)
	{
    td = bd;
    bd = td+ds;
    
	  top_depths.push_back(td);
	  bottom_depths.push_back(bd); 
  }

  // second section of depths
  ds = 0.2;
	for(int i = 0; i<4; i++)
	{
    td = bd;
    bd = td+ds;
    
	  top_depths.push_back(td);
	  bottom_depths.push_back(bd); 
  }
  
  Ceng.set_spacing_line_from_depths(top_depths,bottom_depths);
	
	int n_conditions;
	vector<double> pH_values;

 	list< vector<double> > concentrations;
 	list< vector<double> > mineral_vfracs;
 	list< vector<double> > mineral_bsa;
 	list< vector<double> > reaction_rates;
	Ceng.parse_CRUNCH_files(n_ts, n_conditions, pH_values, spacings,
							top_depths, bottom_depths,
						    concentrations, mineral_vfracs,
							mineral_bsa, reaction_rates);
  cout << "Parsed some crunch files" << endl;


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
