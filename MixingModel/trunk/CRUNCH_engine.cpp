#ifndef CRUNCH_engine_CPP
#define CRUNCH_engine_CPP

#include <iostream>
#include <list>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "CRUNCH_engine.hpp"
#include "mathutil.hpp"
using namespace std;

//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_engine::create()
{
	cout << "you need to initialize with a master file!" << endl;
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function creates a CRUNCH_engine object by reading a master file
// the crunch infile wil be created based on this infile
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_engine::create(string crunch_path, string run_path, string master_filename)
{
  // first make sure the pathname has a slash
  cout << "Crunch path is: " << crunch_path << endl;
  string lchar = crunch_path.substr(crunch_path.length()-1,1);
  string slash = "/";
  cout << "lchar is " << lchar << " and slash is " << slash << endl;
    
  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl;  
    crunch_path = crunch_path+slash;
    cout << "The CRUNCH pathname is: " << crunch_path << endl;
  } 
  CRUNCH_path = crunch_path;
  
  // first make sure the run pathname has a slash
  lchar = run_path.substr(run_path.length()-1,1);
  cout << "run path is: " << run_path << endl;  
  if (lchar != slash)
  {
    cout << "You forgot the frontslash at the end of the path. Appending." << endl;  
    run_path = run_path+slash;
    cout << "The RUN pathname is: " << run_path << endl;
  } 
  RUN_path = run_path;  
  
  cout << "And master fname is: " << master_filename << endl;

	ifstream master_in;
	string this_master = RUN_path+master_filename;
	master_in.open(this_master.c_str());
	string stemp;
	char temp[500];

	list<string> master_read;
	while(master_in.getline(temp,500))
	{
		stemp = temp;
		master_read.push_back(stemp);
	}
	master_in.close();

	master_infile = master_read;

	get_minerals_and_species();
	get_constant_conditions();
	//parse_parent_material_file();
	get_mineral_properties();
  cout << "Finished constructing crunch engine object." << endl;
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this gets a few constant parameters that will go into every condition statement
void CRUNCH_engine::get_constant_conditions()
{
	// first thing: get parameters from master file
	list<string>::iterator l_iter;
	list<string>::iterator condition_iter;
	l_iter = master_infile.begin();
	size_t found = string::npos; 			// integer for testing if the command block has been found
	string solute_condition_name = "Condition SoluteWater";
	// advance through file until you find the solutewater data block
	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(solute_condition_name);
		condition_iter = l_iter;
		l_iter++;
	}
	found = 23;					// reset the found value (23 is arbitrary)

	// now loop through the condition block repeatedly, gathering the correct data
	string end_name = "END";
	string temp_name = "temperature";
	// first get the temperature
	while(l_iter != master_infile.end() && (*l_iter).find(end_name)== string::npos)
	{
		if( (*l_iter).find(temp_name)!=string::npos)
		{
			temperature_line = (*l_iter);
		}
		l_iter++;
	}
	// reset the iterator to the beginning
	l_iter = condition_iter;
	string solid_dens_name = "SolidDensity";
	while(l_iter != master_infile.end() && (*l_iter).find(end_name)== string::npos)
	{
		if( (*l_iter).find(solid_dens_name)!=string::npos)
		{
			density_line = (*l_iter);
		}
		l_iter++;
	}

	// now get the exchanges
	l_iter = condition_iter;
	string ex_name = "X";
	while(l_iter != master_infile.end() && (*l_iter).find(end_name)== string::npos)
	{
		if( (*l_iter).find(ex_name)!=string::npos)
		{
			exchange_lines.push_back(*l_iter);
		}
		l_iter++;
	}

	cout << "temp line is: " << temperature_line << endl;
	cout << "density line is: " << density_line << endl;
	int n_exchange = exchange_lines.size();
	for (int i = 0; i<n_exchange; i++)
	{
		cout << "exchange["<<i<<"]: " << exchange_lines[i] << endl;
	}
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// gets the minerals and species from the master file
void CRUNCH_engine::get_minerals_and_species()
{
	// this function goes through all the crunch files parsing the data so that it
	// can be reinserted into a restarted run
	//
	// n_ts is the timestep number (usually 1)
	//

	// the first step is to caluclate the number of primary species
	// this requires searching through the master file
	// search through list looking for the correct syntax
	list<string>::iterator l_iter;
	l_iter = master_infile.begin();
	size_t found = string::npos; 	// integer for testing if the command block has been found
	string primary_name = "PRIMARY_SPECIES";

	list<string> p_species_names_temp;
	list<string> mineral_names_temp;
	//list<string> mineral_ssa_temp;

	// advance through file until you find the primary species
	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(primary_name);
		//if(found!=string::npos)
		//{
		//	cout << "Found primary species block"<< endl;
		//}
		l_iter++;
	}
	found = 23;					// reset the found value (23 is arbitrary)

	// now iterate through these until you reach the end
	string end_name = "END";
	while(l_iter != master_infile.end() && (*l_iter).find(end_name)== string::npos)
	{
		p_species_names_temp.push_back(*l_iter);
		//cout << "species: " << (*l_iter) << endl;
		l_iter++;
	}

	// now do the same for the mineral names
	l_iter = master_infile.begin();
	found = string::npos; 			// reset found
	string mineral_name = "MINERALS";
	// advance through file until you find the primary species
	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(mineral_name);
		//if(found!=string::npos)
		//{
		//	cout << "Found mineral block"<< endl;
		//}
		l_iter++;
	}
	found = 23;					// reset the found value (23 is arbitrary)

	string delim = " ";
	while(l_iter != master_infile.end() && (*l_iter).find(end_name)== string::npos)
	{
		vector<string> line_words;
		split_string((*l_iter), delim, line_words);
		mineral_names_temp.push_back(line_words[0]);
		//cout << "mineral name is: " << line_words[0] << endl;
		l_iter++;
	}

	// assign the primary species and the mineral species
	// no clue why direct assignment doesn't work but
	// this is the only way I could avoid a stack dump.
	l_iter = p_species_names_temp.begin();
	while(l_iter != p_species_names_temp.end())
	{
		p_species_names.push_back(*l_iter);
		l_iter++;
	}
	l_iter = mineral_names_temp.begin();
	while(l_iter != mineral_names_temp.end())
	{
		mineral_names.push_back(*l_iter);
		l_iter++;
	}

	// now print out the primary species and the mineral species
	//l_iter = p_species_names.begin();
	//cout << "primary species:" << endl;
	//while(l_iter != p_species_names.end())
	//{
	//	cout << (*l_iter) << endl;
	//	l_iter++;
	//}
	//l_iter = mineral_names.begin();
	//cout << "minerals:" << endl;
	//while(l_iter != mineral_names.end())
	//{
	//	cout << (*l_iter) << endl;
	//	l_iter++;
	//}

}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function takes a parent material file an absorbs the information
// at the moment this is not very flexible:
// the minerals have to be the same as those listed in the
// master infile and
// THEY HAVE TO APPEAR IN THE SAME ORDER
void CRUNCH_engine::parse_parent_material_file()
{
	// read in the material file
	ifstream material_in;
	string pmat_full_name = RUN_path+"parent_material.pmat";
	material_in.open(pmat_full_name.c_str());
	if (material_in.fail())
	{
		cout << "parent material file does not exist\n";
		exit(1);
	}

	vector<double> empty_vec;
	init_V_fracs = empty_vec;
	ssa = empty_vec;

	string temp;
	double starting_vfrac;
	double starting_ssa;

	while(material_in >> temp >> starting_vfrac >> starting_ssa)
	{
		init_V_fracs.push_back(starting_vfrac);
		ssa.push_back(starting_ssa);
	}
	material_in.close();
}


//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function allows the user to modify parameter values within
// the crunch master file
void CRUNCH_engine::modify_CRUNCH_value(string value_name, double value)
{
	// get the new name of the line
	string num = dtoa(value);
	string space = " ";
	string new_line = value_name+space+num;
	//cout << "new master line; " << new_line << endl;

	// search through list looking for the correct syntax
	list<string>::iterator l_iter;
	l_iter = master_infile.begin();
	size_t found = string::npos; 	// integer for testing if the command block has been found


	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(value_name);
		//cout << "found is: " << found << endl;
  		if (found!=string::npos)
  		{
    		//cout << "spatial found at: " << int(found) << endl;

			// now replace the line
			(*l_iter) = new_line;
		}
		l_iter++;
	}
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this  reads a series of files containing
// information about pore chemistry, mineral contents, mineral
// surface areas, etc at nodes in the profile and then prints
// these as conditions for the model
//
// the way the particle bins object gets information out is in depth
// intervals, each interval has a top and a bottom. CRUNCH, however, has
// cells where the location is at the centrepoint
//
// we will have both 'soil' zones and 'saprolite' zones. In fact these will be
// zones of pdz and caz.
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_engine::parse_CRUNCH_files(int n_ts, int& n_conditions,
            int this_bin, int n_bins, int n_cells_in_bin
						vector<double>& pH_values, vector<double>& spacings,
						vector<double>& top_depths, vector<double>& bottom_depths,
						list< vector<double> >& concentrations, list< vector<double> >& mineral_vpercents,
						list< vector<double> >& mineral_ssa, list< vector<double> >& reaction_rates)
{
	// empty vectors for resetting data vectors and lists
	vector<double> depth;
	vector<double> empty_vec;
	list< vector<double> > empty_list_vec;
	
	// the mineral vecs are treated differently from concentration vecs since they
	// are indexed by cell.
	// This will need to change later since it is inefficient but no time for that
	// now SMM 10/08/2014
	int total_cells = n_bins*n_cells_in_bin;
	vector<double> mineral_empty_vec(total_cells,0.0);

	pH_values = empty_vec;
	spacings = empty_vec;
	top_depths = empty_vec;
	bottom_depths = empty_vec;
	concentrations = empty_list_vec;
	mineral_vpercents = empty_list_vec;
	mineral_ssa = empty_list_vec;
	reaction_rates = empty_list_vec;

	// initialize the lists
	int n_species = p_species_names.size();
	for(int i = 0; i<n_species; i++)
	{
		concentrations.push_back(empty_vec);
	}
	int n_minerals = mineral_names.size();
	for(int i = 0; i<n_minerals; i++)
	{
		mineral_vpercents.push_back(mineral_empty_vec);
		mineral_ssa.push_back(mineral_empty_vec);
	}

	// first set up the files
	string num = itoa(n_ts);
	string ext = ".out";
	string area_fname = "area";
	area_fname  = RUN_path+area_fname+num+ext;
	ifstream area_in;
	area_in.open(area_fname.c_str());
	string conc_fname = "totcon";
	conc_fname  = RUN_path+conc_fname+num+ext;
	ifstream conc_in;
	conc_in.open(conc_fname.c_str());
	string pH_fname = "pH";
	pH_fname  = RUN_path+pH_fname+num+ext;
	ifstream pH_in;
	pH_in.open(pH_fname.c_str());
	string rate_fname = "rate";
	rate_fname  = RUN_path+rate_fname+num+ext;
	ifstream rate_in;
	rate_in.open(rate_fname.c_str());
	string volume_fname = "volume";
	volume_fname  = RUN_path+volume_fname+num+ext;
	ifstream volume_in;
	volume_in.open(volume_fname.c_str());
	string gas_fname = "gas";
	gas_fname  = RUN_path+gas_fname+num+ext;
	ifstream gas_in;
	gas_in.open(gas_fname.c_str());

	// you then loop through all these files, getting the important information
	// and the spacing
	// for the condition, you need:
	//
	// temperature (from master file)
	// pH (from ph file)
	// concentrations of primary species (from conc file)
	// volume fraction of minerals (from volume file)
	// surface area of minerals (this could be specific or bulk
	// exchange: this from master file
	// solid density: this from the master file

	char data_line[5000];
	string temp_string;
	string delim= " ";
	vector<string> line_words;
	vector<string> empty_str_vec;
	list< vector<double> >::iterator lv_iter;

	// get the pH
	pH_in.getline(data_line,5000);
	pH_in.getline(data_line,5000);
	while (pH_in.getline(data_line,5000))
	{
		temp_string = data_line;
		split_string(temp_string, delim, line_words);
		depth.push_back( atof(line_words[0].c_str() ) );
		pH_values.push_back( atof(line_words[1].c_str() ) );
		line_words = empty_str_vec;
	}
	n_conditions = depth.size();

	// first get the concentrations
	conc_in.getline(data_line,5000);
	conc_in.getline(data_line,5000);
	conc_in.getline(data_line,5000);
	while (conc_in.getline(data_line,5000))
	{
		temp_string = data_line;
		split_string(temp_string, delim, line_words);

		//cout << data_line << endl;

		lv_iter = concentrations.begin();
		int counter = 1;
		while(lv_iter!=concentrations.end() )
		{
			(*lv_iter).push_back( atof(line_words[counter].c_str()) );
			lv_iter++;
			counter++;
		}
		line_words = empty_str_vec;
	}

	// now replace concentrations with gas concentrations
	gas_in.getline(data_line,5000);
	gas_in.getline(data_line,5000);
	gas_in.getline(data_line,5000);
	int gcounter = 0;
	while (gas_in.getline(data_line,5000))
	{
		temp_string = data_line;

		split_string(temp_string, delim, line_words);

		lv_iter = concentrations.begin();
		lv_iter++; 	// advance to CO2
		//cout << "gas is: " << atof( line_words[1].c_str() ) << endl;
		(*lv_iter)[gcounter] = atof( line_words[1].c_str() );
		line_words = empty_str_vec;
		gcounter++;
	}

	//cout << "Line 352 did concentrations\n";

	// now get the mineral volume fractions
	int starting_cell = this_bin*n_cells_in_bin;
	
	volume_in.getline(data_line,5000);
	volume_in.getline(data_line,5000);
	volume_in.getline(data_line,5000);
	while (volume_in.getline(data_line,5000))
	{
		//cout << "data line: " << data_line << endl;
		temp_string = data_line;
		split_string(temp_string, delim, line_words);

    // the iterator has to be advanced to the starting cell
		lv_iter = mineral_vpercents.begin()+starting_cell;
		
		int counter = 1;
		while(lv_iter!=mineral_vpercents.end() )
		{
			//cout << "pushing back: " << atof(line_words[counter].c_str() ) << endl;
			(*lv_iter) = atof(line_words[counter].c_str() );
			lv_iter++;
			counter++;
		}
		line_words = empty_str_vec;
		//cout << endl << endl;
	}
	//cout << "Line 372 did volumes\n";

	// now get the mineral surface areas
	area_in.getline(data_line,5000);
	area_in.getline(data_line,5000);
	area_in.getline(data_line,5000);
	while (area_in.getline(data_line,5000))
	{
		temp_string = data_line;
		split_string(temp_string, delim, line_words);

    // the iterator has to be advanced to the starting cell
		lv_iter = mineral_ssa.begin()+starting_cell;
		
		int counter = 1;
		while(lv_iter!=mineral_ssa.end() )
		{
			(*lv_iter) = atof(line_words[counter].c_str() );
			lv_iter++;
			counter++;
		}
		line_words = empty_str_vec;
	}
	//cout << "Line 392 did mineral bulk surface area\n";

	// now get the rates
	rate_in.getline(data_line,5000);
	rate_in.getline(data_line,5000);
	rate_in.getline(data_line,5000);
	while (rate_in.getline(data_line,5000))
	{
		temp_string = data_line;
		split_string(temp_string, delim, line_words);

		lv_iter = reaction_rates.begin();
		int counter = 2;
		while(lv_iter!=reaction_rates.end() )
		{
			(*lv_iter).push_back( atof(line_words[counter].c_str() ) );
			lv_iter++;
			counter++;
		}
		line_words = empty_str_vec;
	}
	//cout << "Line 392 did mineral reaction rates\n";

	//int sz_pH = pH_values.size();
	//for (int i = 0; i< sz_pH; i++)
	//{
	//	cout << "line 412 pH["<<i<<"]: " << pH_values[i] << endl;
	//}

	// figure out how many spacings there are
	int n_depths = depth.size();
	bottom_depths.push_back(depth[0]*2.0);
	top_depths.push_back(0.0);
	double temp_spacings;
	for(int i = 1; i<n_depths; i++)
	{
		top_depths.push_back(bottom_depths[i-1]);
		temp_spacings = (depth[i]-top_depths[i])*2;
		bottom_depths.push_back(top_depths[i]+temp_spacings);

	}

	//for(int i = 0; i<n_depths; i++)
	//{
	//	cout << "depth: " << depth[i] << " top_depth: " << top_depths[i]
	//	     << " bottom depths: " << bottom_depths[i] << endl;
	//}

	rate_in.close();
	conc_in.close();
	volume_in.close();
	area_in.close();
	gas_in.close();
	pH_in.close();
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function looks at the depths and sets the string for the discretization
// line to reflect the cells defined but the top and bottom depths vectors
// these vectors are generated both by the particle model and by the
// parse_CRUNCH_files
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_engine::set_spacing_line_from_depths(vector<double>& top_depths,
												  vector<double>& bottom_depths)
{
	// first get the number of cells
	int n_cells = top_depths.size();

	vector<int> n_cells_in_spacing;
	vector<double> spacings;

	double temp_spacing;
	// loop through getting spacings
	spacings.push_back(bottom_depths[0]-top_depths[0]);
	n_cells_in_spacing.push_back(1);
	int s_counter = 0;
	for(int i = 1; i<n_cells; i++)
	{
		temp_spacing = bottom_depths[i]-top_depths[i];
		if ( fabs(temp_spacing-spacings[s_counter]) > 1e-9)
		{
			//cout << "temp_spacing: " << temp_spacing << " and old spacing: " 
      //     << spacings[s_counter] << endl;
			spacings.push_back(temp_spacing);
			n_cells_in_spacing.push_back(1);
			s_counter++;
		}
		else
		{
			n_cells_in_spacing[s_counter]++;
		}

	}

	int n_spacings = spacings.size();

	// now write the discretization line
	string discret_name = "xzones";
	string num;
	string space = " ";
	string new_line = discret_name;
	for(int i = 0; i<n_spacings; i++)
	{
		new_line+=space+itoa(n_cells_in_spacing[i])+space+dtoa(spacings[i]);
	}
	//cout << "LINE 504 new discretization line; " << new_line << endl;

	// search through list looking for discretization line
	list<string>::iterator l_iter;
	l_iter = master_infile.begin();
	size_t found = string::npos; 			// integer for testing if the command block has been found

	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(discret_name);
		//cout << "found is: " << found << endl;
  		if (found!=string::npos)
  		{
    		//cout << "spatial found at: " << int(found) << endl;

			// now replace the line
			(*l_iter) = new_line;
		}
		l_iter++;
	}
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function sets the concentration list vecs to some background concentration
//
// the parameters A and  B are for a fit to Ph data where pH = A ln(d) * B
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<double> CRUNCH_engine::set_up_pH_for_particle(
						vector<double>& top_depths,vector<double>& bottom_depths,
						double A, double B)
{
	int n_conditions = top_depths.size();
	vector<double> pH;
	double d_midpoint;
	for (int i = 0; i<n_conditions; i++)
	{
		d_midpoint = top_depths[i]+0.5*(bottom_depths[i]-top_depths[i]);
		pH.push_back(A*log(d_midpoint)+B);
		cout << "depth: " << d_midpoint <<" and pH: " << pH[i] << endl;

	}

	return pH;
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function gets default concentrations
//
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
list< vector<double> > CRUNCH_engine::get_default_concentrations(int n_ts,
							vector<double>& top_depths,vector<double>& bottom_depths)
{
	list< vector<double> > concentrations;
	list< vector<double> > updated_concentrations;
	vector<double> empty_vec;

	char data_line[5000];
	string temp_string;
	string delim= " ";
	vector<string> line_words;
	vector<string> empty_str_vec;
	list< vector<double> >::iterator lv_iter;

	int n_conditions = top_depths.size();


	string num = itoa(n_ts);
	string ext = ".out";
	string conc_fname = "totcon";
	conc_fname  = RUN_path+conc_fname+num+ext;
	ifstream conc_in;
	conc_in.open(conc_fname.c_str());
	string gas_fname = "gas";
	gas_fname  = RUN_path+gas_fname+num+ext;
	ifstream gas_in;
	gas_in.open(gas_fname.c_str());


	// initialize the concentration list
	int n_species = p_species_names.size();
	for(int i = 0; i<n_species; i++)
	{
		concentrations.push_back(empty_vec);
	}

	// first get the concentrations
	conc_in.getline(data_line,5000);
	conc_in.getline(data_line,5000);
	conc_in.getline(data_line,5000);
	
	// this third line contains the species names. You need to make sure that 
	// the CO2(aq) is in the correct place
	temp_string = data_line;
	split_string(temp_string, delim, line_words);	
	int n_line_words = int(line_words.size());
	int CO2_column = 2;
	for(int i = 0; i<n_line_words; i++)
	{
	  string thisword =  line_words[i];
    //string nospace = remove_if(thisword.begin(), thisword.end(), isspace);    
    string CO2word = "CO2(aq)";
    //cout << "i is: " << i << " and species is: " <<   thisword << " and word is: " << CO2word << endl;
    if(thisword == CO2word)
    {
      //cout << "Found aqueous CO2 in column: " << i << endl;
      CO2_column = i; 
    }
    
  }
  
  if (CO2_column != 2)
  {
	  cout << "CRUNCH_engine::get_default_concentrations WARNING CO2 not in 2nd position" 
	       << "check that your primary species has CO2 listed second in the crunch infile" << endl;
	}
	
	
	while (conc_in.getline(data_line,5000))
	{
		temp_string = data_line;
		split_string(temp_string, delim, line_words);

		//cout << data_line << endl;

		lv_iter = concentrations.begin();
		int counter = 1;
		while(lv_iter!=concentrations.end() )
		{
			(*lv_iter).push_back( atof(line_words[counter].c_str()) );
			lv_iter++;
			counter++;
		}
		line_words = empty_str_vec;
	}

	// now replace CO2 concentrations with gas concentrations
	gas_in.getline(data_line,5000);
	gas_in.getline(data_line,5000);
	gas_in.getline(data_line,5000);
	int gcounter = 0;
	while (gas_in.getline(data_line,5000))
	{
		temp_string = data_line;

		split_string(temp_string, delim, line_words);

		lv_iter = concentrations.begin();
		lv_iter++; 	// advance to CO2
		//cout << "gas is: " << atof( line_words[1].c_str() ) << endl;
		(*lv_iter)[gcounter] = atof( line_words[1].c_str() );
		line_words = empty_str_vec;
		gcounter++;
	}


	// now create a new list_vec that just populates a list vec of the appropriate size
	// with the concentrations from the first line.
	double temp_conc;

	lv_iter = concentrations.begin();
	while(lv_iter!=concentrations.end())
	{
		temp_conc = (*lv_iter)[0];
		vector<double> temp_vec;
		for (int i = 0; i<n_conditions; i++)
		{
			temp_vec.push_back(temp_conc);
		}
		updated_concentrations.push_back(temp_vec);
		lv_iter++;
	}

	return updated_concentrations;

}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function creates a CRUNCH infile using a number of different data elements
// these data elements are obtained from either an existing crunch run
// or a combination of a crunch run (to get solute concentrations) and a particle
// run
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_engine::create_CRUNCH_in_file(int& n_conditions, int n_bin,
						int cells_in_bin, vector<double>& pH_values,
						vector<double>& top_depths,vector<double>& bottom_depths,
						list< vector<double> >& concentrations, list< vector<double> >& mineral_vpercents,
						list< vector<double> >& mineral_ssa)
{

	// open the file
	ofstream CRUNCH_write;
	CRUNCH_write.precision(7);			// this precision is specifically chosen to
										// match that of CRUNCH in and out files.
  string bin_name = itoa(n_bin);
  bin_name = "_bin"+bin_name;
	string cmodel_name = RUN_path+"column_model"+bin_name+".in";									
	CRUNCH_write.open(cmodel_name.c_str());

	// you need to sweep through the list of strings associated with the master
	// file and change two things:
	//
	// 1) the spacings
	// 2) the profile information
	//
	// there is a function that changes the line in the infile with the spacings, we use it here
	set_spacing_line_from_depths(top_depths, bottom_depths);

	// now loop through the conditions, setting a new name
	vector<string> condition_names;
	vector<string> CRUNCH_file_condition_lines;
	string num;
	string space = " ";
	string weath_cell = "WeatheringCell";
	for(int i = 1; i<=n_conditions; i++)
	{
		num = itoa(i);
		condition_names.push_back(weath_cell+num);
		CRUNCH_file_condition_lines.push_back(weath_cell+num+space+num);
	}

	// place these name into the infile
	// search through list looking for discretization line
	list<string>::iterator l_iter;
	l_iter = master_infile.begin();
	size_t found = string::npos; 			// integer for testing if the command block has been found
	string icon = "INITIAL_CONDITIONS";
	// first find the intial condition
	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(icon);
		l_iter++;
	}
	found = string::npos; 		// reset found
	// now remove until you reach the end of the datablock
	string endname = "END";
	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(endname);
		if(found==string::npos)
		{
			l_iter = master_infile.erase(l_iter);
		}
	}

	// now start again and find the beginning of the datablock
	l_iter = master_infile.begin();
	found = string::npos;
	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(icon);
		l_iter++;
	}
	// now insert the condition names
	for(int i = 0; i<n_conditions; i++)
	{
		//cout << "i: " << i << " and condition: " << CRUNCH_file_condition_lines[i] << endl;
		l_iter = master_infile.insert(l_iter,CRUNCH_file_condition_lines[i]);
		l_iter++;
	}



	// you also need to update the boundary condition name with the final condition
	l_iter = master_infile.begin();
	found = string::npos;
	string end_bound = "X_end";
	string flux_name = "flux";
	string bound_name = end_bound+space+condition_names[n_conditions-1]+space+flux_name;
	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(end_bound);
		if(found!=string::npos)
		{
			(*l_iter) = bound_name;
		}
		l_iter++;
	}


	// now print the beginning information from the file
	l_iter = master_infile.begin();
	while(l_iter!=master_infile.end())
	{
		//cout << (*l_iter) << endl;
		CRUNCH_write << (*l_iter) << endl;
		l_iter++;
	}
	CRUNCH_write << endl << endl;

	string cond = "Condition ";
	string gas_line = "CO2(aq) CO2(g) ";
	list< vector<double> >::iterator v_iter;
	list< vector<double> >::iterator vssa_iter;
	// now loop through the conditions, printing each one


  int this_cell;
  int starting_cell = n_bin*cells_in_bin;
	for (int i = 0; i<n_conditions; i++)
	{
	
	  this_cell = i+starting_cell;
		CRUNCH_write << "Condition " << condition_names[i] << endl;
		CRUNCH_write << temperature_line << endl;
		CRUNCH_write << density_line << endl;
		CRUNCH_write << "pH " << pH_values[i] << endl;

		// now loop through primary species;
		l_iter = p_species_names.begin();		// the +1 is there because the hydrogen
		v_iter = concentrations.begin();		// ion is reported in the concentration
		l_iter++;								// file but does not go into
		v_iter++;								// the conditions
		CRUNCH_write << gas_line << (*v_iter)[i]  << endl;
		l_iter++;								// these are here because you don't
		v_iter++;								// print the CO2(aq) concnetration but rather
												// the gas concentration
		while( l_iter != p_species_names.end() )
		{
			CRUNCH_write << (*l_iter) << " " << (*v_iter)[i] << endl;

			l_iter++;
			v_iter++;
		}

		// now for the minerals
		// NOTE mineral percents are done by cell, whereas concentrations (above) 
		// are done by the number of the cell withing the column, with inices
		// from 0 to n_cells_in_bin
		l_iter = mineral_names.begin();
		v_iter = mineral_vpercents.begin();
		vssa_iter = mineral_ssa.begin();
		vector<double> m_vec;
		int counter = 0;
		while( l_iter != mineral_names.end())
		{
			// note the volume needs to be divided by 100 because it is reported in %
			CRUNCH_write << (*l_iter) << " " << ((*v_iter)[this_cell])*0.01 << " ssa " << (*vssa_iter)[this_cell] << endl;
			l_iter++;
			v_iter++;
			vssa_iter++;
			counter++;
		}

		// print stuff out:
		//l_iter = mineral_names.begin();
		//v_iter = mineral_vpercents.begin();
		//while( v_iter != mineral_vpercents.end())
		//{
		//	m_vec = (*v_iter);
		//	int sz = m_vec.size();
		//	cout << endl << "mineral is: " << (*l_iter) << endl;
		//	for (int i = 0; i<sz; i++)
		//	{
		//		cout << m_vec[i] << endl;
		//	}
		//	v_iter++;
		//	l_iter++;
		//}


		CRUNCH_write << "END" << endl;
		CRUNCH_write << endl << endl;

	}
	CRUNCH_write.close();
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_engine::call_CRUNCH()
{
	
  // this sets up the file for running crunch. NOTE: crunch has limited 
  // control over where you put the files, so the the pest file needs to be
  // in the same directory as the crunch executable. 
  ofstream pest_out;
	string full_pest_name = "PestControl.ant";
	string pest_name = RUN_path+ "column_model.in";
	pest_out.open(full_pest_name.c_str());
	pest_out << pest_name << endl;
	pest_out.close();

	// check to see if infile exists:
	ifstream c_in;
	c_in.open(pest_name.c_str());
	if (c_in.fail())
	{
		cout << "No infile for crunch!!\n";
		exit(1);
	}
	c_in.close();

	int i;
	cout << "Checking if processor is available..." << endl;
	if (system(NULL)) puts ("Ok");
	else exit (1);
	cout << "Executing CRUNCH...\n";
	
	
	string command_line_str = "cmd /c "+CRUNCH_path+"CrunchFlow2007";
	cout << "command_line_str is: " << command_line_str << endl;

	// this system call is buggy. It works on old laptop (probably an old version of cygwin)
	// with the latest version of cygwin it only works if CrunchFlow2007 is sitting in the
	// cygwin bin directory. However it appears to work with files in the same directory as the
	// volume particle executable
	//i=system("CrunchFlow2007");

	// here is another version. This tells cygwin to use the command interface.
	// it seems to run a bit slower than the above version but is more portable
	i=system(command_line_str.c_str());
	
	int n_ts = 1;
	move_CRUNCH_output_files(n_ts);
}



//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This moves the output files from crunch, which are dumped into the 
// directory of the executable, into the run directory. 
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_engine::move_CRUNCH_output_files(int n_ts)
{

  string num = itoa(n_ts);
	string ext = ".out";
	
	// get the areas
	string area_fname = "area";
	string area_fname_src  = area_fname+num+ext;
	string area_fname_dest  = RUN_path+area_fname+num+ext;
	ifstream src_area(area_fname_src.c_str());
	ofstream dst_area(area_fname_dest.c_str());
	dst_area << src_area.rdbuf();
	
	// get the concentrations
	string conc_fname = "totcon";
	string conc_fname_src  = conc_fname+num+ext;
	string conc_fname_dest  = RUN_path+conc_fname+num+ext;
	ifstream src_conc(conc_fname_src.c_str());
	ofstream dst_conc(conc_fname_dest.c_str());
	dst_conc << src_conc.rdbuf();	

	// get the pH
	string pH_fname = "pH";
	string pH_fname_src  = pH_fname+num+ext;
	string pH_fname_dest  = RUN_path+pH_fname+num+ext;
	ifstream src_pH(pH_fname_src.c_str());
	ofstream dst_pH(pH_fname_dest.c_str());
	dst_pH << src_pH.rdbuf();	

	// get the rate
	string rate_fname = "rate";
	string rate_fname_src  = rate_fname+num+ext;
	string rate_fname_dest  = RUN_path+rate_fname+num+ext;
	ifstream src_rate(rate_fname_src.c_str());
	ofstream dst_rate(rate_fname_dest.c_str());
	dst_rate << src_rate.rdbuf();	

	// get the volume
	string volume_fname = "volume";
	string volume_fname_src  = volume_fname+num+ext;
	string volume_fname_dest  = RUN_path+volume_fname+num+ext;
	ifstream src_volume(volume_fname_src.c_str());
	ofstream dst_volume(volume_fname_dest.c_str());
	dst_volume << src_volume.rdbuf();	
	
	// get the gas
	string gas_fname = "gas";
	string gas_fname_src  = gas_fname+num+ext;
	string gas_fname_dest  = RUN_path+gas_fname+num+ext;
	ifstream src_gas(gas_fname_src.c_str());
	ofstream dst_gas(gas_fname_dest.c_str());
	dst_gas << src_gas.rdbuf();	
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function get the molar weight and molar volume from the database file
void CRUNCH_engine::get_mineral_properties()
{
	//cout << "line 832 getting mineral properties" << endl;

	// reset the molar volume and weight vectors
	vector<double> empty_vec;
	molar_volume = empty_vec;
	molar_weight = empty_vec;

	// initiate some iterators for looking through lists
	list<string>::iterator l_iter;
	list<string>::iterator l_min_name_iter;

	// the first thing you do is get the name of the database file
	l_iter = master_infile.begin();
	int found = string::npos; 			// reset found
	string dbase_name = ".dbs";
	string delim = " ";
	string dbase_fname;

	// advance through file until you find the primary species
	while(l_iter != master_infile.end() && found == string::npos)
	{
		found=(*l_iter).find(dbase_name);
		vector<string> line_words;
		if(found!=string::npos)
		{
			//cout << "line 855 found database " << endl;
			split_string((*l_iter), delim, line_words);
			dbase_fname =line_words[1];
		}
		l_iter++;
	}
	cout << "database name is: " << dbase_fname << endl;
	
	// update the dbase so that it looks in the crunch folder
	string dbase_fname_with_path = dbase_fname;

	// load the database into a list
	list<string> dbase_list;
	ifstream dbase_in;
	dbase_in.open(dbase_fname_with_path.c_str());
	char dbase_line[5000];
	string linestr;
	// then load in the database file into a list of strings
	while(dbase_in.getline(dbase_line,5000))
	{
		linestr = dbase_line;
		dbase_list.push_back(linestr);
	}
	dbase_in.close();

	// first get the number of temperature points
	l_iter = dbase_list.begin();
	found = string::npos;
	int n_temp_points;
	string temp_p = "temperature";
	while(l_iter != dbase_list.end() && found == string::npos)
	{
		found=(*l_iter).find(temp_p);
		if (found!=string::npos)
		{
			vector<string> line_words;
			split_string((*l_iter), delim, line_words);
			n_temp_points = atoi(line_words[2].c_str() );
		}
		l_iter++;
	}
	//cout << "Line 893, number of temperature points is: " << n_temp_points << endl;

	// now loop through mineral names getting the database information for each name in turn
	int n_species;
	l_min_name_iter = mineral_names.begin();
	while(l_min_name_iter != mineral_names.end())
	{
		l_iter = dbase_list.begin();
		found = string::npos;
		while(l_iter != dbase_list.end() && found == string::npos)
		{
			found=(*l_iter).find( (*l_min_name_iter) );
			if (found!=string::npos)
			{
				vector<string> line_words;
				split_string((*l_iter), delim, line_words);

				//cout << "LINE 892 data line is: " << (*l_iter) << endl;
				n_species = atoi(line_words[2].c_str() );
				//cout << "line 912; number of species is: " << n_species << endl;
				int mwgt = 3+2*n_species+n_temp_points;
				molar_weight.push_back( atof(line_words[mwgt].c_str() ));
				molar_volume.push_back( atof(line_words[1].c_str() ) );

			}

			l_iter++;
		}
		l_min_name_iter++;
	}

	// now loop through minerals printing out their properties;
	int counter = 0;
	l_min_name_iter = mineral_names.begin();
	while(l_min_name_iter != mineral_names.end())
	{
		cout << (*l_min_name_iter)
			//<< " ssa: " << ssa[counter]
			<< " m_vol: " << molar_volume[counter]
		     << " m_wght: " << molar_weight[counter] << endl;
		l_min_name_iter++;
		counter++;
	}

}



//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this prints the master file to screen
// (mostly this is used for bug checking)
void CRUNCH_engine::print_master()
{
	list<string>::iterator l_iter;
	l_iter = master_infile.begin();
	int counter = 1;
	while(l_iter != master_infile.end())
	{
		cout << "L" << counter << ": " << (*l_iter) << endl;
		counter++;
		l_iter++;
	}
}
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#endif
