// CRN_tParticle_bins.h
//

#ifndef CRUNCH_engine_H
#define CRUNCH_engine_H

#include <list>
#include <string>
#include <vector>
using namespace std;



class CRUNCH_engine
{
	public:
	CRUNCH_engine()							{ create(); }
	CRUNCH_engine(string crunch_path, string run_path, string master_filename)	
                  { create(crunch_path, run_path, master_filename); }

  /// @brief this removes control characters from the mineral and primary species
  /// names
  /// @author SMM
  /// @date 10/08/2014
  void remove_control_characters_from_species_and_minerals();

	// printing functions
	void print_master();

	// reading data from CRUNCH output
	void parse_CRUNCH_files(int n_ts, int& n_conditions,
	          int this_bin, int n_bins, int n_cells_in_bin,
						vector<double>& pH_values, vector<double>& spacings,
						vector<double>& top_depths, vector<double>& bottom_depths,
						list< vector<double> >& concentrations, list< vector<double> >& mineral_vpercents,
						list< vector<double> >& mineral_ssa, list< vector<double> >& reaction_rates);

	// this changes a specific paramter in the infile
	void modify_CRUNCH_value(string value_name, double new_value);

	// this gets information about the parent material (mostly for getting the ssa)
	void parse_parent_material_file();

  void create_CRUNCH_in_file(int& n_conditions, int n_bin,
						int cells_in_bin, vector<double>& pH_values,
						vector<double>& top_depths,vector<double>& bottom_depths,
						list< vector<double> >& concentrations, list< vector<double> >& mineral_vpercents,
						list< vector<double> >& mineral_ssa);

	// this sets the line in the infile that states the different spacings of the nodes
	void set_spacing_line_from_depths(vector<double>& top_depths,
												  vector<double>& bottom_depths);

	vector<double> set_up_pH_for_particle(
						vector<double>& top_depths,vector<double>& bottom_depths,
						double A, double B);
	list< vector<double> > get_default_concentrations(int n_ts,
							vector<double>& top_depths, vector<double>& bottom_depths);

	// get_mineral_properties: loads things like molar weight, etc;
	void get_mineral_properties();
	
	// gets the names of the primary species
	list<string> get_primary_species_names()   { return p_species_names; }


	// this calls crunch. The infile must exist!
	void call_CRUNCH(int bin_number);
	
	// this little script copies the output of crunch to the run folder
	// n_ts is the number of the timestep
  void move_CRUNCH_output_files(int n_ts);	

	private:
	void create();
	void create(string crunch_path, string run_path, string master_filename);

	// this populates the mineral names and species names
	void get_minerals_and_species();
	void get_constant_conditions();

  string CRUNCH_path;
  string RUN_path;

	list<string> master_infile;
	list<string> p_species_names;
	list<string> mineral_names;

	string temperature_line;
	string density_line;
	vector<string> exchange_lines;

	vector<double> init_V_fracs;
	vector<double> ssa;					// surface area in m^2/m^3 (I think...check)
	vector<double> molar_volume;			// molar volume in cm^3
	vector<double> molar_weight;		// molar weight in g



};


#endif
