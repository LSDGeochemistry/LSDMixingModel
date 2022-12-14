// CRN_tParticle_bins.h
//

#ifndef CRUNCH_engine_H
#define CRUNCH_engine_H
#include<iostream>
#include <list>
#include <string>
#include <vector>
#include "LSDStatsTools.hpp"
#include "VolumeParticleInfo.hpp"
using namespace std;


/// Object to interface with CRUNCHFLOW. It manages the CRUNCH runs, which are separate compiled code. 
/// This involves manupulating parameter files and making system calls
class CRUNCH_engine
{
  public:
	
    /// @brief Default constructor method used to create a CRUNCH engine object
    /// @author SMM
    /// @date 01/01/11
    CRUNCH_engine()							{ create(); }
    
    /// @brief Constructor for CRUNCH engine that reads in path and filenames
    /// @param crunch_path The directory where the crunch files are stored
    /// @param run_path the driectory in which the crunch executable is stored
    /// @param master_filename Then name (without path) of the "master" crunch parameters
    ///   that will be read as default parameters for crunch runs
    /// @author SMM
    /// @date 01/01/11    
	CRUNCH_engine(string crunch_path, string run_path, string master_filename)	
                  { create(crunch_path, run_path, master_filename); }

    /// @brief this removes control characters from the mineral and primary species
    /// names
    /// @author SMM
    /// @date 10/08/2014
    void remove_control_characters_from_species_and_minerals();

	/// @brief This prints all the crunch parameters to a file called
    /// @author SMM
    /// @date 01/01/2011
	void print_master();

	/// @brief Reads data from CRUNCH output
    /// @detail Doesn't return anything but instead overwrites input objects
    /// @param n_ts The timestep number
    /// @param n_conditions the number of species conditions. Is replaces after parsing the files
    /// @param this_bin The number of the particle bin (a crunch engine runs in each vertical bin)
    /// @param n_cells_in_bin The number of vertical cells in the bin: the crunch cells correspond to volume elements in the flow tub
    /// @param pH_value The pH values, duh
    /// @param spacings How far apart the cells are vertically (this can be irregular)
    /// @param top_depths The depth (0 at surface, increasing downward) of the cells
    /// @param botton_depths The bottom depths
    /// @param concentrations The concentrations of all minerals and species in all cells (the list contains the species and the vectors are the concentrations as fxn depth)
    /// @param mineral_vpercents Volumne percentatge of the solid minerals
    /// @param minera_ssa Mineral specific surface areas (in Units?)
    /// @param reaction_rate The reaction rates (in Units?) over the time step
    /// @author SMM
    /// @date 01/01/2011
	void parse_CRUNCH_files(int n_ts, int& n_conditions,
	          int this_bin, int n_bins, int n_cells_in_bin,
						vector<double>& pH_values, vector<double>& spacings,
						vector<double>& top_depths, vector<double>& bottom_depths,
						list< vector<double> >& concentrations, list< vector<double> >& mineral_vpercents,
						list< vector<double> >& mineral_ssa, list< vector<double> >& reaction_rates);

	/// @brief this changes a specific paramter in the infile
    /// @param value_name The name of the value to be changed
    /// @param new_value The quantitiy of the new value
    /// @author SMM
    /// @date 01/01/2011
	void modify_CRUNCH_value(string value_name, double new_value);

	/// this gets information about the parent material (mostly for getting the ssa)
    /// @detail Reads the parent material file and does some volume calcualtions, then updates data members in the crunch engine object
    /// @author SMM
    /// @date 01/01/2011
	void parse_parent_material_file();

    
	/// @brief Creates a crunch input file from various input values
    /// @param n_conditions the number of species conditions. Is replaces after parsing the files
    /// @param n_cells_in_bin The number of vertical cells in the bin: the crunch cells correspond to volume elements in the flow tub
    /// @param pH_value The pH values, duh
    /// @param top_depths The depth (0 at surface, increasing downward) of the cells
    /// @param botton_depths The bottom depths. Note that spacings are derived from the top and bottom depths
    /// @param concentrations The concentrations of all minerals and species in all cells (the list contains the species and the vectors are the concentrations as fxn depth)
    /// @param mineral_vpercents Volumne percentatge of the solid minerals
    /// @param minera_ssa Mineral specific surface areas (in Units?)
    /// @author SMM
    /// @date 01/01/2011       
    void create_CRUNCH_in_file(int& n_conditions, int n_bin,
						int cells_in_bin, vector<double>& pH_values,
						vector<double>& top_depths,vector<double>& bottom_depths,
						list< vector<double> >& concentrations, list< vector<double> >& mineral_vpercents,
						list< vector<double> >& mineral_ssa);

	/// @brief this sets the line in the infile that states the different spacings of the nodes
    /// @param top_depths The depth (0 at surface, increasing downward) of the cells
    /// @param botton_depths The bottom depths. Note that spacings are derived from the top and bottom depths
    /// @author SMM
    /// @date 01/01/2011   
    void set_spacing_line_from_depths(vector<double>& top_depths,
												  vector<double>& bottom_depths);

    /// @brief This uses a simple function (A+log(depth)+B) for establishing pH
    /// @param top_depths The depth (0 at surface, increasing downward) of the cells
    /// @param botton_depths The bottom depths. Note that spacings are derived from the top and bottom depths
    /// @param A parameter in equation (see above)
    /// @param B parameter in equation (see above)
    /// @author SMM
    /// @date 01/01/2011     
	vector<double> set_up_pH_for_particle(
						vector<double>& top_depths,vector<double>& bottom_depths,
						double A, double B);
    
    /// @brief Gets the concentration from a crunchflow run, but doesn't ingest any other parameters
    /// @param n_ts The timestep number
    /// @param top_depths The depth (0 at surface, increasing downward) of the cells
    /// @param botton_depths The bottom depths. Note that spacings are derived from the top and bottom depths
    /// @author SMM
    /// @date 01/01/2011   
	list< vector<double> > get_default_concentrations(int n_ts,
							vector<double>& top_depths, vector<double>& bottom_depths);

	/// @brief get_mineral_properties: loads things like molar weight, etc;
    /// @author SMM
    /// @date 01/01/2011   
	void get_mineral_properties();
	
	/// @brief gets the names of the primary species
    /// @return A list of the primary species
    /// @author SMM
    /// @date 01/01/2011  
	list<string> get_primary_species_names()   { return p_species_names; }


	/// @detail this calls crunch using a system call. The infile must exist!
    /// @param bin_number The bin number
    /// @author SMM
    /// @date 01/01/2011     
	void call_CRUNCH(int bin_number);
	
	/// @brief this little script copies the output of crunch to the run folder
	/// @param n_ts is the number of the timestep
	/// @author SMM
	/// @date 05/08/2014
    void move_CRUNCH_output_files(int n_ts);	

    /// @brief this little script copies the output of crunch to the run folder
    /// it also appends the bin number to the filename
	/// @param n_ts is the number of the timestep
    /// @param bn The bin number
	/// @author SMM
	/// @date 05/08/2014
    void move_CRUNCH_output_files_with_bin_number(int n_ts, int bn);

	private:
	void create();
	void create(string crunch_path, string run_path, string master_filename);

	// this populates the mineral names and species names
	void get_minerals_and_species();
    void get_constant_conditions();

    string CRUNCH_path;
    string RUN_path;

	list<string> master_infile;
	list<string> p_species_names;        // names of the species
	list<string> mineral_names;          // names of the minerals

	string temperature_line;
	string density_line;
	vector<string> exchange_lines;

	vector<double> init_V_fracs;        // dimensionless
	vector<double> ssa;					// surface area in m^2/m^3 (I think...check)
	vector<double> molar_volume;		// molar volume in cm^3
	vector<double> molar_weight;		// molar weight in g



};


#endif
