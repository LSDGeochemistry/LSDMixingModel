// CRUNCH_bins.cpp
// implementation of the CRUNCH_bins object
// it is responsible for aggregating geochemical data
// across multiple columns and also includes printing 
// capabilities


#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include <cctype>
#include "mathutil.hpp"
#include "flowtube.hpp"
#include "CRN_parameters.hpp"
#include "CRN_tParticle_bins.hpp"
#include "CRUNCH_bins.hpp"
#include "VolumeParticleInfo.hpp"
#include "CRN_tParticle_bins.hpp"
#include "CRUNCH_engine.hpp"
using namespace std;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::create()
{
  cout << "You need to at least provide the number of bins and cells. " << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function:
// this just has some information about the cell locations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::create(int n_pdz_per_bin, int n_caz_per_bin, 
                 double b_depth, VolumeParticleInfo this_vpi)
{
  n_pdz_cells_per_bin = n_pdz_per_bin;
  n_caz_cells_per_bin = n_caz_per_bin;
  bottom_depth = b_depth;
  vpi = this_vpi;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// The create function:
// this just has some information about the cell locations
// the data is stored in vectors, you need to get the index into
// the vector of from the bin and cell number
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int CRUNCH_bins::retrieve_index_of_cell(int bin, int cell)
{
  int index_into_cell;
  if (cell < n_caz_cells_per_bin+n_pdz_cells_per_bin)
  {
    index_into_cell = (n_caz_cells_per_bin+n_pdz_cells_per_bin)*bin+cell;
  }
  else
  {
    cout << "Your cell is too deep. Returning top cell "<< endl;
    index_into_cell = (n_caz_cells_per_bin+n_pdz_cells_per_bin)*bin; 
  }
  return index_into_cell;

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function returns prints out the edge coordinates of a cell
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::cell_location_to_screen(int bin, int cell)
{
  // get the cell index
  int cell_index =  retrieve_index_of_cell(bin, cell);
  int this_sz;
  vector<double> verts_z;
  vector<double> verts_s;
  vector<double> verts_d;
  
  if(cell_data_map.find("verts_s") == cell_data_map.end()) 
  {
    cout << "CRUNCH_bins::cell_location_to_screen " 
         << "You haven't got any cell information! " << endl;
  }  
  else
  {
    this_sz =  cell_index_map["cell_node1"].size();
    cout << "Size if cell index is: " << this_sz 
         << " and cell index is: " << cell_index << endl;
    if (cell_index < this_sz)
    {
      verts_z =  cell_data_map["verts_z"];
      verts_s =  cell_data_map["verts_s"];
      verts_d =  cell_data_map["verts_d"];
      
      int ci1 = cell_index_map["cell_node1"][cell_index];
      int ci2 = cell_index_map["cell_node2"][cell_index];
      int ci3 = cell_index_map["cell_node3"][cell_index];
      int ci4 = cell_index_map["cell_node4"][cell_index];
    
      //cout << "bin: " << bin << " cell: " << cell << endl;
      //cout << "s:\t" << verts_s[ci1] << "\t" << verts_s[ci2] << "\t" 
      //               << verts_s[ci3] << "\t" << verts_s[ci4] << endl; 
      //cout << "z:\t" << verts_z[ci1] << "\t" << verts_z[ci2] << "\t" 
      //               << verts_z[ci3] << "\t" << verts_z[ci4] << endl;       
      //cout << "d:\t" << verts_d[ci1] << "\t" << verts_d[ci2] << "\t" 
      //               << verts_d[ci3] << "\t" << verts_d[ci4] << endl;                                       
    }
    else
    {
      cout << "CRUNCH_bins::cell_location_to_screen "
           << "Cell index appears to be out of bounds. " << endl;
    }
      
  }
  
} 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function assumes cells have been assigned to particles
// it then gets the mass fractions of the minerals as well as the 
// depletion fractions
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::populate_cells_with_mfracs_and_depletion(CRN_tParticle_bins& CRN_tPb)
{
  // get the number of bins from the CRN_tParticle_bin_object
  n_bins = CRN_tPb.get_n_bins();

  // calculate the total number of cells 
  total_cells = (n_caz_cells_per_bin+n_pdz_cells_per_bin)*n_bins;  
  
  vector<double> empty_vec;
  vector< list< vector<double> > > empty_vlv(n_bins);
	vec_mineral_mfracs_old = empty_vlv;
 	vec_mineral_depletion_old = empty_vlv;
   
  // get the locations of the cells and the cell indices. 
  // this function also updates the cell indices of the particles
  for (int bn = 0; bn<n_bins; bn++)
  {
    //cout << "size vmvo: " << vec_mineral_vpercents_old.size() << " and bn: " << bn << endl;
    // now collect the data from from the particles
    // this collects data from each cell in this bin
	  list< vector<double> > mineral_mfracs_old;
 	  list< vector<double> > mineral_depletion_old;
	  CRN_tPb.get_mineral_mass_loss_and_mfracs_volumetric(bn, n_pdz_cells_per_bin, 
                           n_caz_cells_per_bin, vpi, mineral_mfracs_old, 
                           mineral_depletion_old);   
    vec_mineral_mfracs_old[bn] = mineral_mfracs_old;
    vec_mineral_depletion_old[bn] = mineral_depletion_old;                 
  }    
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function both assigns cells to the particles and
// calculates several cell based metrics, mainly used in crunch simulations
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::populate_cells_with_geochemical_data_from_CRNtPb(flowtube& ft,
                                                CRN_tParticle_bins& CRN_tPb)
{


  // get the number of bins from the CRN_tParticle_bin_object
  n_bins = CRN_tPb.get_n_bins();
  
  // calculate the total number of cells 
  total_cells = (n_caz_cells_per_bin+n_pdz_cells_per_bin)*n_bins;
  
  // the bottom and top locations of all the cells
  // the reason why there is not one master vector for all the bins
  // is because PDZ thickness can change spatially. 
  vector<double> d_bottom_locs(total_cells);
  vector<double> d_top_locs(total_cells);   
  
  // the bottom and top locations of individual bins
  // these vectors will get replaced in each bin and are temporary
  vector<double> bin_d_bottom_locs;
  vector<double> bin_d_top_locs;

  // some temporary vectors for holding cell information
  vector<double> verts_s; 
  vector<double> verts_z;
  vector<double> verts_d;
  
  // some temporary vectors for holding cell indices
  vector<int> cell_node1;
  vector<int> cell_node2;
  vector<int> cell_node3;
  vector<int> cell_node4;
  
  // First partition the data into cells
  CRN_tPb.update_particles_cell_index(ft,
								n_pdz_cells_per_bin, n_caz_cells_per_bin,
								bottom_depth,
								verts_s, verts_z, verts_d,
							  cell_node1, cell_node2, cell_node3, cell_node4); 			

  cell_data_map["verts_s"] = verts_s;
  cell_data_map["verts_z"] = verts_z;
  cell_data_map["verts_d"] = verts_d;
   
  cell_index_map["cell_node1"] = cell_node1;
  cell_index_map["cell_node2"] = cell_node2;
  cell_index_map["cell_node3"] = cell_node3;
  cell_index_map["cell_node4"] = cell_node4;

  vector<double> empty_vec;
  vector< list< vector<double> > > empty_vlv(n_bins);
	vec_mineral_vpercents_old = empty_vlv;
 	vec_mineral_ssa_old = empty_vlv;
 	vec_mineral_mass_old = empty_vlv;
 	vec_mineral_surface_area_old = empty_vlv;
   
  // get the locations of the cells and the cell indices. 
  // this function also updates the cell indices of the particles
  for (int bn = 0; bn<n_bins; bn++)
  {
    //cout << "size vmvo: " << vec_mineral_vpercents_old.size() << " and bn: " << bn << endl;
    // now collect the data from from the particles
    // this collects data from each cell in this bin
	  list< vector<double> > mineral_vpercents_old;
 	  list< vector<double> > mineral_ssa_old;
 	  list< vector<double> > mineral_mass_old;
 	  list< vector<double> > mineral_surface_area_old;
	  CRN_tPb.get_data_by_cell_volumetric_for_CRUNCH(bn,n_pdz_cells_per_bin, 
                n_caz_cells_per_bin,bottom_depth,
								verts_s, verts_z, verts_d,
								cell_node1, cell_node2, cell_node3, cell_node4, vpi,
								mineral_vpercents_old,mineral_ssa_old,
								mineral_surface_area_old,mineral_mass_old);   
    vec_mineral_vpercents_old[bn] = mineral_vpercents_old;
    vec_mineral_ssa_old[bn] = mineral_ssa_old;
    vec_mineral_mass_old[bn] = mineral_mass_old;
    vec_mineral_surface_area_old[bn] = mineral_surface_area_old;
                  
  }
  
}                                                


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// vtk cell printing
// This function prints the cell data to a vtk file that can be read by prgrams 
// such as paraview
//
// it takes a reference frame switch. If this == 1, then the reference frame
// is fixed to the surface elevation (it is a plot of depths with the
// topography removed)
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRUNCH_bins::vtk_print_cell_header(int reference_frame_switch, 
                                     ofstream& vtk_cell_out)
{
	// now print the vtk file
	// find the number of particles
	
	//cout << "verts size: " << n_verts << endl;

  // check if the data exists
  if(cell_data_map.find("verts_s") == cell_data_map.end()) 
  {
    cout << "CRUNCH_bins::cell_location_to_screen " 
         << "You haven't got any cell information! " << endl;
  }  
  else
  {
    // get the number of vertices
    int n_verts = cell_data_map["verts_z"].size();
    vector<double> verts_s = cell_data_map["verts_s"];
    vector<double> verts_d = cell_data_map["verts_d"];
    vector<double> verts_z = cell_data_map["verts_z"];
    
    // now get the indices into the cell nodes
    vector<int> cell_node1 = cell_index_map["cell_node1"];
    vector<int> cell_node2 = cell_index_map["cell_node2"];
    vector<int> cell_node3 = cell_index_map["cell_node3"];
    vector<int> cell_node4 = cell_index_map["cell_node4"];
    
    
    int n_cells = (n_caz_cells_per_bin+n_pdz_cells_per_bin)*n_bins;
  
    // header lines of the vtk file
  	vtk_cell_out << "# vtk DataFile Version 2.0" << endl 
            << "Unstructured Grid Ptrack_cells"
  	        << endl << "ASCII" << endl << endl 
            << "DATASET UNSTRUCTURED_GRID" << endl
  	        << "POINTS " << n_verts << " float" << endl;
  	if (reference_frame_switch == 1)
  	{
  		for (int i = 0; i< n_verts; i++)
  		{
  			vtk_cell_out << verts_s[i] << " " << -verts_d[i] << " 0.0" <<endl;
  		}
  	}
  	else
  	{
  		for (int i = 0; i< n_verts; i++)
  		{
  			vtk_cell_out << verts_s[i] << " " << verts_z[i] << " 0.0" <<endl;
  		}
  	}

  	vtk_cell_out << endl << "CELLS "<<n_cells << " " << n_cells*5 << endl;
  	for (int i = 0; i< n_cells; i++)
  	{
  		vtk_cell_out << 4 << " " << cell_node1[i] << " " << cell_node2[i] << " "
  		        << cell_node3[i] << " " << cell_node4[i] << endl;
  	}
  
  	vtk_cell_out << endl << "CELL_TYPES " << n_cells << endl;
  	for (int i = 0; i< n_cells; i++)
  	{
  		vtk_cell_out << "9" << endl;
  	}
    
  	vtk_cell_out << endl << "CELL_DATA " << n_cells << endl;
  	vtk_cell_out << "SCALARS SOIL_OR_PARENT int 1" << endl 
                 << "LOOKUP_TABLE default" << endl; 
    int s_or_p;
    int cell_counter = 0;
    for (int bn = 0; bn<n_bins; bn++)
    {
      for(int cib = 0; cib< (n_pdz_cells_per_bin+n_caz_cells_per_bin); cib++)
      {
        if(cib <n_pdz_cells_per_bin)
        {
          s_or_p = 1;
        }  
        else
        {
          s_or_p = 0;
        }
        vtk_cell_out << s_or_p << endl;
      }
      
    }
    
  
  	vtk_cell_out << "SCALARS CELL_NUM int 1" 
                 << endl << "LOOKUP_TABLE default" << endl;
  	for (int i = 0; i< n_cells; i++)
  	{
  		vtk_cell_out << i <<endl;
  	}
  }
}
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This prints data from a map to a vtk file
//
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::vtk_print_cell_from_map_of_vectors(ofstream& vtk_cell_out, 
                                       map<string, vector<double> >& data_map)
{
  int n_cells_in_bin = n_pdz_cells_per_bin+n_caz_cells_per_bin;
  
  map<string, vector<double> >::iterator map_iter;  // the iterator for the map object
  
  map_iter = data_map.begin();
  
  while(map_iter != data_map.end())
  {
    string this_key = map_iter->first;
    vector<double> this_data = map_iter->second;
    
    if(total_cells != int(this_data.size()))
    {
      cout << "LINE 350 trying to print map data but the data is not the same "
          << "size as the number of cells" << endl;
    }
    else
    {
      vtk_cell_out << "SCALARS " << this_key << " float 1" << endl 
                   << "LOOKUP_TABLE default" << endl;
       for(int cn = 0; cn<total_cells; cn++)
       {
         vtk_cell_out << this_data[cn] << endl;
       }                                
    }
    map_iter++;
  }   
}                                                    
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function parses specific vec list vecs
// for vtk printing
//
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::vtk_print_cell_mineral_solid_state(ofstream& vtk_cell_out)
{
  // test if the geochem data has been derived
  //int n_cells_in_bin = n_pdz_cells_in_bin+n_caz_cells_in_bin;
  //list< vector<double> > this_lv = vec_mineral_vpercents_old[0];
  //vector<double> = 
  
  // get the names of the minerals
  list<string> mineral_names = get_names_of_minerals();
  
  // get the maps
  string mvname = "Mineral_volume_percents"; 
  map< string, vector<double> > mvpercents =  parse_vec_list_vec_to_vec_map(mvname, 
                            mineral_names, vec_mineral_vpercents_old);
  string mssaname = "Mineral_specific_surface_area_m2perg"; 
  map< string, vector<double> > mssafracs =  parse_vec_list_vec_to_vec_map(mssaname, 
                            mineral_names, vec_mineral_ssa_old);
  string mmname = "Mineral_mass_kg"; 
  map< string, vector<double> > mmass =  parse_vec_list_vec_to_vec_map(mmname, 
                            mineral_names, vec_mineral_mass_old);
  string msaname = "Mineral_surface_area_m2"; 
  map< string, vector<double> > msafracs =  parse_vec_list_vec_to_vec_map(msaname, 
                            mineral_names, vec_mineral_surface_area_old);                                                                                    

  // now print the map data to the vtk file
  vtk_print_cell_from_map_of_vectors(vtk_cell_out, mvpercents);
  vtk_print_cell_from_map_of_vectors(vtk_cell_out, mssafracs);
  vtk_print_cell_from_map_of_vectors(vtk_cell_out, mmass);
  vtk_print_cell_from_map_of_vectors(vtk_cell_out, msafracs);
}
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function parses specific vec list vecs
// for vtk printing
//
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::vtk_print_cell_mineral_mfracs_and_depletion(ofstream& vtk_cell_out)
{
 
  // get the names of the minerals
  list<string> mineral_names = get_names_of_minerals();
  
  // get the maps
  string mfracsname = "Mineral_mass_fraction"; 
  map< string, vector<double> > mfracs =  parse_vec_list_vec_to_vec_map(mfracsname, 
                            mineral_names, vec_mineral_mfracs_old );
  string mdepletionname = "Mineral_depletion_ratio"; 
  map< string, vector<double> > mdepletion =  parse_vec_list_vec_to_vec_map(mdepletionname, 
                            mineral_names, vec_mineral_depletion_old);
                                                                                   
  // now print the map data to the vtk file
  vtk_print_cell_from_map_of_vectors(vtk_cell_out, mfracs);
  vtk_print_cell_from_map_of_vectors(vtk_cell_out, mdepletion);
 
}
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function prints the mass fractions to a vtk file
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::vtk_print_cell_mineral_fractions(ofstream& vtk_cell_out, 
                    map< string, vector<double> >& raw_data)
{
  int n_cells_in_bin = n_pdz_cells_per_bin+n_caz_cells_per_bin;
  
  map<string, vector<double> >::iterator map_iter;  // the iterator for the map object
  
  map_iter = raw_data.begin();
  vector<double> this_data = map_iter->second;
  int n_nodes = int(this_data.size());
  vector<double> total_data(n_nodes,0.0);
  
  // first loop through the data collecting the total
  while(map_iter != raw_data.end())
  {
    string this_key = map_iter->first;
    vector<double> this_data = map_iter->second;
    
    if(total_cells != int(this_data.size()))
    {
      cout << "LINE 426 trying to print map data but the data is not the same "
          << "size as the number of cells" << endl;
    }
    else
    {
      for(int cell = 0; cell<total_cells; cell++)
      {
        total_data[cell]+=this_data[cell];
      }                         
    }
    map_iter++;
  }   

  map_iter = raw_data.begin();
  while(map_iter != raw_data.end())
  {
    string this_key = map_iter->first;
    vector<double> this_data = map_iter->second;
    this_key = this_key+"_dimensionless_fraction";
    
    if(total_cells != int(this_data.size()))
    {
      cout << "LINE 350 trying to print map data but the data is not the same "
          << "size as the number of cells" << endl;
    }
    else
    {
      vtk_cell_out << "SCALARS " << this_key << " float 1" << endl 
                   << "LOOKUP_TABLE default" << endl;
       for(int cn = 0; cn<total_cells; cn++)
       {
         if (total_data[cn] == 0)
         {
           vtk_cell_out << 0 << endl;  
         }
         else
         {
           vtk_cell_out << this_data[cn]/total_data[cn] << endl;
         }
       }                                
    }
    map_iter++;
  }   
}                   
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function parses specific vec list vecs
// for vtk printing
//
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::vtk_print_cell_CRUNCH_data(ofstream& vtk_cell_out, 
                     CRUNCH_engine& Ceng)
{
  // test if the geochem data has been derived
  //int n_cells_in_bin = n_pdz_cells_in_bin+n_caz_cells_in_bin;
  //list< vector<double> > this_lv = vec_mineral_vpercents_old[0];
  //vector<double> = 
  
  // get the names of the minerals
  list<string> mineral_names = get_names_of_minerals();
  
  // get the names of the primary species
  list<string> pspecies_names = get_names_of_primary_species(Ceng);
  
  //list<string>::iterator siter;
  
  //siter = mineral_names.begin();
  //while(siter!=mineral_names.end())
  //{
  //  cout << "Mineral is: " << (*siter) << endl;
  //  siter++;
  //}
  
  //siter = pspecies_names.begin();
  //while(siter!=pspecies_names.end())
  //{
  //  cout << "Species is: " << (*siter) << endl;
  //  siter++;  
  //}
  
  
  // get the maps
  string scname = "Solute_Concentration"; 
  map< string, vector<double> > solute_conc =  
            parse_CRUNCH_vec_list_vec_to_vec_map(scname, 
                            pspecies_names, vec_new_conc);
  
  
  //cout << "getting reaction rates " << endl;
  string rxnname = "Reaction_rates_Moles_per_L_porusmedium_per_second"; 
  map< string, vector<double> > rxn_rates =  
            parse_CRUNCH_vec_list_vec_to_vec_map(rxnname, 
                            mineral_names, vec_new_rxn_rates); 
  
  //list< vector<double> > this_rxn = vec_new_rxn_rates[0];                          
  //list< vector<double> >::iterator l_iter;
  //l_iter = this_rxn.begin();
  //vector<double> these_rates = (*l_iter);
  //for(int i = 0; i< int(these_rates.size()); i++)
  //{
  //  cout << "reaction["<<i<<"]: "<< these_rates[i] << endl;
  //}
                             
                                                                                                             
  //cout << "got reaction rates " << endl;
  // now print the map data to the vtk file
  vtk_print_cell_from_map_of_vectors(vtk_cell_out, solute_conc);
  
  vtk_print_cell_from_map_of_vectors(vtk_cell_out, rxn_rates);
  
  
}
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function bundles several of the vtk printing functions so
// you get the correct header, an appropriate filename, and both solute and
// mineral data
//
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::vtk_cell_bundler(double t_ime, int reference_switch, 
                                   string vtk_cell_fname, CRUNCH_engine& Ceng, 
                                   CRN_tParticle_bins& CRN_tPb)
{ 
  // set up the filename
  string time_bit = itoa( int(t_ime+0.5) );
	string vtk_ext = ".vtk";
	string fname = vtk_cell_fname+time_bit+vtk_ext;
	ofstream vtk_cell_out;
	vtk_cell_out.open(fname.c_str());
	
	// get the mass fractions: we do this here since it is only used for printing 
	// so this saves computational expense. 
	populate_cells_with_mfracs_and_depletion(CRN_tPb);
  
  // print to file
  vtk_print_cell_header(reference_switch, vtk_cell_out);
  vtk_print_cell_mineral_solid_state(vtk_cell_out);
  vtk_print_cell_CRUNCH_data(vtk_cell_out, Ceng); 
  vtk_print_cell_mineral_mfracs_and_depletion(vtk_cell_out);
  
  vtk_cell_out.close(); 
}

//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function takes the vec list vec generated from the 
// populate_cells_with_geochemical_data_from_CRNtPb
// and generates the crunch in files
//
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::generate_CRUNCH_in_files(CRUNCH_engine& Ceng,
                                   flowtube& ft, 
                                   CRN_tParticle_bins& CRN_tPb)
{

  int cells_in_bin = n_pdz_cells_per_bin+n_caz_cells_per_bin;

	vector<double> d_top_locs;
	vector<double> d_bottom_locs;
	
  // loop through each bin, generating a CrunchFlow in file
  for(int bn = 0; bn<n_bins; bn++)
  {
    // get the top and bottom depths in this bin
    CRN_tPb.partition_bins_into_cells(bn, ft,n_pdz_cells_per_bin, 
										n_caz_cells_per_bin,bottom_depth,d_top_locs,d_bottom_locs);

	  // get the pH vector (note: this is a stand in: will get it from the 
    // parsed data files later. SMM 10/08/2014
 	  vector<double> pH_vec;
	  double A = 0.2083;			// fitted from kate's data
	  double B = 6.2221;			//
	  pH_vec = Ceng.set_up_pH_for_particle(d_top_locs,d_bottom_locs, A, B);

	  int n_ts = 1;
	  list < vector<double> > default_concentrations 
            = Ceng.get_default_concentrations(n_ts,d_top_locs,d_bottom_locs);

    // create the in files. 
    // the first time you do this you need to have some default concentrations. 
    // after that you use the concentrations from the last run. 
    Ceng.create_CRUNCH_in_file(cells_in_bin, bn, cells_in_bin, pH_vec,
						d_top_locs, d_bottom_locs, default_concentrations, 
            vec_mineral_vpercents_old[bn], vec_mineral_ssa_old[bn]);
  }

}
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function calls CRUNCH for each bin and gets the data from the resulting
// files
//
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::call_CRUNCH_and_parse_data(CRUNCH_engine& Ceng)
{

  vector< list< vector<double> > > reset_vlv(n_bins);
  vec_new_conc = reset_vlv;
  vec_mineral_vpercents_new = reset_vlv;
  vec_new_min_ssa = reset_vlv;
  vec_new_rxn_rates = reset_vlv;


  // some vectors that are replaced during the parsing process
  vector<double> new_pH_vec;
  vector<double> spacings;
  vector<double> CRUNCH_tdepths;
  vector<double> CRUNCH_bdepths;  
  
  // some listvecs that will be replaced
  list< vector<double> > new_conc;
  list< vector<double> > mineral_vperc_new;
  list< vector<double> > new_min_ssa;
  list< vector<double> > new_rxn_rates;

  // call crunch for each bin
  for(int bn = 0; bn<n_bins; bn++)
  {
    // call crunch
    
    cout << "LINE 472, calling CrunchFlow in bin " << bn << endl;
    Ceng.call_CRUNCH(bn);
    
    // for bug checking, move the crunch files with bin number names
    int n_ts = 2;
    Ceng.move_CRUNCH_output_files_with_bin_number(n_ts,bn);
    //cout << "LINE 477, moved output files " << bn << endl;
        
    // now read in the resulting data
    int number_timestep = 1;
    int n_cells = n_pdz_cells_per_bin+n_caz_cells_per_bin;
    
    //cout << "LINE 483 Parsing the crunch data" << endl;
    // parse the resulting data. 
    Ceng.parse_CRUNCH_files(number_timestep, n_cells, bn, n_bins, n_cells,
						new_pH_vec, spacings, CRUNCH_tdepths, CRUNCH_bdepths,
						new_conc, mineral_vperc_new, new_min_ssa, new_rxn_rates);
		 //cout << "LINE 488 Parsed the crunch data, moving on to updating the vlvs" << endl;
    				
     // add the listvecs to the vlv
     vec_new_conc[bn] = new_conc;
     vec_mineral_vpercents_new[bn] = mineral_vperc_new;
     vec_new_min_ssa[bn] = new_min_ssa;
     vec_new_rxn_rates[bn] = new_rxn_rates;

  }
}
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This weathers particles in each bin
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::weather_particles_based_on_CRUNCH(CRN_tParticle_bins& CRN_tPb)
{
  // loop through the bins
  for (int bn = 0; bn<n_bins; bn++)
  {
	  CRN_tPb.weather_particles_from_CRUNCH(bn,n_pdz_cells_per_bin,
                    n_caz_cells_per_bin,bottom_depth,
										cell_data_map["verts_s"], cell_data_map["verts_z"], 
                    cell_data_map["verts_d"], cell_index_map["cell_node1"],
										cell_index_map["cell_node2"], cell_index_map["cell_node3"], 
                    cell_index_map["cell_node4"], vpi, vec_mineral_vpercents_old[bn],
										vec_mineral_vpercents_new[bn], vec_mineral_surface_area_old[bn],
										vec_mineral_mass_old[bn]);    
  }

}
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This function bundles all functions related to running CRUNCHflow 
//
//
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRUNCH_bins::run_CRUNCH_timestep(CRN_tParticle_bins& CRN_tPb, 
                  flowtube& ft, CRUNCH_engine& Ceng)
{
  // get the geochemical data from the particles
  populate_cells_with_geochemical_data_from_CRNtPb(ft, CRN_tPb);

  // generate the .in files for CRUNCH
  generate_CRUNCH_in_files(Ceng, ft, CRN_tPb);
  
  // call and parse CRUNCH output
  call_CRUNCH_and_parse_data(Ceng);

  // weather the particles
  weather_particles_based_on_CRUNCH(CRN_tPb);
}
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                  
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses the crunch output. The organisation of this is stupid and needs a 
// rewrite, but at the moment this isn't possible. 
// So the interface with CRN_tParticle_bins produces list vectors where the vectors
// have the size of total bins in the profile, but only the elements in the
// bin are ever used. So they are full of useless data. 
// 
// Crunch parsing, as written in CRUNCH_engine, doesn't know about bins. 
// So we need to map the appropriate bins to the crunch engine. 
//
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
map< string, vector<double> > CRUNCH_bins::parse_CRUNCH_vec_list_vec_to_vec_map(string master_name, 
                            list<string> element_list,
                            vector< list < vector<double> > >& vlv)
{
  // this is the map that stores the data
  map<string, vector<double> > data_map;
  
  
  if (int(vlv.size()) != n_bins)
  {
    cout << "Parsing vlv to map, your vlv size does not correspont to the number "
         << "of bins" << endl; 
  }
  else
  {
    // get the numer of elements in the lists
    int n_names = int(element_list.size());
    int n_elements_in_list = int(vlv[0].size());
    
    list< string >::iterator str_iter;	// list iterator for the string list
       
    // make sure the naming vector has the right number of elements
    string uscore = "_";   
    if (n_names != n_elements_in_list)
    {
      cout << "Parsing vlv to map, list names don't seem to correspond to " 
           << "elements in list." << endl;
      cout << "creating a stand in list" << endl;
      list<string> number_list;
      for(int i = 0; i<n_elements_in_list; i++)
      {
        string number_for_list = itoa(i);
        number_for_list =  master_name+uscore+number_for_list;
        number_list.push_back( number_for_list );
      }     
      element_list = number_list;
    }
    else
    {
    
      str_iter = element_list.begin();
      while(str_iter != element_list.end())
      {
        string this_name = master_name+uscore+*str_iter;
        *str_iter = this_name;
        str_iter++;
      }      
    }
    
    // now you need to populate the map with empty vectors
    
    vector<double> empty_vec;
    str_iter = element_list.begin();
    while(str_iter != element_list.end())
    {
      //cout << "CRUNCH_bins, adding key: " << *str_iter << endl;
      data_map[ *str_iter ] =  empty_vec;
      str_iter++;
    }
       
    // now go through the lvl, appending the vectors
	  list< vector<double> >::iterator vec_iter;	// list iterator for the vector
	  vector<double>::iterator element_iter_start;      // this is used to modify 
                                               //elements of the resulting vectors
    vector<double>::iterator element_iter_end; 
    
    // loop through the bins    
    for (int bn = 0; bn < n_bins; bn++ )
    {
    
      // get some information about the number of cells
    	int starting_cell = bn*(n_pdz_cells_per_bin+n_caz_cells_per_bin);

      // now loop through the elements
      vec_iter = vlv[bn].begin();
      str_iter = element_list.begin();
      while(vec_iter != vlv[bn].end())
      {
        // append the vector to the vector in the map
        vector<double> this_vec = *vec_iter;
        element_iter_start =  this_vec.begin();
        element_iter_end = this_vec.end();
        vector<double> vec_with_correct_cells;
        vec_with_correct_cells.assign(element_iter_start,element_iter_end);              
        //cout << "bn: " << bn << " mineral: " <<  *str_iter << " vector size: " 
        //     << vec_with_correct_cells.size() << endl;
        vector<double> map_vec = data_map[ *str_iter ];
        map_vec.insert(map_vec.end(), vec_with_correct_cells.begin(), 
                              vec_with_correct_cells.end());
        data_map[ *str_iter ] = map_vec;
        
        // increment the iterators
        vec_iter++;
        str_iter++;
      } 
    }                         
  }

  // check the contents of the map
  map< string, vector<double> >::iterator map_iter;
  
  map_iter = data_map.begin();
  while(map_iter != data_map.end())
  {
    vector<double> thisvec = map_iter->second;
    //cout << "key is: " << map_iter->first 
    //     << " and size is " <<  int(thisvec.size())  
    //     << " and n total cells are: " << total_cells << endl;  
    map_iter++;
  }

  return data_map;
}                            

//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function takes a vector list vector and converts its contents
// to a map of vectors. Each vector contains data from the lists
// with elements correspond to a cell.
//
// This function is mostly used with the vtk printing so you can see all the
// elements that are printed to the vtk files.
//
// The rationale for using a map is to be able to print to 
// the vtk file with a name of a field 
// 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
map< string, vector<double> > CRUNCH_bins::parse_vec_list_vec_to_vec_map(string master_name, 
                            list<string> element_list,
                            vector< list < vector<double> > >& vlv)
{
  // this is the map that stores the data
  map<string, vector<double> > data_map;
  
  
  if (int(vlv.size()) != n_bins)
  {
    cout << "Parsing vlv to map, your vlv size does not correspont to the number "
         << "of bins" << endl; 
  }
  else
  {
    // get the numer of elements in the lists
    int n_names = int(element_list.size());
    int n_elements_in_list = int(vlv[0].size());
    
    list< string >::iterator str_iter;	// list iterator for the string list
       
    // make sure the naming vector has the right number of elements
    string uscore = "_";   
    if (n_names != n_elements_in_list)
    {
      cout << "Parsing vlv to map, list names don't seem to correspond to " 
           << "elements in list." << endl;
      cout << "creating a stand in list" << endl;
      list<string> number_list;
      for(int i = 0; i<n_elements_in_list; i++)
      {
        string number_for_list = itoa(i);
        number_for_list =  master_name+uscore+number_for_list;
        number_list.push_back( number_for_list );
      }     
      element_list = number_list;
    }
    else
    {
    
      str_iter = element_list.begin();
      while(str_iter != element_list.end())
      {
        string this_name = master_name+uscore+*str_iter;
        *str_iter = this_name;
        str_iter++;
      }      
    }
    
    // now you need to populate the map with empty vectors
    
    vector<double> empty_vec;
    str_iter = element_list.begin();
    while(str_iter != element_list.end())
    {
      //cout << "CRUNCH_bins, adding key: " << *str_iter << endl;
      data_map[ *str_iter ] =  empty_vec;
      str_iter++;
    }
       
    // now go through the lvl, appending the vectors
	  list< vector<double> >::iterator vec_iter;	// list iterator for the vector
	  vector<double>::iterator element_iter_start;      // this is used to modify 
                                               //elements of the resulting vectors
    vector<double>::iterator element_iter_end; 
    
    // loop through the bins    
    for (int bn = 0; bn < n_bins; bn++ )
    {
    
      // get some information about the number of cells
    	int starting_cell = bn*(n_pdz_cells_per_bin+n_caz_cells_per_bin);

      // now loop through the elements
      vec_iter = vlv[bn].begin();
      str_iter = element_list.begin();
      while(vec_iter != vlv[bn].end())
      {
        // append the vector to the vector in the map
        vector<double> this_vec = *vec_iter;
        element_iter_start =  this_vec.begin()+starting_cell;
        element_iter_end = this_vec.begin()+starting_cell+n_pdz_cells_per_bin+n_caz_cells_per_bin;
        vector<double> vec_with_correct_cells;
        vec_with_correct_cells.assign (element_iter_start,element_iter_end);
               
        //cout << "bn: " << bn << " mineral: " <<  *str_iter << " vector size: " 
        //     << vec_with_correct_cells.size() << endl;
        vector<double> map_vec = data_map[ *str_iter ];
        map_vec.insert(map_vec.end(), vec_with_correct_cells.begin(), 
                              vec_with_correct_cells.end());
        data_map[ *str_iter ] = map_vec;
        
        // increment the iterators
        vec_iter++;
        str_iter++;
      } 
    }                         
  }

  // check the contents of the map
  map< string, vector<double> >::iterator map_iter;
  
  map_iter = data_map.begin();
  while(map_iter != data_map.end())
  {
    vector<double> thisvec = map_iter->second;
    //cout << "key is: " << map_iter->first 
    //     << " and size is " <<  int(thisvec.size())  
    //     << " and n total cells are: " << total_cells << endl;  
    map_iter++;
  }

  return data_map;
}                            
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This returns a list of the particle names
//  used for printing to vtk
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
list<string> CRUNCH_bins::get_names_of_minerals()
{
  int n_types = vpi.get_n_types();
  
  list<string> names_of_minerals;
  
  for(int i = 0; i<n_types; i++)
  {
    //cout << "Line 463, mineral name is: " <<  vpi.get_type_name(i) << endl;
    names_of_minerals.push_back( vpi.get_type_name(i));
  }
  return names_of_minerals;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This returns a list of the primary species names
//  used for printing to vtk
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
list<string> CRUNCH_bins::get_names_of_primary_species(CRUNCH_engine& cEng)
{

  // get the list. 
  list<string> names_of_pspecies = cEng.get_primary_species_names();
  
  // loop through this list removing any control characters
  list<string>::iterator liter;
  liter = names_of_pspecies.begin();
  while(liter!=names_of_pspecies.end())
  {
    string this_string = *liter;
    
    int len =  this_string.length();
    if(len != 0)
    {
      if (iscntrl(this_string[len-1]))
      {
        cout << "the last item in the species is a control character!" << endl;
        this_string.erase(len-1);
      }
    }
    else
    {
      cout << "Warning, getting species from CRUNCH_bins, but species "
           << "list contains and empty string." << endl;
    } 
    *liter = this_string;
    liter++;           
  }
    
  
  // note, the first species should always be H+. 
  // check this
  liter = names_of_pspecies.begin();
  string this_thing = *liter;
  string Hplus = "H+";
  if (this_thing.compare(Hplus) != 0)
  {
    cout << "WARNING!!, first species should be " << Hplus
         << ", but instead it is: " << this_thing << " this_thing" << endl;
  }
  
  // now we get rid of the H+ name since the H+ is reported in the pH 
  // parameter of CRUNCHflow
  //names_of_pspecies.pop_front();
  
  // now print the primary species for bug checking
  //liter = names_of_pspecies.begin();
  //while (liter != names_of_pspecies.end())
  //{
  //  cout << "Line 713, species name is: " <<  *liter << endl;
  //  liter++;
  //}
  return names_of_pspecies;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
