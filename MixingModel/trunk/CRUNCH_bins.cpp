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
    
      cout << "bin: " << bin << " cell: " << cell << endl;
      cout << "s:\t" << verts_s[ci1] << "\t" << verts_s[ci2] << "\t" 
                     << verts_s[ci3] << "\t" << verts_s[ci4] << endl; 
      cout << "z:\t" << verts_z[ci1] << "\t" << verts_z[ci2] << "\t" 
                     << verts_z[ci3] << "\t" << verts_z[ci4] << endl;       
      cout << "d:\t" << verts_d[ci1] << "\t" << verts_d[ci2] << "\t" 
                     << verts_d[ci3] << "\t" << verts_d[ci4] << endl;                                       
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
// The create function:
// this just has some information about the cell locations
// the data is stored in vectors, you need to get the index into
// the vector of from the bin and cell number
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
	 vec_mineral_vfracs_old = empty_vlv;
 	 vec_mineral_ssa_old = empty_vlv;
 	 vec_mineral_mass_old = empty_vlv;
 	 vec_mineral_surface_area_old = empty_vlv;
   
  // get the locations of the cells and the cell indices. 
  // this function also updates the cell indices of the particles
  for (int bn = 0; bn<n_bins; bn++)
  {
    //cout << "size vmvo: " << vec_mineral_vfracs_old.size() << " and bn: " << bn << endl;
    // now collect the data from from the particles
    // this collects data from each cell in this bin
	  list< vector<double> > mineral_vfracs_old;
 	  list< vector<double> > mineral_ssa_old;
 	  list< vector<double> > mineral_mass_old;
 	  list< vector<double> > mineral_surface_area_old;
	  CRN_tPb.get_data_by_cell_volumetric_for_CRUNCH(bn,n_pdz_cells_per_bin, 
                n_caz_cells_per_bin,bottom_depth,
								verts_s, verts_z, verts_d,
								cell_node1, cell_node2, cell_node3, cell_node4, vpi,
								mineral_vfracs_old,mineral_ssa_old,
								mineral_surface_area_old,mineral_mass_old);   
    vec_mineral_vfracs_old[bn] = mineral_vfracs_old;
    vec_mineral_ssa_old[bn] = mineral_ssa_old;
    vec_mineral_mass_old[bn] = mineral_mass_old;
    vec_mineral_surface_area_old[bn] = mineral_surface_area_old;
                
    // each particle type has a vector, so we need to get each particle
    // type out.
                 
								
    // now you need to distribute data to the main vectors								
  
  
  }

    // first partition the bins into cells. This just gives the top and bottom
    // depths in the centrepoint of each cell. 
    //CRN_tPb.partition_bins_into_cells(bn, ft, n_PDZ_cells_per_bin, n_CAZ_cells_per_bin,
		//								bottom_depth,bin_d_top_locs,bin_d_bottom_locs);
										  
  
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
  //list< vector<double> > this_lv = vec_mineral_vfracs_old[0];
  //vector<double> = 
  
  // get the names of the minerals
  list<string> mineral_names = get_names_of_minerals();
  
  // get the maps
  string mvname = "Mineral_volume_fractions"; 
  map< string, vector<double> > mvfracs =  parse_vec_list_vec_to_vec_map(mvname, 
                            mineral_names, vec_mineral_vfracs_old);
  string mssaname = "Mineral_specific_surface_area_m2perg"; 
  map< string, vector<double> > mssafracs =  parse_vec_list_vec_to_vec_map(mssaname, 
                            mineral_names, vec_mineral_ssa_old);
  string mmname = "Mineral_mass_kg"; 
  map< string, vector<double> > mmfracs =  parse_vec_list_vec_to_vec_map(mmname, 
                            mineral_names, vec_mineral_mass_old);
  string msaname = "Mineral_surface_area_m2"; 
  map< string, vector<double> > msafracs =  parse_vec_list_vec_to_vec_map(msaname, 
                            mineral_names, vec_mineral_surface_area_old);                                                                                    

  // now print the map data to the vtk file
  vtk_print_cell_from_map_of_vectors(vtk_cell_out, mvfracs);
  
  
  

}
//==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


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
      cout << "CRUNCH_bins, adding key: " << *str_iter << endl;
      data_map[ *str_iter ] =  empty_vec;
      str_iter++;
    }
       
    // now go through the lvl, appending the vectors
	  list< vector<double> >::iterator vec_iter;	// list iterator for the vector
    
    // loop through the bins    
    for (int bn = 0; bn < n_bins; bn++ )
    {
      // now loop through the elements
      vec_iter = vlv[bn].begin();
      str_iter = element_list.begin();
      while(vec_iter != vlv[bn].end())
      {
        // append the vector to the vector in the map
        vector<double> this_vec = *vec_iter;
        vector<double> map_vec = data_map[ *str_iter ];
        map_vec.insert(map_vec.end(), this_vec.begin(), this_vec.end());
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
    cout << "key is: " << map_iter->first 
         << " and size is " <<  int(thisvec.size()) << endl;  
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
