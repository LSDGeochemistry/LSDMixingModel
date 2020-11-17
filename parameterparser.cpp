//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//parameter parser
// An object that holds properties of particles
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group mixing model
//  for exploring hillslope mixing and particle weathering
//
// Developed by:
//  Simon M. Mudd
//
// Copyright (C) 2018 Simon M. Mudd 2018
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <fstream>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include "TNT/tnt.h"
#include "LSDStatsTools.hpp"
#include "parameterparser.hpp"
using namespace std;
using namespace TNT;

#ifndef parameterparser_CPP
#define parameterparser_CPP

// This basic function is not part of the parameter parser object
// but instead is the function called at the beginning of any
// LSDTT driver call
vector<string> DriverIngestor(int nNumberofArgs,char *argv[])
{
  cout << "=========================================================" << endl;
  cout << "|| You have called an LSDMixingModel program.            ||" << endl;
  cout << "|| Prepare to explore hillslope data!                ||" << endl;
  cout << "=========================================================" << endl;

  string path_name = ".";
  string file_name = "LSDMM_parameters.driver";

  //Test for correct input arguments
  if (nNumberofArgs == 1)
  {
    cout << "You have not given me any arguments. I am going to look" << endl;
    cout << "in this directory for a file with the extension .driver" << endl;
    cout << "I'll use the first one I find (in alphabetical ordering)." << endl;
    cout << "If I don't find one I am going to exit." << endl;
  }
  if (nNumberofArgs == 2)
  {
    cout << "I have one argument. I don't know if this is a directory path" << endl;
    cout << "or a driver filename. I am going to assume it is a directory path" << endl;
    cout << "if it containes the character . or /" << endl;

    string temp_arg = argv[1];
    string s_dot = ".";
    string sl_dot = "/";

    bool this_is_a_path = false;
    bool path_find = false;
    if (temp_arg.find(s_dot) != std::string::npos)
    {
      path_find = true;
    }
    if (temp_arg.find(sl_dot) != std::string::npos)
    {
      path_find = true;
    }

    if (this_is_a_path)
    {
      path_name = temp_arg;
    }
    else
    {
      file_name = temp_arg;
    }
  }
  //Test for correct input arguments
  if (nNumberofArgs==3)
  {
    cout << "I am reading the two arguments you gave me as the path name and the file name." << endl;
    path_name = argv[1];
    file_name = argv[2];
  }
  if (nNumberofArgs>=3)
  {
    cout << "You have provided more than two arguments. " << endl;
    cout << "I only expect 2. I am going to assume you meant" << endl;
    cout << "to give me the first two." << endl;
    path_name = argv[1];
    file_name = argv[2];

  }

  vector<string> path_and_file;
  path_and_file.push_back(path_name);
  path_and_file.push_back(file_name);

  return path_and_file;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Create functions
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void parameterparser::create()
{
  cout << "I have created an empty parameter parser object. " << endl;
  cout << "Surely you want to give it a filename?" << endl;
}

// This creates using a path and a filename
void parameterparser::create(string PathName, string FileName)
{

  // Make sure the path has an extension
  PathName = FixPath(PathName);
  string FullName = PathName+FileName;

  param_file_path = PathName;
  param_fname = FileName;

  ifstream file_info_in;
  file_info_in.open(FullName.c_str());

  // check if the parameter file exists
  if( file_info_in.fail() )
  {
    cout << "\nFATAL ERROR: The parameter file \"" << FullName
         << "\" doesn't exist" << endl;
    exit(EXIT_FAILURE);
  }

  // now ingest the parameters
  cout << "Parsing the file" << endl;
  LSDPP_parse_file_into_parameter_map(FullName);
  parse_file_IO();

  // make sure the files are okay
  check_boundary_conditions();
  check_file_extensions_and_paths();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Gets a line of the parameter file. Has a long buffer so you can add long path
// names.
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void parameterparser::LSDPP_parse_line(ifstream &infile, string &parameter, string &value)
{
  char c;
  char buff[1024];
  int pos = 0;
  int word = 0;

  while ( infile.get(c) )
  {
    if (pos >= 1024)
    {
      cout << "Buffer overrun, word too long in parameter line: " << endl;
      string line;
      getline(infile, line);
      cout << "\t" << buff << " ! \n" << line << endl;
      exit(1);
    }
    // preceeding whitespace
    if (c == '#')
    {
      if (word == 0)
      {
        parameter = "NULL";
        value = "NULL";
      }
      if (word == 1)
        value = "NULL";
      word = 2;
    }

    if ((c == ' ' || c == '\t') && pos == 0)
      continue;
    else if ( (c == ':' && word == 0) || ( (c == ' ' || c == '\n' || c == '\t') && word == 1))
    {
      while (buff[pos-1] == ' ' || buff[pos-1] == '\t')
        --pos;    // Trailing whitespace
      buff[pos] = '\0';  // Append Null char
      if (word == 0)
        parameter = buff;  // Assign buffer contents
      else if (word == 1)
        value = buff;
      ++word;
      pos = 0;    // Rewind buffer
    }
    else if ( c == '\n' && word == 0 )
    {
      parameter = "NULL";
      buff[pos] = '\0';
      value = buff;
      ++word;
    }
    else if (word < 2)
    {
      buff[pos] = c;
      ++pos;
    }

    if (c == '\n')
      break;
  }
  if (word == 0)
  {
    parameter = "NULL";
    value = "NULL";
  }
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This reads the parameter file, placing all parameters into a map
// with string values and string key. As you give the parameter parser
// default maps, it will scan these sting and convert them into the correct data type
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void parameterparser::LSDPP_parse_file_into_parameter_map(string FullName)
{
  ifstream infile;
  infile.open(FullName.c_str());
  string parameter, value, lower, lower_val;
  string bc;

  cout << "Hello, I am going to parse your LSDTopoTools parameter file for you. " << endl;
  cout << "The parameter filename is: " << FullName << endl;

  // this will hold all the parameter values.
  map<string,string> temp_parameters;

  // now ingest parameters
  while (infile.good())
  {
    LSDPP_parse_line(infile, parameter, value);
    lower = parameter;
    //if (parameter == "NULL")
    //  continue;
    //for (unsigned int i=0; i<parameter.length(); ++i)
    //{
    //  lower[i] = tolower(parameter[i]);
    //}

    cout << "parameter is: " << lower << " and value is: " << value << endl;

    // get rid of control characters
    value = RemoveControlCharactersFromEndOfString(value);

    temp_parameters[lower] = value;
  }

  parameter_map = temp_parameters;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This uses the parameter map to get file input and output
void parameterparser::parse_file_IO()
{
  cout << endl << endl << endl << "----------------------" << endl;
  cout << "Parsing the file I/O" << endl;

  if(parameter_map.find("write path") != parameter_map.end())
  {
    write_path = parameter_map["write path"];
    // get rid of any control characters from the end (if param file was made in DOS)
    write_path = RemoveControlCharactersFromEndOfString(write_path);
  }
  if(parameter_map.find("write fname") != parameter_map.end())
  {
    write_fname = parameter_map["write fname"];
    // get rid of any control characters from the end (if param file was made in DOS)
    write_fname = RemoveControlCharactersFromEndOfString(write_fname);
    //cout << "Got the write name, it is: "  << write_fname << endl;
  }
  if(parameter_map.find("read path") != parameter_map.end())
  {
    read_path = parameter_map["read path"];
    // get rid of any control characters from the end (if param file was made in DOS)
    read_path = RemoveControlCharactersFromEndOfString(read_path);
    //cout << "Got the write name, it is: "  << write_fname << endl;
  }
  else
  {
    cout << "I did not find a read path so I am assuming the file is in this current directory." << endl;
    read_path = "./";
  }

  if(parameter_map.find("read fname") != parameter_map.end())
  {
    read_fname = parameter_map["read fname"];
    // get rid of any control characters from the end (if param file was made in DOS)
    read_fname = RemoveControlCharactersFromEndOfString(read_fname);
    //cout << "Got the write name, it is: "  << write_fname << endl;
  }


}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This forces parsing of everything in the file
// Used for copying parameter files
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void parameterparser::force_parse()
{
  cout << "Forcing parsing of the parameter file" << endl;
  for( map<string, string >::iterator it = parameter_map.begin(); it != parameter_map.end(); ++it)
  {
    string key = it->first;
    cout << "Key is: " <<it->first << "\n";
    if(key != "CHeads file" && key != "read fname" && key != "write fname" &&
        key != "read path" && key != "write path" && key != "read extension" &&
        key != "write extension" && key != "NULL")
    {
      parameters_read_map[it->first] = it->second;
    }
  }

}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses all the default parameter maps.
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void parameterparser::parse_all_parameters(map<string,float> default_map_f,
                      map<string,int> default_map_i, map<string,bool> default_map_b,
                      map<string,string> default_map_s)
{
  parse_float_parameters(default_map_f);
  parse_int_parameters(default_map_i);
  parse_bool_parameters(default_map_b);
  parse_string_parameters(default_map_s);

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses all the default parameter maps.
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void parameterparser::parse_all_parameters(map<string,float> default_map_f,
                      map<string,int> default_map_i, map<string,bool> default_map_b,
                      map<string,string> default_map_s,
                      map<string,double> default_map_d)
{
  parse_float_parameters(default_map_f);
  parse_int_parameters(default_map_i);
  parse_bool_parameters(default_map_b);
  parse_string_parameters(default_map_s);
  parse_double_parameters(default_map_d);

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// These two functions takes a map of defualt parameters and returns the parameters for the
// current implementation
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void parameterparser::parse_float_parameters(map<string,float> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      cout << "Found key " << these_keys[i];

      // convert the value to float
      float_parameters[these_keys[i]] = atof(parameter_map[these_keys[i]].c_str());
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
      cout << " it is: " << parameter_map[these_keys[i]] << " check: " << float_parameters[these_keys[i]] << endl;

    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      float_parameters[these_keys[i]] = default_map[these_keys[i]];
      defaults_used_map[these_keys[i]] = dtoa(default_map[these_keys[i]]);

    }

  }
}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Parse double parameters (for coords) FJC 20/11/17
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void parameterparser::parse_double_parameters(map<string,double> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      cout << "Found key " << these_keys[i];

      // convert the value to float
      double_parameters[these_keys[i]] = atof(parameter_map[these_keys[i]].c_str());
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
      cout << " it is: " << parameter_map[these_keys[i]] << " check: " << double_parameters[these_keys[i]] << endl;

    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      double_parameters[these_keys[i]] = default_map[these_keys[i]];
      defaults_used_map[these_keys[i]] = dtoa(default_map[these_keys[i]]);

    }

  }
}

void parameterparser::parse_int_parameters(map<string,int> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      //cout << "Found the key" << endl;
      // convert the value to float
      int_parameters[these_keys[i]] = atoi(parameter_map[these_keys[i]].c_str());
      //cout << " it is: " << parameter_map[these_keys[i]] << " check: " << int_parameters[these_keys[i]] << endl;
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      int_parameters[these_keys[i]] = default_map[these_keys[i]];
      defaults_used_map[these_keys[i]] = itoa(default_map[these_keys[i]]);
    }
  }
}


void parameterparser::parse_bool_parameters(map<string,bool> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      // convert the value to bool
      string value = parameter_map[these_keys[i]];
      bool temp_bool = (value == "true" || value== "True" || value == "TRUE" || value== "T" || value== "t") ? true : false;
      bool_parameters[these_keys[i]] = temp_bool;
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      bool_parameters[these_keys[i]] = default_map[these_keys[i]];
      if (default_map[these_keys[i]] == true)
      {
        defaults_used_map[these_keys[i]] = "true";
      }
      else
      {
        defaults_used_map[these_keys[i]] = "false" ;
      }

    }
  }
}

void parameterparser::parse_string_parameters(map<string,string> default_map)
{
  // the idea is to look through the default map, getting the keys, and then
  // looking for the keys in the parameter maps
  vector<string> these_keys = extract_keys(default_map);

  // loop through the keys
  int n_keys = int(these_keys.size());
  for(int i = 0; i<n_keys; i++)
  {
    cout << "Key is: " << these_keys[i] << endl;

    // If the key is contained in the parsed parameters, use the parsed parameter
    if(parameter_map.find(these_keys[i]) != parameter_map.end())
    {
      string_parameters[these_keys[i]] = parameter_map[these_keys[i]];
      parameters_read_map[these_keys[i]] = parameter_map[these_keys[i]];
    }
    else  // the key is not in the parsed parameters. Use the default.
    {
      string_parameters[these_keys[i]] = default_map[these_keys[i]];
      defaults_used_map[these_keys[i]] = default_map[these_keys[i]];
    }
  }
}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses a vector of strings
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<string> parameterparser::parse_string_vector(string key)
{
  string this_string = string_parameters[key];

  // reset the string vec
  vector<string> this_string_vec;

  // create a stringstream
  stringstream ss(this_string);

  // import the data, using a comma to separate
  while( ss.good() )
  {
    string substr;
    getline( ss, substr, ',' );

    // remove the spaces
    substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

    // remove control characters
    substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

    // add the string to the string vec
    this_string_vec.push_back( substr );
  }

  return this_string_vec;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses a vector of ints
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<int> parameterparser::parse_int_vector(string key)
{
  string this_string = string_parameters[key];

  // reset the string vec
  vector<int> this_int_vec;

  // create a stringstream
  stringstream ss(this_string);

  // import the data, using a comma to separate
  while( ss.good() )
  {
    string substr;
    getline( ss, substr, ',' );

    // remove the spaces
    substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

    // remove control characters
    substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

    // add the string to the string vec
    this_int_vec.push_back( atoi(substr.c_str()) );
  }

  return this_int_vec;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This parses a vector of ints
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
vector<float> parameterparser::parse_float_vector(string key)
{
  cout << "I am going to parse a float vector for you!" << endl;

  string this_string = string_parameters[key];

  // reset the string vec
  vector<float> this_float_vec;

  // create a stringstream
  stringstream ss(this_string);

  // import the data, using a comma to separate
  while( ss.good() )
  {
    string substr;
    getline( ss, substr, ',' );

    // remove the spaces
    substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

    // remove control characters
    substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

    // add the string to the string vec
    this_float_vec.push_back( atof(substr.c_str()) );
  }

  return this_float_vec;

}


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function checks filenames to see if they include a path. 
// If not it add the read path
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
string parameterparser::check_for_path_and_add_read_path_if_required(string this_string)
{
  string sl_dot = "/";
  string new_string;
  
  cout << "The string to check is: " << this_string << endl;
  if (this_string == "NULL" || this_string == "null" || this_string == "Null")
  {
    new_string = "NULL";
  }
  else
  {
    if (this_string.find(sl_dot) != std::string::npos)
    {
      cout << "This filename includes a path. I am not going to modify it." << endl;
      new_string = this_string;
    } 
    else
    {
      cout << "This finlename doesn't have a path. I am adding the read path." << endl;
      new_string = read_path+this_string;
      cout << "The new filename is: " << new_string << endl;
    }
  } 

  return new_string;

}

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This function strips the text after the final dot in a string
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
string parameterparser::get_string_before_dot(string this_string)
{
  string cut_string;
  unsigned found = this_string.find_last_of(".");
  cut_string = this_string.substr(0,found);
  return cut_string;
}
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This prints parameters read to file, so you can make sure your parameters
// have ingested properly
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void parameterparser::print_parameters()
{
  string fname = write_path+write_fname+"_ingestedParam.param";
  ofstream params_out;
  params_out.open(fname.c_str());

  params_out << "# Here are the paramters ingested and set by the parameter parser" << endl;
  params_out << "# The file names and paths are: " << endl;
  params_out << "read path: " << read_path << endl;
  params_out << "read fname: " << read_fname << endl;
  params_out << "write path: " << write_path << endl;
  params_out << "write fname: " << write_fname << endl;
  
  params_out << "# ===================================="  << endl;
  params_out << "# Now for parameters read from file." << endl;
  params_out << "# If an expected parameter is not here check your spelling." << endl;

  vector<string> empty_vec;
  vector<string> these_keys = extract_keys(parameters_read_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    params_out << these_keys[i] << ": " << parameters_read_map[these_keys[i]]  << endl;
  }


  params_out << endl << "# ===================================="  << endl;
  params_out << "# Now for the default parameters." << endl;

  these_keys = empty_vec;
  these_keys = extract_keys(defaults_used_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    params_out << these_keys[i] << ": " << defaults_used_map[these_keys[i]]  << endl;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This prints parameters read to file, so you can make sure your parameters
// have ingested properly
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void parameterparser::print_parameters(string fname_prefix)
{
  string fname = write_path+fname_prefix;
  ofstream params_out;
  params_out.open(fname.c_str());

  params_out << "# Here are the paramters ingested and set by the parameter parser" << endl;
  params_out << "# The file names and paths are: " << endl;
  params_out << "read path: " << read_path << endl;
  params_out << "read fname: " << read_fname << endl;
  params_out << "write path: " << write_path << endl;
  params_out << "write fname: " << write_fname << endl;
  

  params_out << "# ===================================="  << endl;
  params_out << "# Now for parameters read from file." << endl;
  params_out << "# If an expected parameter is not here check your spelling." << endl;

  vector<string> empty_vec;
  vector<string> these_keys = extract_keys(parameters_read_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    params_out << these_keys[i] << ": " << parameters_read_map[these_keys[i]]  << endl;
  }


  params_out << endl << "# ===================================="  << endl;
  params_out << "# Now for the default parameters." << endl;

  these_keys = empty_vec;
  these_keys = extract_keys(defaults_used_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    params_out << these_keys[i] << ": " << defaults_used_map[these_keys[i]]  << endl;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This prints parameters read to file, so you can make sure your parameters
// have ingested properly
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void parameterparser::replace_and_print_parameter_file(string parameter_fname,
                                     string new_read_path, string new_read_fname,
                                     string new_write_path, string new_write_fname,
                                     map<string,string> replace_parameters)
{
  string fname = write_path+parameter_fname;
  ofstream params_out;
  params_out.open(fname.c_str());

  params_out << "# This is an adjusted parameter file" << endl;
  params_out << "# The file names and paths are: " << endl;
  params_out << "read path: " << new_read_path << endl;
  params_out << "read fname: " << new_read_fname << endl;
  params_out << "write path: " << new_write_path << endl;
  params_out << "write fname: " << new_write_fname << endl;
  

  params_out << "# ===================================="  << endl;
  params_out << "# Now for parameters read from file." << endl;
  params_out << "# If an expected parameter is not here check your spelling." << endl;

  vector<string> empty_vec;
  vector<string> these_keys = extract_keys(parameters_read_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    // If the key is contained in the replace parameters, use the replace parameter
    if(replace_parameters.find(these_keys[i]) != replace_parameters.end())
    {
      cout << "I found a replace parameter!" << endl;
      params_out << these_keys[i] << ": " << replace_parameters[these_keys[i]]  << endl;
    }
    else
    {
      params_out << these_keys[i] << ": " << parameters_read_map[these_keys[i]]  << endl;
    }
  }

  params_out << endl << "# ===================================="  << endl;
  params_out << "# Now for the default parameters." << endl;

  these_keys = empty_vec;
  these_keys = extract_keys(defaults_used_map);
  for(int i = 0; i< int(these_keys.size()); i++)
  {
    // If the key is contained in the replace parameters, use the replace parameter
    if(replace_parameters.find(these_keys[i]) != replace_parameters.end())
    {
      cout << "I found a replace parameter!" << endl;
      params_out << these_keys[i] << ": " << replace_parameters[these_keys[i]]  << endl;
    }
    else
    {
      params_out << these_keys[i] << ": " << defaults_used_map[these_keys[i]]  << endl;
    }
  }
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#endif
