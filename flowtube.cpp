//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// flowtub
// An object that controls a 1-D hillslope, keeps track of hillslope
// evolution, sediment transport and interfaces with particles
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

#include <iostream>
#include <string>
#include <math.h>
#include "flowtube.hpp"
#include "tParticle.hpp"
#include "CRN_tParticle_bins.hpp"
using namespace std;

// the create function. currently empty
void flowtube::create()
{
	cout << "you need to initialize with a parameter file!" << endl;
}

// create algorithm that takes an infile, it reads the
// tube properties and loads them into the data vactors
void flowtube::create(ifstream& paramfile_in)
{
	int n_box_bdrys;
	string s_test;

	// load the number of nodes
	paramfile_in >> s_test >> n_box_bdrys
			   >> s_test >> s_test >> s_test >> s_test >> s_test;
	n_nodes = (n_box_bdrys-1)/2;

	// temporary vectors, the data from the infile is loaded into
	// these vectors
	vector<double> temp_A(n_nodes);
	vector<double> distance_h(n_nodes);
	vector<double> distance_b(n_nodes+1);
	vector<double> temp_b(n_nodes+1);
	vector<double> meas_zeta(n_nodes);
	vector<double> meas_h(n_nodes);
	vector<double> meas_eta(n_nodes);
	vector<double> temp_DeltaXh(n_nodes-1);
	vector<double> temp_DeltaXb(n_nodes);
	vector<double> temp_A_bins( 2*n_nodes );
	vector<double> temp_bin_edge_loc( 2*n_nodes + 1);
	//cout << "flowtube.cpp LINE 83 size temp bel: " << temp_bin_edge_loc.size() << endl;


	// now load the data in from the infile
	double A1,A2,ds_dist1,ds_dist2,width1,width2,measz1,measz2,meas_h1,meas_h2;
	int counter = 0;
	int counter_bin = 0;
	double total_A = 0;
	for (int i = 0; i<n_nodes; i++)
	{
		paramfile_in >> ds_dist1 >> width1 >> A1 >> measz1 >> meas_h1
					 >> ds_dist2 >> width2 >> A2 >> measz2 >> meas_h2;

		//cout << "LINE 56 flowtube.cpp A1: " << A1 << " and A2: " << A2 << endl;

		temp_bin_edge_loc[counter_bin] = ds_dist1;
		temp_A_bins[counter_bin] = A1;
		counter_bin++;
		temp_bin_edge_loc[counter_bin] = ds_dist2;
		temp_A_bins[counter_bin] = A2;
		counter_bin++;

		temp_A[counter] = A1+A2;
		total_A += temp_A[counter];
		distance_h[counter] = ds_dist2;
		distance_b[counter] = ds_dist1;
		meas_zeta[counter] = measz2;
		meas_h[counter] = meas_h2;
		meas_eta[counter] = meas_zeta[counter]-meas_h[counter];
		temp_b[counter] = width1;
		//if (i > 0)
		// temp_b[counter-1] = width1;
		//cout << distance_b[counter] << endl;
		counter++;
		//cout << "flowtube.cpp LINE 83 size temp bel: " << temp_bin_edge_loc.size() << endl;
	}

	//cout << "LINE 84 flowtube.cpp" << endl;
	//for (int i =  0; i<2*n_nodes; i++)
	//{
	//	cout << "temp_A_bins [" << i << "]: " << temp_A_bins[i] << endl;
	//}

	paramfile_in >> ds_dist1 >> width1 >> A1 >> measz1 >> meas_h1;
	temp_bin_edge_loc[counter_bin] = ds_dist1;
	temp_b[counter] = width1;
	distance_b[counter] = ds_dist1;

	for (int i =0; i<n_nodes-1; i++)
	{
		temp_DeltaXh[i] = distance_h[i+1]-distance_h[i];
	}
	for (int i =0; i<n_nodes; i++)
	{
		temp_DeltaXb[i] = distance_b[i+1]-distance_b[i];
	}

	//cout << "flowtube.cpp LINE 83 size temp bel: " << temp_bin_edge_loc.size() << endl;

	// now save the data from the infile to the private data members
	A = temp_A;
	b = temp_b;
	A_bins = temp_A_bins;
	s_h = distance_h;
	s_b = distance_b;
	zeta = meas_zeta;
	h = meas_h;
	eta = meas_eta;
	DeltaXh = temp_DeltaXh;
	DeltaXb = temp_DeltaXb;
	bin_edge_loc = temp_bin_edge_loc;

	old_h = h;
	old_zeta = zeta;
	old_eta = eta;
	intermediate_h = h;
	intermediate_zeta = zeta;
	pre_surface_h = h;
	pre_surface_zeta = zeta;

	vector<double> fluff_temp (n_nodes,0.0);
	fluff = fluff_temp;
	vector<double> MF_temp(n_nodes+1,0.0);
	Mass_Flux = MF_temp;

	//cout << "flowtube.cpp LINE 100 size bel: " << bin_edge_loc.size() << endl;

//	cout << "LINE 120 flowtube.cpp" << endl;
//	for (int i =  0; i<2*n_nodes; i++)
//	{
//		cout << "zeta [" << i << "]: " << zeta[i] << endl;
//	}

}


void flowtube::create(int sn_nodes, double sS_c, double sK_h, double sW_0, double sgamma,
			 double srho_s, double srho_r, double sN, double sN_0, double sN_m,
			 double sbeta, double sK_g, vector<double> sA, vector<double> sA_bins,
			 vector<double> sb, vector<double> szeta, vector<double> seta,
			 vector<double> sh, vector<double> ss_h, vector<double> ss_b,
			 vector<double> sDeltaXh, vector<double> sDeltaXb,
			 vector<double> sbin_edge_loc,vector<double> sold_h,
			 vector<double> sold_eta, vector<double> sold_zeta,
			 vector<double> sintermediate_zeta, vector<double> sintermediate_h,
			 vector<double> spre_surface_zeta, vector<double> spre_surface_h,
			 vector<double> sMass_Flux, vector<double> sfluff)
{
	n_nodes = sn_nodes; S_c = sS_c; K_h = sK_h; W_0 = sW_0; gamma = sgamma;
	rho_s = srho_s; rho_r = srho_r; N = sN; N_0 = sN_0; N_m = sN_m; beta = sbeta;
	K_g = sK_g; A = sA; A_bins = sA_bins; b = sb; zeta = szeta; eta = seta;
	h = sh; s_h = ss_h; s_b = ss_b; DeltaXh = sDeltaXh; DeltaXb = sDeltaXb;
	bin_edge_loc = sbin_edge_loc; old_h = sold_h; old_eta = sold_eta; old_zeta = sold_zeta;
	intermediate_zeta = sintermediate_zeta; intermediate_h = sintermediate_h;
	pre_surface_zeta = spre_surface_zeta; pre_surface_h = spre_surface_h;
	Mass_Flux = sMass_Flux; fluff = sfluff;
}


flowtube& flowtube::operator=(flowtube& rhs)
 {
  if (&rhs != this)
   {
    create(rhs.get_n_nodes(),rhs.get_S_c(),rhs.get_K_h(),rhs.get_W_0(),
    rhs.get_gamma(),rhs.get_rho_s(),rhs.get_rho_r(),rhs.get_N(),
    rhs.get_N_0(),rhs.get_N_m(),rhs.get_beta(),rhs.get_K_g(),rhs.get_A(),
	rhs.get_A_bins(),rhs.get_b(),rhs.get_zeta(),rhs.get_eta(),rhs.get_h(),
	rhs.get_s_h(),rhs.get_s_b(),rhs.get_DeltaXh(),rhs.get_DeltaXb(),
	rhs.get_bin_edge_loc(),rhs.get_old_h(),rhs.get_old_zeta(),rhs.get_old_eta(),
	rhs.get_intermediate_zeta(),rhs.get_intermediate_h(),rhs.get_pre_surface_zeta(),
	rhs.get_pre_surface_h(),rhs.get_Mass_Flux(),rhs.get_fluff());
   }
  return *this;
 }



// this function loads a profile
void flowtube::load_profile(ifstream& initial_in)
{
	// set up the inital profile
	vector<double> zi(n_nodes);
	vector<double> ei(n_nodes);
	vector<double> hi(n_nodes);

	double placeholder;
	for (int i = 0; i<n_nodes; i++)
	{
		initial_in >> placeholder >> placeholder >> placeholder
				 >> zi[i] >> hi[i];
	  	ei[i] = zi[i] - hi[i];
	  	//cout << "z["<<i<<"]: " << zi[i] << " ei: " << ei[i] << endl;
	}
	zeta = zi;
	eta = ei;
	h = hi;
}


void flowtube::export_input_profile(ofstream& outfile)
{
	for (int i = 0; i<n_nodes; i++)
	{
		outfile << s_h[i]  << " " << zeta[i] << " " << h[i] << " "
				 << zeta[i] << " " << h[i] << endl;
	}
}


// this function raises the profile such that the downslope
// boundary is at ds_elev. Currently the ds_elev is set at 100 within the code. Not entirely sure why atm
void flowtube::raise_zeta_eta_ds_bound(double ds_elev)
{
	double ds_bound = zeta[n_nodes-1];
	double raise = ds_elev-ds_bound;

	for (int i = 0; i<n_nodes; i++)
	{
		zeta[i] = zeta[i]+raise;
		eta[i] = eta[i]+raise;
	}
}

// this function raises the profile such that the mean elevation
// is at mean_elev. Currently not called anywhere in the models
void flowtube::raise_zeta_eta_mean(double mean_elev)
{
	// get the mean elevation
	double mean = 0;
	for (int i = 0; i< n_nodes; i++)
	{
		mean+=zeta[i];
	}
	mean = mean/double(n_nodes);

	// now raise the profile
	double raise = mean_elev-mean;
	for (int i = 0; i<n_nodes; i++)
	{
		zeta[i] = zeta[i]+raise;
		eta[i] = eta[i]+raise;
	}

}

// this resets the flowtibe to have a flat profile
void flowtube::set_const_zeta_eta_h(double zeta_flat, double h_flat)
{
	for (int i = 0; i<n_nodes; i++)
	{
		zeta[i] = zeta_flat;
		h[i] = h_flat;
		eta[i] = zeta[i]-h[i];
	}
}

// this function sets up an interpolation routine
// you give the function the x locations of the values to be
// interpolated. The function then loops through the locations
// of the nodes in order to find the two adjacent nodes of each
// point. The function then updates two vectors,
// ds_interp_node_num and us_interp_node_num
// for each x-location of the value to be updated, the node numbers
// are stored
//
// IMPORTANT: the interpolated nodes are organized such that the
//				upslope node appears first.
void flowtube::initialize_interpolation_nodes(vector<double> interp_x_loc,
						vector<int>& ds_interp_node_num,
						vector<int>& us_interp_node_num,
						vector<double>& interpolated_distance_fraction)
{
	int n_interp_nodes = interp_x_loc.size();		// number of nodes in the vector to be
													// interpolated
	int h_node_number = 0;				// the node number of the downslope distance
	int interp_node_number = 0;			// node number of the interpolated vector
	vector<int> temp_ds_node_num;		// a temporary vector for storing distance_h nodes
										// adjacent to interpolated nodes (downslope).
	vector<int> temp_us_node_num;		// a temporary vector for storing distance_h nodes
										// adjacent to interpolated nodes (upslope).
	vector<double> temp_idf;			// a vector that holds the farction of distance
										// that an interpolated value is between nodes
	double idf;
	// loop through interpolating vector. The loop stops when it runs out of either
	// interpolating nodes or distance_h_nodes
	while (h_node_number < n_nodes-1 && interp_node_number <n_interp_nodes)
	{
		//cout << "interp_nn: " << interp_node_number << " and h_nn: " << h_node_number << endl;
		//cout << "interp_s_loc: " << interp_x_loc[interp_node_number] << endl;
		//cout << "us s: " << s_h[h_node_number] << " and ds_s: " << s_h[h_node_number+1] << endl;
		// if the interpolated node lies between distance_h nodes...
		if (s_h[h_node_number+1] > interp_x_loc[interp_node_number] &&
		    s_h[h_node_number] <= interp_x_loc[interp_node_number])
		{
			// store the index of the adjacent nodes
			temp_ds_node_num.push_back(h_node_number+1);
			temp_us_node_num.push_back(h_node_number);



			// get the fractional distance of teh interpolated value
			idf = (interp_x_loc[interp_node_number]-s_h[h_node_number])/
				  (s_h[h_node_number+1]-s_h[h_node_number]);

			//cout << "interp x: " << interp_x_loc[interp_node_number] << endl;
			//cout << "us s: " << s_h[h_node_number] << " and ds_s: " << s_h[h_node_number+1] << endl;
			//cout << "idf: " << idf << endl;

			temp_idf.push_back(idf);
			interp_node_number++;		// and increment the interpolating node number
		}
		else
		{
			h_node_number++;				// now increment the distance_h node number
											// the else statement is necessary for the case
											// when two interpolated node lie between the
											// same two downslope distance nodes.
		}
	}

	// once finished, update the interpolating referernce vectors
	ds_interp_node_num = temp_ds_node_num;
	us_interp_node_num = temp_us_node_num;
	interpolated_distance_fraction = temp_idf;

	// bug check:
	//cout << "n_interp_nodes: " << n_interp_nodes
	//     << " and n_indeces: " << ds_interp_node_num.size() << endl;
	//int n_interp = ds_interp_node_num.size();
	//for (int i = 0; i<ds_interp_node_num.size(); i++)
	//{
	//	cout << "interp_dist: " << interp_x_loc[i]
	//	     << " us_dist: " << s_h[ us_interp_node_num[i] ]
	//	     << " ds_dist: " << s_h[ ds_interp_node_num[i] ]
	//	     << endl;
	//}
}


// interpolating function:
// takes information generated using the initialize_interpolation_nodes
// function and returns the interpolaed values of h
// using linear interpolation
vector<double> flowtube::interpolated_modeled_h(vector<int> ds_interp_node_num,
						vector<int> us_interp_node_num,
						vector<double> interpolated_fractional_distance )
{
	vector<double> temp_interp_h;
	int n_interp_nodes = ds_interp_node_num.size();
	for(int i = 0; i<n_interp_nodes; i++)

	{
		//cout << "h_us: " << h[ us_interp_node_num[i]]
		//	 << " and h_ds: " << h[ ds_interp_node_num[i] ]
		//	 << " and frac_dist: " << interpolated_fractional_distance[i] << endl;

		temp_interp_h.push_back( (h[ ds_interp_node_num[i] ]
		                    -h[ us_interp_node_num[i] ])
		                    *interpolated_fractional_distance[i]
		                    +h[ us_interp_node_num[i] ] );
	}
	return temp_interp_h;
}

// interpolating function:
// takes information generated using the initialize_interpolation_nodes
// function and returns the interpolaed values of zeta
// using linear interpolation
vector<double> flowtube::interpolated_modeled_zeta(vector<int> ds_interp_node_num,
						vector<int> us_interp_node_num,
						vector<double> interpolated_fractional_distance )
{
	vector<double> temp_interp_zeta;
	int n_interp_nodes = ds_interp_node_num.size();
	for(int i = 0; i<n_interp_nodes; i++)
	{
		temp_interp_zeta.push_back( (zeta[ ds_interp_node_num[i] ]
		                    -zeta[ us_interp_node_num[i] ])
		                    *interpolated_fractional_distance[i]
		                    +zeta[ us_interp_node_num[i] ] );
	}
	return temp_interp_zeta;
}



// set the transport parameters
void flowtube::set_transport_params(double temp_S_c, double temp_K_h,
								double temp_W_0, double temp_gamma,
							  	double temp_rho_s, double temp_rho_r)
{
	S_c = temp_S_c;
	K_h = temp_K_h;
	W_0 = temp_W_0;
	gamma = temp_gamma;
	rho_s = temp_rho_s;
	rho_r = temp_rho_r;
}

void flowtube::set_transport_params(double temp_S_c, double temp_K_h,
								double temp_W_0, double temp_gamma,
							  	double temp_rho_s, double temp_rho_r,
							  	double temp_N, double temp_N_0,
							  	double temp_N_m,
							  	double temp_beta, double temp_K_g)
{
	S_c = temp_S_c;
	K_h = temp_K_h;
	W_0 = temp_W_0;
	gamma = temp_gamma;
	rho_s = temp_rho_s;
	rho_r = temp_rho_r;
	N =	temp_N;
	N_0 = temp_N_0;
	N_m = temp_N_m;
	beta = temp_beta;
	K_g = temp_K_g;
}
// Calculates a steady state flux for a erosion rate and upslope flux
double flowtube::calculate_steady_flux(double br_erosion_rate, double flux_us)
{
	double A_tot = 0;
	double flux_ds;
	for(int i = 0; i<n_nodes; i++)
	{
		A_tot+= A[i];
	}
	flux_ds = flux_us+br_erosion_rate*rho_r*A_tot;
	return flux_ds;
}


// this function calcualtes the evolution of the hillslope profile using
// a flux boundary condition at both ends
void flowtube::flux_timestep_flux_bc(double dt,
							double flux_us, double flux_ds,
							int flux_switch, int prod_switch,
							vector<double> surface_change_rate)
{
	// do a timestep
	double mean_thick;				// mean soil thickness
	double dx,dy,denom,sap_lowering_rate,mass_present;


	// reset the 'old' values
	old_h = h;
	old_eta = eta;
	old_zeta = zeta;

	// some data
	vector<double> prod(n_nodes);
	vector<double> sap_lowering(n_nodes);
	vector<double> dh_MassFlux(n_nodes);

	// the second set are stored at the half spaces
	// because there are two flux conditions, the half nodes are not stored at the edges,
	// so there are n_nodes-1 half nodes
	vector<double> slopes(n_nodes-1,0.0);
	vector<double> cos_theta(n_nodes-1,0.0);

	// the third set is just the fluxes, these occur at the box boundaries, including the
	// top and bottom boxes, so this vector has n_nodes+1 data elements
	vector<double> fluxes(n_nodes+1,0.0);	// fluxes in kg/yr

	// calculate fluxes
	// the two boundary fluxes are assinged
	fluxes[0] = flux_us;
	fluxes[n_nodes] = flux_ds;

	// calculate the slopes along the profile
	for (int i = 0; i<n_nodes-1; i++)
	{
		slopes[i] = (zeta[i+1]-zeta[i])/(DeltaXh[i]);
		dx = DeltaXh[i];
		dy = zeta[i+1]-zeta[i];
		cos_theta[i] = sqrt( dx*dx / (dx*dx + dy*dy) );
		///Print statements used to test if working
        // cout << "h: " << h[i] << " dx = " << dx << " dy: " << dy
		    // << " slopes: " << slopes[i] << " K: " << K_h << endl ;
	}

	// calcualte the MASS fluxes along the profile
	// these fluxes are in units of mass per time
	for (int i = 0; i<n_nodes-1; i++)
	{
		mean_thick = 0.5*(h[i+1]+h[i]);
		switch ( flux_switch )
		{
			case 1 :
			fluxes[i+1] = -rho_s*K_h*b[i+1]*mean_thick*slopes[i];
			break;
			case 2 :
			fluxes[i+1] =  -rho_s*K_h*b[i+1]*mean_thick*slopes[i]*cos_theta[i];
			break;
			case 3 :
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			fluxes[i+1] =  -rho_s*K_h*b[i+1]*mean_thick*slopes[i]*denom;
			break;
			case 4 :
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			fluxes[i+1] =  -rho_s*K_h*b[i+1]*mean_thick*slopes[i]
						  *denom*cos_theta[i];
			break;
			case 5 :
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*cos_theta[i]/beta));
			fluxes[i+1] = -K_g*N*rho_s*b[i+1]*slopes[i]*denom;
			break;
			case 6:
			N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*cos_theta[i]/beta));
			fluxes[i+1] = -K_g*N*rho_s*b[i+1]*slopes[i];
			break;
			case 7:
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*cos_theta[i]/beta));
			fluxes[i+1] = -K_g*N*rho_s*b[i+1]*slopes[i]*denom*mean_thick;
			break;
			default:
			fluxes[i+1] = -rho_s*K_h*b[i+1]*mean_thick*slopes[i];

		}
	}



	// update the hillslope characteristics
	for (int i = 0; i<n_nodes; i++)
	{
		switch (prod_switch)
		{
			case 1:
			sap_lowering_rate = W_0*exp(-h[i]/gamma);
			break;
			case 2:
			if (i == 0)
			sap_lowering_rate = (W_0/cos_theta[i])*exp(-h[i]*cos_theta[i]/gamma);
			else
			sap_lowering_rate = (W_0/cos_theta[i-1])*exp(-h[i]*cos_theta[i-1]/gamma);
			break;
			default :
			sap_lowering_rate = W_0*exp(-h[i]/gamma);
		}

		sap_lowering[i] = dt*sap_lowering_rate;
		eta[i] = old_eta[i] - dt*sap_lowering_rate;
		prod[i] = rho_r*A[i]*sap_lowering_rate;
		// cout << "sediment produced: " << prod[i] << endl;
		//cout << " h: " << h[i]
		//     << " p: " << sap_lowering_rate
		//     << " p[i]: " << sap_lowering[i]
		//     << " cos(theta) "<< cos_theta[0]
		//     << " dt: " << dt
		//     << " W_0: " << W_0
		//     << " gamma: " << gamma << endl;
        ///Calculates the mass present within the node box, if the mass present is less than the fluxes arriving downslope then modifies the downslope flux to account for this.
		mass_present = dt*fluxes[i] + dt*prod[i] + rho_s*A[i]*old_h[i];
		if (mass_present < fluxes[i+1]*dt)
			fluxes[i+1] = mass_present/dt;


		// calculate soil thickness and surface elevation from
		// creep like sediment transport only
		dh_MassFlux[i] = dt*(fluxes[i]-fluxes[i+1])/(rho_s*A[i]);
		intermediate_h[i] = dh_MassFlux[i] + old_h[i];
		intermediate_zeta[i] = old_eta[i]+intermediate_h[i];

		// calculate the fluffing factor (becuase material changes in
		// density when converted from bedrock to soil, there is a
		// vertical component of particle motion in the soil
		fluff[i] = sap_lowering_rate*dt*((rho_r/rho_s) - 1);

		// calculate surface elevation, and thickness after production and
		// sediment transport
		pre_surface_h[i] = dh_MassFlux[i] + dt*prod[i]/(rho_s*A[i]) + old_h[i];
		if (pre_surface_h[i] < 0)
			pre_surface_h[i] = 0;
		pre_surface_zeta[i] = pre_surface_h[i]+eta[i];

		// now calculate the soil thickness and surface elevation
		// after deposition or erosion on the surface
		h[i] = dt*surface_change_rate[i] + pre_surface_h[i];
		if (h[i] < 0)
			h[i] = 0;
		zeta[i] = h[i]+eta[i];
	}



	Mass_Flux = fluxes;			// this is a vector of mass fluxes evaluated at
//                                // the node boundaries
  // for (int i = 0; i<n_nodes+1; i++)
	// {
	// 	cout <<"fluxes at box edge "<< i << " is: "<< fluxes[i] << endl;
	//
	// }
//     cout <<"flux us:" << fluxes[0] << endl
//          <<"flux ds:" << fluxes[1] << endl
//          <<"flux ds:" << fluxes[2];
}

// this function updates the hillslope profile using a fixed elevation
// boudnary condition. The fixed elevation is ds_elevation
// ds_elevation is a dummy node, that lies beyond the final node of
// the flowtube. The mean soil thickness related to this dummy node
// is simply the soil thickness at the end of the flowtube
void flowtube::flux_timestep_elev_bc(double dt,
							double flux_us, double ds_elevation,
							int flux_switch, int prod_switch,
							vector<double> surface_change_rate)
{
	// do a timestep
	double mean_thick;				// mean soil thickness
	double dx,dy,denom,sap_lowering_rate,mass_present;
	double flux_ds;

	// reset the 'old' values
	old_h = h;
	old_eta = eta;
	old_zeta = zeta;

	// some data
	vector<double> prod(n_nodes);
	vector<double> sap_lowering(n_nodes);
	vector<double> dh_MassFlux(n_nodes);

	// the second set are stored at the half spaces
	// because there are two flux conditions, the half nodes are not stored at the edges,
	// so there are n_nodes-1 half nodes
	vector<double> slopes(n_nodes-1,0.0);
	vector<double> cos_theta(n_nodes-1,0.0);

	// the third set is just the fluxes, these occur at the box boundaries, including the
	// top and bottom boxes, so this vector has n_nodes+1 data elements
	vector<double> fluxes(n_nodes+1,0.0);	// fluxes in kg/yr

	// calculate fluxes
	// the two boundary fluxes are assinged
	fluxes[0] = flux_us;

	//cout << "LINE 619 FLUX_us: " << fluxes[0] << " and flux_ds: " << flux_ds << endl;
    /// So deltaXh runs from 0 to less than n_nodes-1 which means calling deltaXh[n_nodes-1] in code results in getting a value outside
    /// since deltaXh actually stores values from 0 to n_nodes-2. At least changing line 648 results in the code working seemingly fine
    /// Double check this with Simon since otherwise not sure.
	double ds_slope = (ds_elevation-zeta[n_nodes-1])/DeltaXh[n_nodes-2];
	double dx_ds = DeltaXh[n_nodes-2];
	double dy_ds = ds_elevation-zeta[n_nodes-1];
	double ds_cos_theta = sqrt( dx_ds*dx_ds / (dx_ds*dx_ds + dy_ds*dy_ds) );

    ///This is the product of far too much time identifying errors
    // cout << "zeta is: " << zeta[n_nodes-1] << " ds elevation is: " << ds_elevation << endl;
    // cout << "delta Xh is: " << DeltaXh[n_nodes+1324] <<endl;
    //  cout << "ds_slope: " << ds_slope
    //     << "dx_ds: " << dx_ds
    //      << "dy_ds: " << dy_ds
    //      << "ds_cos_theta: " << ds_cos_theta << endl;
	double mean_thick_ds = h[n_nodes-1];
	switch ( flux_switch )
	{
		case 1 :
		flux_ds = -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope;
		break;
		case 2 :
		flux_ds =  -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope*ds_cos_theta;
		break;
		case 3 :
		denom = 1/(1-ds_slope*ds_slope/(S_c*S_c));
		flux_ds =  -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope*denom;
		break;
		case 4 :
		denom = 1/(1-ds_slope*ds_slope/(S_c*S_c));
		flux_ds =  -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope
					  *denom*ds_cos_theta;
		break;
		case 5 :
		denom = 1/(1-ds_slope*ds_slope/(S_c*S_c));
		N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*ds_cos_theta/beta));
		flux_ds = -K_g*N*rho_s*b[n_nodes-1]*ds_slope*denom;
		break;
		case 6:
		N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*ds_cos_theta/beta));
		flux_ds = -K_g*N*rho_s*b[n_nodes-1]*ds_slope;
		break;
		case 7:
		denom = 1/(1-ds_slope*ds_slope/(S_c*S_c));
		N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*ds_cos_theta/beta));
		flux_ds = -K_g*N*rho_s*b[n_nodes-1]*ds_slope*denom*mean_thick_ds;
		break;
		default:
		flux_ds = -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope;
	}
	fluxes[n_nodes] = flux_ds;
    // cout <<"look at me I'm a potential error: " << flux_ds << endl;
	// calculate the slopes along the profile
	for (int i = 0; i<n_nodes-1; i++)
	{
		slopes[i] = (zeta[i+1]-zeta[i])/(DeltaXh[i]);
		dx = DeltaXh[i];
        //cout << "Dx is"<< dx << endl;
		dy = zeta[i+1]-zeta[i];
        ///Print statements used to see if working properly
        //cout << "zeta1: "<< zeta[i+1] << "zeta2:" << zeta[i] << "Dy is " << dy << endl;
		cos_theta[i] = sqrt( dx*dx / (dx*dx + dy*dy) );
        //cout << "Cos theta is" << cos_theta[i] << endl;
	}

	// calcualte the MASS fluxes along the profile
	// these fluxes are in units of mass per time
	for (int i = 0; i<n_nodes-1; i++)
	{
		mean_thick = 0.5*(h[i+1]+h[i]);
		switch ( flux_switch )
		{
			case 1 :
			fluxes[i+1] = -rho_s*K_h*b[i+1]*mean_thick*slopes[i];
			break;
			case 2 :
			fluxes[i+1] =  -rho_s*K_h*b[i+1]*mean_thick*slopes[i]*cos_theta[i];
			break;
			case 3 :
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			fluxes[i+1] =  -rho_s*K_h*b[i+1]*mean_thick*slopes[i]*denom;
			break;
			case 4 :
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			fluxes[i+1] =  -rho_s*K_h*b[i+1]*mean_thick*slopes[i]
						  *denom*cos_theta[i];
			break;
			case 5 :
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*cos_theta[i]/beta));
			fluxes[i+1] = -K_g*N*rho_s*b[i+1]*slopes[i]*denom;

            //cout << "N is: " << N << endl;
            //cout << "Cos theta is: " << cos_theta[i] << endl;
			break;
			case 6:
			N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*cos_theta[i]/beta));
			fluxes[i+1] = -K_g*N*rho_s*b[i+1]*slopes[i];
			break;
			case 7:
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*cos_theta[i]/beta));
			fluxes[i+1] = -K_g*N*rho_s*b[i+1]*slopes[i]*denom*mean_thick;
			break;
			default:
			fluxes[i+1] = -rho_s*K_h*b[i+1]*mean_thick*slopes[i];
		}
	}

	// update the hillslope characteristics
	for (int i = 0; i<n_nodes; i++)
	{
		switch (prod_switch)
		{
			case 1:
			sap_lowering_rate = W_0*exp(-h[i]/gamma);
			break;
			case 2:
			if (i == 0)
			sap_lowering_rate = (W_0/cos_theta[i])*exp(-h[i]*cos_theta[i]/gamma);
			else
			sap_lowering_rate = (W_0/cos_theta[i-1])*exp(-h[i]*cos_theta[i-1]/gamma);
			break;
			default :
			sap_lowering_rate = W_0*exp(-h[i]/gamma);
		}

		sap_lowering[i] = dt*sap_lowering_rate;
		eta[i] = old_eta[i] - dt*sap_lowering_rate;
		prod[i] = rho_r*A[i]*sap_lowering_rate;

		mass_present = dt*fluxes[i] + dt*prod[i] + rho_s*A[i]*old_h[i];
		if (mass_present < fluxes[i+1]*dt)
			fluxes[i+1] = mass_present/dt;


		// calculate soil thickness and surface eleavtion from
		// creep like sediment transport only
		dh_MassFlux[i] = dt*(fluxes[i]-fluxes[i+1])/(rho_s*A[i]);
		intermediate_h[i] = dh_MassFlux[i] + old_h[i];
		intermediate_zeta[i] = old_eta[i]+intermediate_h[i];

		// calculate the fluffing factor (becuase material changes in
		// density when converted from bedrock to soil, there is a
		// vertical component of particle motion in the soil
		fluff[i] = sap_lowering_rate*dt*((rho_r/rho_s) - 1);

		// calculate surface elevation, and thickness after production and
		// sediment transport
		pre_surface_h[i] = dh_MassFlux[i] + dt*prod[i]/(rho_s*A[i]) + old_h[i];
		if (pre_surface_h[i] < 0)
			pre_surface_h[i] = 0;
		pre_surface_zeta[i] = pre_surface_h[i]+eta[i];
		// now calculate the soil thickness and surface elevation
		// after deposition or erosion on the surface
		h[i] = dt*surface_change_rate[i] + pre_surface_h[i];
		if (h[i] < 0)
			h[i] = 0;
		zeta[i] = h[i]+eta[i];
        //cout << "eta: " << eta[i] << endl;
        //cout << "zeta: " << zeta[i] << endl;
	}

	Mass_Flux = fluxes;			// this is a vector of mass fluxes evaluated at
								// the node boundaries
	// for (int i = 0; i<n_nodes+1; i++)
	// 		{
	// 		cout <<"fluxes at box edge "<< i << " is: "<< fluxes[i] << endl;
	//
	// 		}
 }

///Third boundary condition
void flowtube::flux_timestep_varying_elev_bc(double dt,
							double flux_us, double ds_elevation,
							int flux_switch, int prod_switch,
							vector<double> surface_change_rate, double constant_surface_change_rate, double base_level_change)
{
	// do a timestep
	double mean_thick;				// mean soil thickness
	double dx,dy,denom,sap_lowering_rate,mass_present;
	double flux_ds;

	// reset the 'old' values
	old_h = h;
	old_eta = eta;
	old_zeta = zeta;

	// some data
	vector<double> prod(n_nodes);
	vector<double> sap_lowering(n_nodes);
	vector<double> dh_MassFlux(n_nodes);

	// the second set are stored at the half spaces
	// because there are two flux conditions, the half nodes are not stored at the edges,
	// so there are n_nodes-1 half nodes
	vector<double> slopes(n_nodes-1,0.0);
	vector<double> cos_theta(n_nodes-1,0.0);

	// the third set is just the fluxes, these occur at the box boundaries, including the
	// top and bottom boxes, so this vector has n_nodes+1 data elements
	vector<double> fluxes(n_nodes+1,0.0);	// fluxes in kg/yr

	// calculate fluxes
	// the two boundary fluxes are assinged
	fluxes[0] = flux_us;

	//cout << "LINE 619 FLUX_us: " << fluxes[0] << " and flux_ds: " << flux_ds << endl;
    /// So deltaXh runs from 0 to less than n_nodes-1 which means calling deltaXh[n_nodes-1] in code results in getting a value outside
    /// since deltaXh actually stores values from 0 to n_nodes-2. At least changing line 648 results in the code working seemingly fine
    /// Double check this with Simon since otherwise not sure.
	double ds_slope = (ds_elevation-zeta[n_nodes-1])/DeltaXh[n_nodes-2];
	double dx_ds = DeltaXh[n_nodes-2];
	double dy_ds = ds_elevation-zeta[n_nodes-1];
	double ds_cos_theta = sqrt( dx_ds*dx_ds / (dx_ds*dx_ds + dy_ds*dy_ds) );

    // /This is the product of far too much time identifying errors
    // cout << "zeta is: " << zeta[n_nodes-1] << " ds elevation is: " << ds_elevation << endl;
    // cout << "delta Xh is: " << DeltaXh[n_nodes-2] <<endl;
    //  cout << "ds_slope: " << ds_slope
    //     << "dx_ds: " << dx_ds
    //      << "dy_ds: " << dy_ds
    //      << "ds_cos_theta: " << ds_cos_theta << endl;
	double mean_thick_ds = h[n_nodes-1];
	switch ( flux_switch )
	{
		case 1 :
		flux_ds = -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope;
		break;
		case 2 :
		flux_ds =  -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope*ds_cos_theta;
		break;
		case 3 :
		denom = 1/(1-ds_slope*ds_slope/(S_c*S_c));
		flux_ds =  -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope*denom;
		break;
		case 4 :
		denom = 1/(1-ds_slope*ds_slope/(S_c*S_c));
		flux_ds =  -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope
					  *denom*ds_cos_theta;
		break;
		case 5 :
		denom = 1/(1-ds_slope*ds_slope/(S_c*S_c));
		N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*ds_cos_theta/beta));
		flux_ds = -K_g*N*rho_s*b[n_nodes-1]*ds_slope*denom;
		break;
		case 6:
		N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*ds_cos_theta/beta));
		flux_ds = -K_g*N*rho_s*b[n_nodes-1]*ds_slope;
		break;
		case 7:
		denom = 1/(1-ds_slope*ds_slope/(S_c*S_c));
		N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*ds_cos_theta/beta));
		flux_ds = -K_g*N*rho_s*b[n_nodes-1]*ds_slope*denom*mean_thick_ds;
		break;
		default:
		flux_ds = -rho_s*K_h*b[n_nodes-1]*mean_thick_ds*ds_slope;
	}
	fluxes[n_nodes] = flux_ds;
	// cout << "looks at me I'm a potenital error: " << flux_ds << endl;
	// double ds_elev_temp = update_ds_elev(dt, base_level_change, prod_switch,
	// 			 constant_surface_change_rate, ds_elevation);
	// 				ds_elevation = ds_elev_temp;
	// 					cout << "New ds_elev is: "<< ds_elev_temp << endl;

	// calculate the slopes along the profile
	for (int i = 0; i<n_nodes-1; i++)
	{
		slopes[i] = (zeta[i+1]-zeta[i])/(DeltaXh[i]);
		dx = DeltaXh[i];
        //cout << "Dx is"<< dx << endl;
		dy = zeta[i+1]-zeta[i];
        ///Print statements used to see if working properly
        //cout << "zeta1: "<< zeta[i+1] << "zeta2:" << zeta[i] << "Dy is " << dy << endl;
		cos_theta[i] = sqrt( dx*dx / (dx*dx + dy*dy) );
        //cout << "Cos theta is" << cos_theta[i] << endl;
	}

	// calcualte the MASS fluxes along the profile
	// these fluxes are in units of mass per time
	for (int i = 0; i<n_nodes-1; i++)
	{
		mean_thick = 0.5*(h[i+1]+h[i]);
		switch ( flux_switch )
		{
			case 1 :
			fluxes[i+1] = -rho_s*K_h*b[i+1]*mean_thick*slopes[i];
			break;
			case 2 :
			fluxes[i+1] =  -rho_s*K_h*b[i+1]*mean_thick*slopes[i]*cos_theta[i];
			break;
			case 3 :
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			fluxes[i+1] =  -rho_s*K_h*b[i+1]*mean_thick*slopes[i]*denom;
			break;
			case 4 :
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			fluxes[i+1] =  -rho_s*K_h*b[i+1]*mean_thick*slopes[i]
						  *denom*cos_theta[i];
			break;
			case 5 :
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*cos_theta[i]/beta));
			fluxes[i+1] = -K_g*N*rho_s*b[i+1]*slopes[i]*denom;

            //cout << "N is: " << N << endl;
            //cout << "Cos theta is: " << cos_theta[i] << endl;
			break;
			case 6:
			N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*cos_theta[i]/beta));
			fluxes[i+1] = -K_g*N*rho_s*b[i+1]*slopes[i];
			break;
			case 7:
			denom = 1/(1-slopes[i]*slopes[i]/(S_c*S_c));
			N = N_m/(1 + (N_m/N_0-1)*exp(-mean_thick*cos_theta[i]/beta));
			fluxes[i+1] = -K_g*N*rho_s*b[i+1]*slopes[i]*denom*mean_thick;
			break;
			default:
			fluxes[i+1] = -rho_s*K_h*b[i+1]*mean_thick*slopes[i];
		}
	}

	// update the hillslope characteristics
	for (int i = 0; i<n_nodes; i++)
	{
		switch (prod_switch)
		{
			case 1:
			sap_lowering_rate = W_0*exp(-h[i]/gamma);
			break;
			case 2:
			if (i == 0)
			sap_lowering_rate = (W_0/cos_theta[i])*exp(-h[i]*cos_theta[i]/gamma);
			else
			sap_lowering_rate = (W_0/cos_theta[i-1])*exp(-h[i]*cos_theta[i-1]/gamma);
			break;
			default :
			sap_lowering_rate = W_0*exp(-h[i]/gamma);
		}

		sap_lowering[i] = dt*sap_lowering_rate;
		eta[i] = old_eta[i] - dt*sap_lowering_rate;
		prod[i] = rho_r*A[i]*sap_lowering_rate;

		mass_present = dt*fluxes[i] + dt*prod[i] + rho_s*A[i]*old_h[i];
		if (mass_present < fluxes[i+1]*dt)
			fluxes[i+1] = mass_present/dt;


		// calculate soil thickness and surface eleavtion from
		// creep like sediment transport only
		dh_MassFlux[i] = dt*(fluxes[i]-fluxes[i+1])/(rho_s*A[i]);
		intermediate_h[i] = dh_MassFlux[i] + old_h[i];
		intermediate_zeta[i] = old_eta[i]+intermediate_h[i];

		// calculate the fluffing factor (becuase material changes in
		// density when converted from bedrock to soil, there is a
		// vertical component of particle motion in the soil
		fluff[i] = sap_lowering_rate*dt*((rho_r/rho_s) - 1);

		// calculate surface elevation, and thickness after production and
		// sediment transport
		pre_surface_h[i] = dh_MassFlux[i] + dt*prod[i]/(rho_s*A[i]) + old_h[i];
		if (pre_surface_h[i] < 0)
			pre_surface_h[i] = 0;
		pre_surface_zeta[i] = pre_surface_h[i]+eta[i];
		// now calculate the soil thickness and surface elevation
		// after deposition or erosion on the surface
		h[i] = dt*surface_change_rate[i] + pre_surface_h[i];
		if (h[i] < 0)
			h[i] = 0;
		zeta[i] = h[i]+eta[i];
        //cout << "eta: " << eta[i] << endl;
        //cout << "zeta: " << zeta[i] << endl;
	}
	// for (int i = 0; i<n_nodes+1; i++)
 	// 		{
 	// 		cout <<"zeta change  "<< i << " is: "<< old_zeta[i]-zeta[i] << endl;
	//
 	// 		}


  old_h = h;
	Mass_Flux = fluxes;			// this is a vector of mass fluxes evaluated at
								// the node boundaries
	// for (int i = 0; i<n_nodes+1; i++)
	// 		{
	// 		cout <<"fluxes at box edge "<< i << " is: "<< fluxes[i] << endl;
	//
	// 		}

}
///Function to update the lowering rate of ds-elev based on the soil thickness of the final node.
double flowtube::update_ds_elev(double dt,double base_level_change,int prod_switch,
			double constant_surface_change_rate, double ds_elev,vector<double> bc_h)
{
	///The ds lowering rate variable
	// double ds_lowering_rate;
	double ds_elev_new;
	// double ds_slope = (ds_elev-zeta[n_nodes-1])/DeltaXh[n_nodes-2];
	// double dx_ds = DeltaXh[n_nodes-2];
	// double dy_ds = ds_elev-zeta[n_nodes-1];
	// double ds_cos_theta = sqrt( dx_ds*dx_ds / (dx_ds*dx_ds + dy_ds*dy_ds) );
	// cout << "The H used will be: " << bc_h[n_nodes-1] << endl;
	// cout << "zeta is: " << zeta[n_nodes-1] << " ds elevation is: " << ds_elev << endl;
	// cout << "old h is:" << old_h[n_nodes-1] << endl;
	// cout << "delta Xh is: " << DeltaXh[n_nodes-2] <<endl;
	//  cout << "ds_slope: " << ds_slope
	// 		<< "dx_ds: " << dx_ds
	// 		 << "dy_ds: " << dy_ds
	// 		 << "ds_cos_theta: " << ds_cos_theta << endl;
		// switch (prod_switch)
		// {
		// 	case 1:
		// 	ds_lowering_rate = W_0*exp(-bc_h[n_nodes-1] /gamma);
		// 	break;
		// 	case 2:
		//
		// 	ds_lowering_rate = (W_0/ds_cos_theta)*exp(-bc_h[n_nodes-1] *ds_cos_theta/gamma);
		// 	break;
		// 	default :
		// 	ds_lowering_rate = W_0*exp(-bc_h[n_nodes-1] /gamma);
		// }

	 // ds_elev -= dt*(base_level_change-constant_surface_change_rate+ds_lowering_rate);
	 ds_elev += (zeta[n_nodes-1]-bc_h[n_nodes-1])-(dt*base_level_change);

	 ds_elev_new = ds_elev;
		return ds_elev_new;
}
//=-=-=-=-=-=-=-=-=-=-=
// printing functions
//=-=-=-=-=-=-=-=-=-=-=
void flowtube::print_h_s(ofstream& outf)
{
	outf << "-99";
	for (int i = 0; i< n_nodes; i++)
	{
		outf << " " << s_h[i];
	}
	outf << endl;
}

void flowtube::print_zeta(double t_ime, ofstream& zeta_out)
{
	zeta_out << t_ime;
	for (int i = 0; i< n_nodes; i++)
	{
		zeta_out << " " << zeta[i];
	}
	zeta_out << endl;
}

void flowtube::print_relative_zeta(double t_ime, ofstream& zeta_out)
{
	double mean_zeta =0;
	for (int i = 0; i< n_nodes; i++)
	{
		mean_zeta+=zeta[i];
	}
	mean_zeta = mean_zeta/double(n_nodes);

	zeta_out << t_ime;
	for (int i = 0; i< n_nodes; i++)
	{
		zeta_out << " " << zeta[i]-mean_zeta;
	}
	zeta_out << endl;
}


void flowtube::print_eta(double t_ime, ofstream& eta_out)
{
	eta_out << t_ime;
	for (int i = 0; i< n_nodes; i++)
	{
		eta_out << " " << eta[i];
	}
	eta_out << endl;
}

void flowtube::print_h(double t_ime, ofstream& h_out)
{
	h_out << t_ime;
	for (int i = 0; i< n_nodes; i++)
	{
		h_out << " " << h[i];
	}
	h_out << endl;
}

void flowtube::print_ft_vecs_to_screen()
{
	cout << endl;
	for (int i = 0; i< n_nodes; i++)
	{
			cout << "i: " << i << " eta: " << eta[i] << " zeta: "
				 << zeta[i] << " h: " << h[i] << endl;
	}
	cout << endl;
}

// prints the flowtube properties to file
// note that the flowtube widths are offset from the
// other data members, the b location are at intermediate
// half nodes
void flowtube::print_ft_properties(ofstream& outfile)
{
	outfile << "s b A zeta eta h" << endl;
	for (int i = 0; i<n_nodes; i++)
	{
		outfile << s_h[i] << " " << b[i] << " " << A[i] << " "
				<< zeta[i] << " " << eta[i] << " " << h[i] << endl;
	}
	outfile << "-99 " << b[n_nodes] << " -99 -99 -99 -99" << endl;
}
