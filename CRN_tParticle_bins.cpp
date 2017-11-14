#ifndef CRN_tParticle_bins_CPP
#define CRN_tParticle_bins_CPP

#include <iostream>
#include <vector>
#include <list>
#include <time.h>
#include <math.h>
#include "tParticle.hpp"
#include "chronos_particle_info.hpp"

#include "flowtube.hpp"
#include "CRN_tParticle_bins.hpp"
#include "VolumeParticleInfo.hpp"

using namespace std;


void CRN_tParticle_bins::create()
{
	cout << "you need to initialize with a vector!" << endl;
}

void CRN_tParticle_bins::create(vector<double> bel,vector<double> s_h,
						vector<double> s_b, vector<double> b,
						vector<double> A_bin_temp)
{
	A_bins = A_bin_temp;
	bin_edge_loc = bel;
	int sz_bel = bel.size();
	n_bins = sz_bel-1;
	vector< list<CRN_tParticle> > temp_pbins(n_bins);
	particle_bins = temp_pbins;
	//cout << "CRN_tParticle_bins.h LINE 20 n_bins = " << particle_bins.size() << endl;

	vector<double> bv_temp(n_bins);
	vector<int> bv_nodes_temp(n_bins);

	h_node_us = bv_nodes_temp;
	h_node_ds = bv_nodes_temp;
	b_node_us = bv_nodes_temp;
	b_node_ds = bv_nodes_temp;

	dx_h = bv_temp;
	dx_b = bv_temp;
	s_us_h = bv_temp;
	s_us_b = bv_temp;
	b_us = bv_temp;
	Slope_b = bv_temp;
	bin_width = bv_temp;

	int n_h = s_h.size();

	//for(int i = 0; i< s_b.size(); i++)
	// cout << b[i] << endl;
	//cout << endl;

	for (int bn = 0; bn< n_bins; bn++)
	{
		bin_width[bn] = bin_edge_loc[bn+1]-bin_edge_loc[bn];
		if (bn == 0)
		{
			h_node_us[bn] = 0;
			h_node_ds[bn] = 0;
			dx_h[bn] = s_h[1]-s_h[0];
			s_us_h[bn] = bin_edge_loc[0];
		}
		else if (bn == n_bins-1)
		{
			h_node_us[bn] = n_h-1;
			h_node_ds[bn] = n_h-1;
			dx_h[bn] = s_h[n_h-1]-s_h[n_h-2];
			s_us_h[bn] = bin_edge_loc[n_bins-1];
		}
		else
		{
			h_node_us[bn] = (bn-1)/2;
			h_node_ds[bn] = (bn-1)/2+1;
			dx_h[bn] = s_h[h_node_ds[bn]]-s_h[h_node_us[bn]];
			s_us_h[bn] = s_h[h_node_us[bn]];
		}

		b_node_us[bn] = bn/2;
		b_node_ds[bn] = bn/2+1;
		dx_b[bn] = s_b[b_node_ds[bn]]-s_b[b_node_us[bn]];
		b_us[bn] = b[b_node_us[bn]];

		//cout << "LINE 81 CRN_tP_bins.cpp bin number: " << bn << " b_ds: " << b[b_node_ds[bn]]
		//     << " b_us: " << b[b_node_us[bn]] << endl;

		Slope_b[bn] = (b[b_node_ds[bn]]-b[b_node_us[bn]])/dx_b[bn];
		s_us_b[bn] = s_b[b_node_us[bn]];
	}
}

// this creates the bins using a flowtube object
void CRN_tParticle_bins::create(flowtube ft)
{
	//vector<double> bin_edge_loc = ft.get_bin_edge_loc();
	bin_edge_loc = ft.get_bin_edge_loc();
	vector<double> b = ft.get_b();
	vector<double> s_h = ft.get_s_h();
	vector<double> s_b = ft.get_s_b();
	A_bins = ft.get_A_bins();

	int sz_bel = bin_edge_loc.size();
	n_bins = sz_bel-1;
	vector< list<CRN_tParticle> > temp_pbins(n_bins);
	particle_bins = temp_pbins;
	//cout << "CRN_tParticle_bins.h LINE 20 n_bins = " << particle_bins.size() << endl;

	vector<double> bv_temp(n_bins);
	vector<int> bv_nodes_temp(n_bins);

	h_node_us = bv_nodes_temp;
	h_node_ds = bv_nodes_temp;
	b_node_us = bv_nodes_temp;
	b_node_ds = bv_nodes_temp;

	dx_h = bv_temp;
	dx_b = bv_temp;
	s_us_h = bv_temp;
	s_us_b = bv_temp;
	b_us = bv_temp;
	Slope_b = bv_temp;
	bin_width = bv_temp;

	int n_h = s_h.size();

	//for(int i = 0; i< s_b.size(); i++)
	// cout << b[i] << endl;
	//cout << endl;

	for (int bn = 0; bn< n_bins; bn++)
	{
		bin_width[bn] = bin_edge_loc[bn+1]-bin_edge_loc[bn];
		if (bn == 0)
		{
			h_node_us[bn] = 0;
			h_node_ds[bn] = 0;
			dx_h[bn] = s_h[1]-s_h[0];
			s_us_h[bn] = bin_edge_loc[0];
		}
		else if (bn == n_bins-1)
		{
			h_node_us[bn] = n_h-1;
			h_node_ds[bn] = n_h-1;
			dx_h[bn] = s_h[n_h-1]-s_h[n_h-2];
			s_us_h[bn] = bin_edge_loc[n_bins-1];
		}
		else
		{
			h_node_us[bn] = (bn-1)/2;
			h_node_ds[bn] = (bn-1)/2+1;
			dx_h[bn] = s_h[h_node_ds[bn]]-s_h[h_node_us[bn]];
			s_us_h[bn] = s_h[h_node_us[bn]];
		}

		b_node_us[bn] = bn/2;
		b_node_ds[bn] = bn/2+1;
		dx_b[bn] = s_b[b_node_ds[bn]]-s_b[b_node_us[bn]];
		b_us[bn] = b[b_node_us[bn]];

		//cout << "bin number: " << bn << " b_ds: " << b[b_node_ds[bn]]
		//     << " b_us: " << b[b_node_us[bn]] << endl;

		Slope_b[bn] = (b[b_node_ds[bn]]-b[b_node_us[bn]])/dx_b[bn];
		s_us_b[bn] = s_b[b_node_us[bn]];
	}

	//cout << "YOYOYOYOYOYOYO bel sz = " << bin_edge_loc.size() << endl;
}


// NOTE: delta lowered is the same size as eta, each h node has one lowering thickness
// (thus affecting two bins each)
// this inserts particles into the simulation. One keeps track of soil production,
// and after some time T the total lowering of the soil bedrock boundary is calculated,
// this lowering is Delta_lowered = integral from T to T+Delta T of p dt where p is
// soil production in m/yr. Thus Delta_lowered is in metres, and a potentially different
// value exists for each h node.
// The particles are introduced in a box whose upper boundary is the start depth
// below the _current_ eta elevation, and whose bottom boundary is this top zone minus
// the distance Delta_lowered.
// part_conc is the number of particles per kg of rock
// rho_s and rho_r are used to calculate both the number of particles introduced and
// the effective depth of the particles.
//
// the pID is a particle identification number, each particle is tagged with a different integer
// so in plotting algorithms one particle at a time can be followed
// the function returns and integer, that is the value of the
// pID for the first particle of the next inserting sequence.
int CRN_tParticle_bins::insert_particles(flowtube ft,
									vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double part_conc, vector<int> starting_pID, vector<double> starting_p_mfrac)
{

	vector<double> new_bottom_depth = old_bottom_depth;
	vector<double> zeta = ft.get_zeta();		// surface elevation (m)
	vector<double> eta = ft.get_eta();			// soil-bedrock boundary of the bin (m)


	double rho_r = ft.get_rho_r();
	double rho_s = ft.get_rho_s();

	//cout << endl << "LINE 176 CRN_tP_bins" << endl;
	//for (int i = 0; i< zeta.size(); i++)
	// cout << "zeta["<<i<<"]: " << zeta[i] << " and eta[" <<i<< "]: " << eta[i] << endl;
	//cout << endl;

	long seed = time(NULL);             // seed for random number generator
	int n_parts_inserted;				// number of particles inserted into a bin
	double insert_zone_top;				// top of insertion zone (m)
	double insert_zone_bottom;			// bottom of insertion zone (m)
	double ran_sl;						// downslope distance of the randomly inserted particle (m)
	double ran_zl;						// the elevation of the randomly inserted particle (m)
	//double theta;						//
	int eta_node;						// teh node index for eta, h, and zeta

	// interpolated hillslope properties
	double interpolated_zeta;			// surface elevation (m)
	double interpolated_eta;			// soil-bedrock boundary (m)
	double interpolated_h;				// soil thickness (m)
	double d_rock;						// depth of rock overlying particle (m)
	double d;							// depth of particle (m)
	double eff_d;						// effective depth of particle (g/cm^2)

	// during a timestep some particles are inserted into a zone in the bedrock
	// the number of particles is defined by the volume of the insertion zone times
	// the particle concnetration. The volume times the concentration will typicaly
	// be an integer, so we need a way to maintin some kind of rasonably constant
	// and unbiased particle concentration (i.e., avoiding rounding bias)
	// the approach is to get thh fractional part of the number of particles
	// and then select a random number from a uniform distribution on 0 to 1
	// , if the ranom number is less than the fraction a particle is generated,
	// if not no particle. The below variables are for making this calculation
	double n_parts_doub;
	double n_parts_int;
	double n_parts_frac;

	// this function loops through the bins, inserting particles based on how
	// much the soil/saprolite boundary has lowered and the particle concnetration,
	// which is given in particles/kg
	//
	// in general we consider the surfaces zeta and eta to be sloping, but for
	// calculation of the insertion zone the eta location is considered horizontal
	// this is because if the insertion zone was considered sloping then transient changes
	// in the bedrock might lead to ovelapping insertion zones. A horizontal insertion
	// zone ensures that particles are always inserted into 'fresh' rock
	for (int bn = 0; bn<n_bins; bn++)
	{
		// get the number of particles inserted
		n_parts_doub = A_bins[bn]*Delta_lowered[bn/2]*rho_r*part_conc;
		n_parts_frac = modf(n_parts_doub,&n_parts_int);
		n_parts_inserted = int(n_parts_int);
		//cout << "LINE 211 CRN_tP_bins np_doub: " << n_parts_doub
		//	 << " np_frac: " << n_parts_frac
		//     << " np_ins: " << n_parts_inserted << endl
		//     << "	A_bin: "<< A_bins[bn] << " and Dl: " << Delta_lowered[bn/2]
		//     << " rho_r: " << rho_r << " and part_conc: " << part_conc << endl;
		if( ran3(&seed) <= n_parts_frac )
		 n_parts_inserted++;
		//cout << " np_ins: " << n_parts_inserted << endl;

		eta_node = bn/2;
		insert_zone_top = old_bottom_depth[eta_node];
		insert_zone_bottom = insert_zone_top - Delta_lowered[eta_node];
		new_bottom_depth[eta_node] = insert_zone_bottom;

		//cout << "LINE 217 CRN_tP_bins, bin_number: "<< bn
		//     << " eta: " << eta[eta_node]
		//	 << " zeta: " << zeta[eta_node]
		//	 << " Delta_eta: " << Delta_lowered[eta_node] << endl
		//	 << " Ins zone top: " << insert_zone_top
		//	 << " Ins zone bot: " << insert_zone_bottom << endl;

		for (int part = 0; part<n_parts_inserted; part++)
		{
			// randomly generate the s and z locations
			ran_sl = (bin_edge_loc[bn+1] - bin_edge_loc[bn])*ran3(&seed) +
					  bin_edge_loc[bn];

			ran_zl = (insert_zone_top-insert_zone_bottom)*ran3(&seed) +
					 insert_zone_bottom;

			//cout << "Line 240 CRN_tP_bins, ran_sl: " << ran_sl
			//     << " ran_zl " << ran_zl << endl;

			// now create and place the particles
			// for a given z and s location, you must calculate the depth and effective
			// d location
			// the effectve d location is in units g/cm^2, everything else is
			// in SI

			// first interpolate zeta and eta:
			interpolated_zeta = ((zeta[ h_node_ds[bn] ] - zeta[ h_node_us[bn] ])/
								dx_h[bn])*(ran_sl- s_us_h[bn]) + zeta[ h_node_us[bn] ];
			interpolated_eta = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*(ran_sl- s_us_h[bn]) + eta[ h_node_us[bn] ];

			//cout << "	interpolated zeta: " << interpolated_zeta
			//     << " and interpolated_eta: " << interpolated_eta << endl;

			//cout << endl << endl;
			//cout << h_node_ds[bn] << " " << h_node_us[bn] << " "
			//     << zeta[ h_node_ds[bn] ] << " " << zeta[ h_node_us[bn] ] << endl;
			//cout << endl << endl;

			interpolated_h = interpolated_zeta-interpolated_eta;
			d = interpolated_zeta-ran_zl;

			// if the depth is greater than the soil thickness,
			// the effective depth must take into account
			// the density change
			if (ran_zl<interpolated_eta)
			{
				d_rock = interpolated_eta-ran_zl;
				eff_d = 0.1*(rho_r*d_rock+rho_s*interpolated_h);
						// the 0.1 factor is to convert from
						// kg/m^2 to g/cm^2
			}
			// if the depth is less than the soil thickness,
			// the effective depth only accounts for the density
			// of the soil
			else
			{
				eff_d = 0.1*(rho_s*d);
						// the 0.1 factor is to convert from
						// kg/m^2 to g/cm^2
				//cout << "LINE 325 the inserted particle is in the soil! ";

			}

			// now create the particle:
			// determine the starting type
			double type_prob = ran3(&seed);
			int n_ptypes = starting_pID.size();
			int startType = starting_pID[n_ptypes-1];
			double cum_prob = 0.0;
			double last_prob = 0.0;
			for (int ptype_index = 0; ptype_index<n_ptypes; ptype_index++)
			{
				last_prob = cum_prob;
				cum_prob = cum_prob+starting_p_mfrac[ptype_index];

				if (type_prob >= last_prob && type_prob < cum_prob)
				{
					startType = starting_pID[ptype_index];
				}
			}

			CRN_tParticle ins_part(startType, ran_sl,d,eff_d, ran_zl);
			if (ran_zl>=interpolated_eta)
			{
				ins_part.SoilAgeExpose();
			}
			particle_bins[bn].push_back(ins_part);
			startType++; 				// increment the particle type

		}

	}
	old_bottom_depth = new_bottom_depth;

	//for (int i = 0; i< n_bins; i++)
	//{
	//	cout << "LINE 347 CRN_tParticle_bins.cpp; n_parts of bin " << i << " = " << particle_bins[i].size() << endl;
	//}
	return 1;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// NOTE: delta lowered is the same size as eta, each h node has one lowering thickness
// (thus affecting two bins each)
// this inserts particles into the simulation. One keeps track of soil production,
// and after some time T the total lowering of the soil bedrock boundary is calculated,
// this lowering is Delta_lowered = integral from T to T+Delta T of p dt where p is
// soil production in m/yr. Thus Delta_lowered is in metres, and a potentially different
// value exists for each h node.
// The particles are introduced in a box whose upper boundary is the start depth
// below the _current_ eta elevation, and whose bottom boundary is this top zone minus
// the distance Delta_lowered.
// part_conc is the number of particles per kg of rock
// rho_s and rho_r are used to calculate both the number of particles introduced and
// the effective depth of the particles.
//
// the pID is a particle identification number, each particle is tagged with a different integer
// so in plotting algorithms one particle at a time can be followed
// the function returns and integer, that is the value of the
// pID for the first particle of the next inserting sequence.
//
// this overloaded version of the function sets the scaling and parameter values for all of the particles
// CRN_muon_param_switch:
//		0: Granger
//		1: Schaller
//
// TO DO: need to implement an insertiona algorithm that can cope with a 'deep' 
// insertion wherein both the
// soil column and parent material are populated.
// probably the best approach is this: check to see if the insertion zone
// goes across the boundary. If it does,
// then split the insertion zones into two parts
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int CRN_tParticle_bins::insert_particles_volumetric(flowtube ft,
									vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double C_10Be, double C_26Al, double C_36Cl, double C_14C,
                  double C_21Ne, double C_3He,
									VolumeParticleInfo vpi)
{

	vector<double> new_bottom_depth = old_bottom_depth;
	vector<double> zeta = ft.get_zeta();		// surface elevation (m)
	vector<double> eta = ft.get_eta();			// soil-bedrock boundary of the bin (m)

	double rho_r = ft.get_rho_r();
	double rho_s = ft.get_rho_s();

	long seed = time(NULL);             // seed for random number generator
	int n_parts_inserted;				// number of particles inserted into a bin
	double insert_zone_top;				// top of insertion zone (m)
	double insert_zone_bottom;			// bottom of insertion zone (m)
	double ran_sl;						// downslope distance of the randomly inserted particle (m)
	double ran_zl;						// the elevation of the randomly inserted particle (m)
	//double theta;						//
	int eta_node;						// the node index for eta, h, and zeta

	// interpolated hillslope properties
	double interpolated_zeta;			// surface elevation (m)
	double interpolated_eta;			// soil-bedrock boundary (m)
	double interpolated_h;				// soil thickness (m)
	double d_rock;						// depth of rock overlying particle (m)
	double d;							// depth of particle (m)
	double eff_d;						// effective depth of particle (g/cm^2)

	// during a timestep some particles are inserted into a zone in the bedrock
	// the number of particles is determined through a target volume
	double n_parts_double;
	double n_parts_frac;
	double particle_mass;
	double soil_particle_mass;
	double soil_particle_surface_area;
	double density_ratio = rho_s/rho_r;
	double insert_particle_mass;
	double insert_particle_surface_area;
	double particle_surface_area;
	double cell_volume;
	double cell_mass;

	// get the number of types and sizes of the particles
	int n_types = vpi.get_n_types();
	int n_sizes = vpi.get_n_sizes();

	// the mass of each type/size
	double mass_of_type_size;

	// this function loops through the bins, inserting particles based on how
	// much the soil/saprolite boundary has lowered and the particle concentration,
	// which is given in particles/kg
	//
	// in general we consider the surfaces zeta and eta to be sloping, but for
	// calculation of the insertion zone the eta location is considered horizontal
	// this is because if the insertion zone was considered sloping then transient changes
	// in the bedrock might lead to overlapping insertion zones. A horizontal insertion
	// zone ensures that particles are always inserted into 'fresh' rock
	for (int bn = 0; bn<n_bins; bn++)
	{
		// get the number of particles inserted
		cell_volume = A_bins[bn]*Delta_lowered[bn/2];
		cell_mass = cell_volume*rho_r;					// mass calucalted based on rock density

		eta_node = bn/2;
		insert_zone_top = old_bottom_depth[eta_node];
		insert_zone_bottom = insert_zone_top - Delta_lowered[eta_node];
		new_bottom_depth[eta_node] = insert_zone_bottom;

		cout << "particle bins line 479 cell volume = " << cell_volume << " and mass: " << cell_mass << endl;
		cout << "zone top: " << insert_zone_top << " and bottom: " << insert_zone_bottom << endl;

    //cout << "Line 478, bel: " << bin_edge_loc[bn] << " and ds bel: "  << bin_edge_loc[bn+1] << endl;

		// loop through each particle type and each size class
		for (int type = 0; type<n_types; type++)
		{
			for (int sizet = 0; sizet< n_sizes; sizet++)
			{
				// get the mass of these particle types within the insertion zone
				mass_of_type_size = cell_mass*vpi.return_mass_fraction(type,sizet);

				// only generate particles if there is mass
				if (mass_of_type_size > 0)
				{
					// now get the number of particles of this type/size
					n_parts_double = mass_of_type_size/vpi.get_ParticleTargetMass();

					// now get the integer version of this:
					n_parts_inserted = int(n_parts_double);
					
					//cout << "type: " << type << " size " << sizet << " and inserted: " << n_parts_inserted << endl;

					// a bit of logic in case one of the particles has too little mass
					// to make up an individual particle near the target mass
					if (n_parts_inserted == 0)
					{
						n_parts_inserted = 1;
					}

					// now get the actual mass of each particle
					// this is exact: the mass in the cell will exactly match the
					// mass in the cell.
					particle_mass = mass_of_type_size/double(n_parts_inserted);
					soil_particle_mass = particle_mass*density_ratio;

					// get the particle surface area
					particle_surface_area = vpi.return_surface_area(type,sizet,particle_mass);
					soil_particle_surface_area = vpi.return_surface_area(type,sizet,soil_particle_mass);

					for (int part = 0; part<n_parts_inserted; part++)
					{
						// randomly generate the s and z locations
						ran_sl = (bin_edge_loc[bn+1] - bin_edge_loc[bn])*ran3(&seed) +
								  bin_edge_loc[bn];

						ran_zl = (insert_zone_top-insert_zone_bottom)*ran3(&seed) +
								 insert_zone_bottom;

						// now create and place the particles
						// for a given z and s location, you must calculate the depth and effective
						// d location
						// the effectve d location is in units g/cm^2, everything else is
						// in SI

						// first interpolate zeta and eta:
						interpolated_zeta = ((zeta[ h_node_ds[bn] ] - zeta[ h_node_us[bn] ])/
											dx_h[bn])*(ran_sl- s_us_h[bn]) + zeta[ h_node_us[bn] ];
						interpolated_eta = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
											dx_h[bn])*(ran_sl- s_us_h[bn]) + eta[ h_node_us[bn] ];

						interpolated_h = interpolated_zeta-interpolated_eta;
						d = interpolated_zeta-ran_zl;

						// if the depth is greater than the soil thickness,
						// the effective depth must take into account
						// the density change
						if (ran_zl<interpolated_eta)
						{
							d_rock = interpolated_eta-ran_zl;
							eff_d = 0.1*(rho_r*d_rock+rho_s*interpolated_h);
									// the 0.1 factor is to convert from
									// kg/m^2 to g/cm^2
							// create a particle
							insert_particle_mass = particle_mass;
							insert_particle_surface_area = particle_surface_area;

						}
						// if the depth is less than the soil thickness,
						// the effective depth only accounts for the density
						// of the soil
						else
						{
							eff_d = 0.1*(rho_s*d);
									// the 0.1 factor is to convert from
									// kg/m^2 to g/cm^2
							insert_particle_mass = soil_particle_mass;
							insert_particle_surface_area = soil_particle_surface_area;
						}

						// create a particle
						CRN_tParticle ins_part(type, sizet, ran_sl, d,eff_d, ran_zl,
												insert_particle_mass,insert_particle_surface_area);

						//cout << "s_loc: LINE 560: " << ins_part.getxLoc() << endl;

						// now update the initial cosmo concentrations
						ins_part.update_cosmo_conc_const(C_10Be, C_26Al, C_36Cl,
											 C_14C, C_21Ne, C_3He);

						if (ran_zl>=interpolated_eta)
						{
							ins_part.SoilAgeExpose();
						}
						particle_bins[bn].push_back(ins_part);

					}		// end particle inserion loop
				}			// end of logic: only insert particles if there is finite mass
			}				// end particle size loop
		}					// end particle type loop
	}						// end bin loop

	old_bottom_depth = new_bottom_depth;

	//for (int i = 0; i< n_bins; i++)
	//{
	//	cout << "LINE 588 CRN_tParticle_bins.cpp; n_parts of bin " << i << " = " 
  //         << particle_bins[i].size() << endl;
	//}
	return 1;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// NOTE: delta lowered is the same size as eta, each h node has one lowering thickness
// (thus affecting two bins each)
// this inserts particles into the simulation. One keeps track of soil production,
// and after some time T the total lowering of the soil bedrock boundary is calculated,
// this lowering is Delta_lowered = integral from T to T+Delta T of p dt where p is
// soil production in m/yr. Thus Delta_lowered is in metres, and a potentially different
// value exists for each h node.
// The particles are introduced in a box whose upper boundary is the start depth
// below the _current_ eta elevation, and whose bottom boundary is this top zone minus
// the distance Delta_lowered.
// part_conc is the number of particles per kg of rock
// rho_s and rho_r are used to calculate both the number of particles introduced and
// the effective depth of the particles.
//
// the pID is a particle identification number, each particle is tagged with a different integer
// so in plotting algorithms one particle at a time can be followed
// the function returns and integer, that is the value of the
// pID for the first particle of the next inserting sequence.
//
// this overloaded version of the function sets the scaling and parameter values for all of the particles
// CRN_muon_param_switch:
//		0: Granger
//		1: Schaller
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
int CRN_tParticle_bins::insert_particles(flowtube ft,
									vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double part_conc, vector<int> starting_pID, vector<double> starting_p_mfrac,
									double C_10Be, double C_26Al, double C_36Cl, double C_14C, double C_21Ne, double C_3He)
{

	vector<double> new_bottom_depth = old_bottom_depth;
	vector<double> zeta = ft.get_zeta();		// surface elevation (m)
	vector<double> eta = ft.get_eta();			// soil-bedrock boundary of the bin (m)

	double rho_r = ft.get_rho_r();
	double rho_s = ft.get_rho_s();

	long seed = time(NULL);             // seed for random number generator
	int n_parts_inserted;				// number of particles inserted into a bin
	double insert_zone_top;				// top of insertion zone (m)
	double insert_zone_bottom;			// bottom of insertion zone (m)
	double ran_sl;						// downslope distance of the randomly inserted particle (m)
	double ran_zl;						// the elevation of the randomly inserted particle (m)
	//double theta;						//
	int eta_node;						// the node index for eta, h, and zeta

	// interpolated hillslope properties
	double interpolated_zeta;			// surface elevation (m)
	double interpolated_eta;			// soil-bedrock boundary (m)
	double interpolated_h;				// soil thickness (m)
	double d_rock;						// depth of rock overlying particle (m)
	double d;							// depth of particle (m)
	double eff_d;						// effective depth of particle (g/cm^2)

	// during a timestep some particles are inserted into a zone in the bedrock
	// the number of particles is defined by the volume of the insertion zone times
	// the particle concnetration. The volume times the concentration will typicaly
	// be an integer, so we need a way to maintin some kind of rasonably constant
	// and unbiased particle concentration (i.e., avoiding rounding bias)
	// the approach is to get thh fractional part of the number of particles
	// and then select a random number from a uniform distribution on 0 to 1
	// , if the ranom number is less than the fraction a particle is generated,
	// if not no particle. The below variables are for making this calculation
	double n_parts_doub;
	double n_parts_int;
	double n_parts_frac;

	// this function loops through the bins, inserting particles based on how
	// much the soil/saprolite boundary has lowered and the particle concnetration,
	// which is given in particles/kg
	//
	// in general we consider the surfaces zeta and eta to be sloping, but for
	// calculation of the insertion zone the eta location is considered horizontal
	// this is because if the insertion zone was considered sloping then transient changes
	// in the bedrock might lead to ovelapping insertion zones. A horizontal insertion
	// zone ensures that particles are always inserted into 'fresh' rock
	for (int bn = 0; bn<n_bins; bn++)
	{
		// get the number of particles inserted
		n_parts_doub = A_bins[bn]*Delta_lowered[bn/2]*rho_r*part_conc;
		n_parts_frac = modf(n_parts_doub,&n_parts_int);
		n_parts_inserted = int(n_parts_int);
		if( ran3(&seed) <= n_parts_frac )
		 n_parts_inserted++;

		eta_node = bn/2;
		insert_zone_top = old_bottom_depth[eta_node];
		insert_zone_bottom = insert_zone_top - Delta_lowered[eta_node];
		new_bottom_depth[eta_node] = insert_zone_bottom;



		for (int part = 0; part<n_parts_inserted; part++)
		{
			// randomly generate the s and z locations
			ran_sl = (bin_edge_loc[bn+1] - bin_edge_loc[bn])*ran3(&seed) +
					  bin_edge_loc[bn];
			//cout << "bin number: " << bn << " us: " << bin_edge_loc[bn] << " ds: " << bin_edge_loc[bn+1]
			//     << " s_loc: " << ran_sl << endl;


			ran_zl = (insert_zone_top-insert_zone_bottom)*ran3(&seed) +
					 insert_zone_bottom;

			// now create and place the particles
			// for a given z and s location, you must calculate the depth and effective
			// d location
			// the effectve d location is in units g/cm^2, everything else is
			// in SI

			// first interpolate zeta and eta:
			interpolated_zeta = ((zeta[ h_node_ds[bn] ] - zeta[ h_node_us[bn] ])/
								dx_h[bn])*(ran_sl- s_us_h[bn]) + zeta[ h_node_us[bn] ];
			interpolated_eta = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*(ran_sl- s_us_h[bn]) + eta[ h_node_us[bn] ];

			interpolated_h = interpolated_zeta-interpolated_eta;
			d = interpolated_zeta-ran_zl;

			// if the depth is greater than the soil thickness,
			// the effective depth must take into account
			// the density change
			if (ran_zl<interpolated_eta)
			{
				d_rock = interpolated_eta-ran_zl;
				eff_d = 0.1*(rho_r*d_rock+rho_s*interpolated_h);
						// the 0.1 factor is to convert from
						// kg/m^2 to g/cm^2
			}
			// if the depth is less than the soil thickness,
			// the effective depth only accounts for the density
			// of the soil
			else
			{
				eff_d = 0.1*(rho_s*d);
						// the 0.1 factor is to convert from
						// kg/m^2 to g/cm^2
				//cout << "LINE 325 the inserted particle is in the soil! ";

			}

			// now create the particle:
			// determine the starting type
			double type_prob = ran3(&seed);
			int n_ptypes = starting_pID.size();
			int startType = starting_pID[n_ptypes-1];
			double cum_prob = 0.0;
			double last_prob = 0.0;
			for (int ptype_index = 0; ptype_index<n_ptypes; ptype_index++)
			{
				last_prob = cum_prob;
				cum_prob = cum_prob+starting_p_mfrac[ptype_index];

				if (type_prob >= last_prob && type_prob < cum_prob)
				{
					startType = starting_pID[ptype_index];
				}
			}

			// create a particle
			CRN_tParticle ins_part(startType, ran_sl,d,eff_d, ran_zl);

			// now update the initial cosmo concentrations
			ins_part.update_cosmo_conc_const(C_10Be, C_26Al, C_36Cl,
								 C_14C, C_21Ne, C_3He);

			if (ran_zl>=interpolated_eta)
			{
				ins_part.SoilAgeExpose();
			}
			//cout << "LINE 549 type: " << ins_part.getType() << " zeta: " << ins_part.get_zetaLoc()
			//	 << " d: " << ins_part.getdLoc() << " and eff_d: " << ins_part.geteffective_dLoc() << endl;
			particle_bins[bn].push_back(ins_part);

			startType++; 				// increment the particle type
		}
	}

	old_bottom_depth = new_bottom_depth;

	//for (int i = 0; i< n_bins; i++)
	//{
	//	cout << "LINE 347 CRN_tParticle_bins.cpp; n_parts of bin " << i << " = " << particle_bins[i].size() << endl;
	//}
	return 1;
}

// NOTE: delta lowered is the same size as eta, each h node has one lowering thickness
// (thus affecting two bins each)
// this inserts particles into the simulation. One keeps track of soil production,
// and after some time T the total lowering of the soil bedrock boundary is calculated,
// this lowering is Delta_lowered = integral from T to T+Delta T of p dt where p is
// soil production in m/yr. Thus Delta_lowered is in metres, and a potentially different
// value exists for each h node.
// The particles are introduced in a box whose upper boundary is the start depth
// below the _current_ eta elevation, and whose bottom boundary is this top zone minus
// the distance Delta_lowered.
// part_conc is the number of particles per kg of rock
// rho_s and rho_r are used to calculate both the number of particles introduced and
// the effective depth of the particles.
//
// the pID is a particle identification number, each particle is tagged with a different integer
// so in plotting algorithms one particle at a time can be followed
// the function returns and integer, that is the value of the
// pID for the first particle of the next inserting sequence.
//
// this overloaded version of the function sets the scaling and parameter values for all of the particles
// CRN_muon_param_switch:
//		0: Granger
//		1: Schaller
// it also initiates the profile with steady state concentrations
int CRN_tParticle_bins::insert_particles(flowtube ft, vector<double> Delta_lowered, vector<double>& old_bottom_depth,
									double part_conc, vector<int> starting_pID, vector<double> starting_p_mfrac,
									CRN_parameters& CRNp,
									double erosion_rate_in_mass_per_time_per_area)
{

	vector<double> new_bottom_depth = old_bottom_depth;
	vector<double> zeta = ft.get_zeta();		// surface elevation (m)
	vector<double> eta = ft.get_eta();			// soil-bedrock boundary of the bin (m)

	double rho_r = ft.get_rho_r();
	double rho_s = ft.get_rho_s();

	long seed = time(NULL);             // seed for random number generator
	int n_parts_inserted;				// number of particles inserted into a bin
	double insert_zone_top;				// top of insertion zone (m)
	double insert_zone_bottom;			// bottom of insertion zone (m)
	double ran_sl;						// downslope distance of the randomly inserted particle (m)
	double ran_zl;						// the elevation of the randomly inserted particle (m)
	//double theta;						//
	int eta_node;						// teh node index for eta, h, and zeta

	// interpolated hillslope properties
	double interpolated_zeta;			// surface elevation (m)
	double interpolated_eta;			// soil-bedrock boundary (m)
	double interpolated_h;				// soil thickness (m)
	double d_rock;						// depth of rock overlying particle (m)
	double d;							// depth of particle (m)
	double eff_d;						// effective depth of particle (g/cm^2)

	// during a timestep some particles are inserted into a zone in the bedrock
	// the number of particles is defined by the volume of the insertion zone times
	// the particle concnetration. The volume times the concentration will typicaly
	// be an integer, so we need a way to maintin some kind of rasonably constant
	// and unbiased particle concentration (i.e., avoiding rounding bias)
	// the approach is to get thh fractional part of the number of particles
	// and then select a random number from a uniform distribution on 0 to 1
	// , if the ranom number is less than the fraction a particle is generated,
	// if not no particle. The below variables are for making this calculation
	double n_parts_doub;
	double n_parts_int;
	double n_parts_frac;

	// this function loops through the bins, inserting particles based on how
	// much the soil/saprolite boundary has lowered and the particle concnetration,
	// which is given in particles/kg
	//
	// in general we consider the surfaces zeta and eta to be sloping, but for
	// calculation of the insertion zone the eta location is considered horizontal
	// this is because if the insertion zone was considered sloping then transient changes
	// in the bedrock might lead to ovelapping insertion zones. A horizontal insertion
	// zone ensures that particles are always inserted into 'fresh' rock
	for (int bn = 0; bn<n_bins; bn++)
	{
		// get the number of particles inserted
		n_parts_doub = A_bins[bn]*Delta_lowered[bn/2]*rho_r*part_conc;
		n_parts_frac = modf(n_parts_doub,&n_parts_int);
		n_parts_inserted = int(n_parts_int);
		if( ran3(&seed) <= n_parts_frac )
		 n_parts_inserted++;

		eta_node = bn/2;
		insert_zone_top = old_bottom_depth[eta_node];
		insert_zone_bottom = insert_zone_top - Delta_lowered[eta_node];
		new_bottom_depth[eta_node] = insert_zone_bottom;

		for (int part = 0; part<n_parts_inserted; part++)
		{
			// randomly generate the s and z locations
			ran_sl = (bin_edge_loc[bn+1] - bin_edge_loc[bn])*ran3(&seed) +
					  bin_edge_loc[bn];

			ran_zl = (insert_zone_top-insert_zone_bottom)*ran3(&seed) +
					 insert_zone_bottom;

			// now create and place the particles
			// for a given z and s location, you must calculate the depth and effective
			// d location
			// the effectve d location is in units g/cm^2, everything else is
			// in SI

			// first interpolate zeta and eta:
			interpolated_zeta = ((zeta[ h_node_ds[bn] ] - zeta[ h_node_us[bn] ])/
								dx_h[bn])*(ran_sl- s_us_h[bn]) + zeta[ h_node_us[bn] ];
			interpolated_eta = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*(ran_sl- s_us_h[bn]) + eta[ h_node_us[bn] ];

			interpolated_h = interpolated_zeta-interpolated_eta;
			d = interpolated_zeta-ran_zl;

			// if the depth is greater than the soil thickness,
			// the effective depth must take into account
			// the density change
			if (ran_zl<interpolated_eta)
			{
				d_rock = interpolated_eta-ran_zl;
				eff_d = 0.1*(rho_r*d_rock+rho_s*interpolated_h);
						// the 0.1 factor is to convert from
						// kg/m^2 to g/cm^2
			}
			// if the depth is less than the soil thickness,
			// the effective depth only accounts for the density
			// of the soil
			else
			{
				eff_d = 0.1*(rho_s*d);
						// the 0.1 factor is to convert from
						// kg/m^2 to g/cm^2
				//cout << "LINE 325 the inserted particle is in the soil! ";

			}

			// now create the particle:
			// determine the starting type
			double type_prob = ran3(&seed);
			int n_ptypes = starting_pID.size();
			int startType = starting_pID[n_ptypes-1];
			double cum_prob = 0.0;
			double last_prob = 0.0;
			for (int ptype_index = 0; ptype_index<n_ptypes; ptype_index++)
			{
				last_prob = cum_prob;
				cum_prob = cum_prob+starting_p_mfrac[ptype_index];

				if (type_prob >= last_prob && type_prob < cum_prob)
				{
					startType = starting_pID[ptype_index];
				}
			}

			// create a particle
			CRN_tParticle ins_part(startType, ran_sl,d,eff_d, ran_zl);

			// now update the initial cosmo concentrations
			ins_part.update_all_CRN_SSfull(erosion_rate_in_mass_per_time_per_area, CRNp);

			if (ran_zl>=interpolated_eta)
			{
				ins_part.SoilAgeExpose();
			}
			particle_bins[bn].push_back(ins_part);
			startType++; 				// increment the particle type

		}

	}
	old_bottom_depth = new_bottom_depth;

	//for (int i = 0; i< n_bins; i++)
	//{
	//	cout << "LINE 347 CRN_tParticle_bins.cpp; n_parts of bin " << i << " = " << particle_bins[i].size() << endl;
	//}
	return 1;
}


// this updates fallout radionuclides
// NOTE!!!!
// the units of M_supply are in atoms/cm^2
// and conc is in atoms per gram
// but conversion of rho and deltad goes on
// withing the CRN_tParticle object
// the units of k_f10Be here are cm^2/g
void CRN_tParticle_bins::update_fallout_10Be_bins(double dt, double M_supply_surface,
					double rho_s, double k_f10Be, double deltad, CRN_parameters& CRNp)
{
	list<CRN_tParticle>::iterator part_iter;	// list iterator

	for (int bn = 0; bn<n_bins; bn++)
	{
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			(*part_iter).update_fallout10Be_simple_density(dt, M_supply_surface,
					      rho_s, k_f10Be, deltad, CRNp);
		    part_iter++;
		}
	}
}

// overloaded function
void CRN_tParticle_bins::update_fallout_10Be_bins(double dt, double M_supply_surface,
					double rho_s, double k1_f10Be, double k2_f10Be, double chi_f10Be,
					double deltad, CRN_parameters& CRNp)
{
	list<CRN_tParticle>::iterator part_iter;	// list iterator

	for (int bn = 0; bn<n_bins; bn++)
	{
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			(*part_iter).update_fallout10Be_simple_density_2exp(dt, M_supply_surface,
					      rho_s, k1_f10Be, k2_f10Be, chi_f10Be, deltad, CRNp);
		    part_iter++;
		}
	}
}

// this function takes several data elements from the flowtube object
// that have been
vector< list<CRN_tParticle> > CRN_tParticle_bins::particle_motion(double dt, flowtube ft,
										double Omega,double vert_vel_fluctuating,
										double horiz_vel_fluctuating,
										const int CRN_switch,
										CRN_parameters& CRNp)
{

	vector< list<CRN_tParticle> > eroded_bins(n_bins+1);
	vector< list<CRN_tParticle> > moved_bins(n_bins);

	long seed = time(NULL);             // seed for random number generator

	// soil and rock densities
	double rho_s = ft.get_rho_s();
	double rho_r = ft.get_rho_r();

	// hillslope properties
	vector<double> h = ft.get_h();
	vector<double> intermediate_h = ft.get_intermediate_h();
	vector<double> old_h = ft.get_old_h();
	vector<double> zeta = ft.get_zeta();
	vector<double> intermediate_zeta = ft.get_intermediate_zeta();
	vector<double> old_zeta = ft.get_old_zeta();
	vector<double> eta = ft.get_eta();
	vector<double> old_eta = ft.get_old_eta();
	int n_h_nodes = h.size();
	vector<double> Mass_Flux = ft.get_Mass_Flux();
	vector<double> fluff = ft.get_fluff();

	// local parameters (these are all linearly interpolated
	double z_loc;
	double s_loc;
	double eta_new_local;
	double eta_old_local;
	double h_old_local;
	double zeta_new_local;
	double zeta_old_local;
	double intermediate_zeta_local;
	double intermediate_eta_local;
	double fluff_local;
	double h_local;
	double F_local;
	double b_local;

	//double test_z_loc;					// for bug checking
	//double test_s_loc;					// for bug checking
	//double test_d_loc;					// for bug checking

	double d_rock;						// depth of rock overlying particle (m)
	double d;							// depth of particle (m)
	double old_eff_d;
	double eff_d;						// effective depth of particle (g/cm^2)
	double delta_eff_d;					// change in effective depth
	double eff_erosion_rate;

	double psi;							// a dimensionless number that
										// is used to represent the fraction
										// elevation in either the fluff zone
										// or the soil column

	double v_avg;						// the average downslope velocity of a particle
	double v_ds;						// downslope velocity
	int post_move_bn = 0;				// the bin number of a particle moving downslope
	double reflect_distance;
	double type;

	list<CRN_tParticle>::iterator part_iter;	// list iterator
	list<CRN_tParticle>::iterator remove_iter;	// iterator that points to removed
												// particles


  int n_eroded = 0;
  
	// loop through all the bins
	for (int bn = 0; bn< n_bins; bn++)
	{
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			// get the z_location of the particle
			z_loc = (*part_iter).get_zetaLoc();
			s_loc = (*part_iter).getxLoc();
			type  = (*part_iter).getType();

			// get the old effective depth
			old_eff_d = (*part_iter).geteffective_dLoc();

			//cout << "LINE 843 Type: " << type << " z_loc is: " << z_loc << " and old_eff_d is: " << old_eff_d << endl;



			// there are three options for the particle, all
			// based on its elevation relative to the
			// soil-bedrock interface
			// 1. In bedrock the entire time: z_loc < eta_new
			// 2. In fluff zone: z_loc > eta_new && < eta_old
			// 3. In soil: z_loc > eta_old
			// To determine this we must know both
			// eta_new_local and eta_old_local
			eta_new_local = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + eta[ h_node_us[bn] ];
			eta_old_local = ((old_eta[ h_node_ds[bn] ] - old_eta[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + old_eta[ h_node_us[bn] ];



			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			// PARTICLES THAT REMAIN IN BEDROCK
			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			// if the particle is in bedrock throughout the timestep, there is no
			// update of the z and s locations, but we must update the d location,
			// as well as the effective d location
			// to do this we need the new zeta location
			if (z_loc < eta_new_local)
			{
				// bin number does not change
				post_move_bn = bn;

				// get the new local elevation of the soil surface
				zeta_new_local = ((zeta[ h_node_ds[bn] ] - zeta[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + zeta[ h_node_us[bn] ];

				// now calculate the new depth and d_effective
				h_local = zeta_new_local-eta_new_local;
				d = zeta_new_local-z_loc;

				d_rock = eta_new_local-z_loc;
				eff_d = 0.1*(rho_r*d_rock+rho_s*h_local);
							// the 0.1 factor is to convert from
							// kg/m^2 to g/cm^2

				// now get the erosion rate
				delta_eff_d = eff_d - old_eff_d;
				eff_erosion_rate = delta_eff_d/dt;

				// update the depth and effective depth of the particle
				(*part_iter).update_depths(d, eff_d);

				// update the CRN concentrations
				switch (CRN_switch)
				{
					case 0:
						break;
					case 1:
						(*part_iter).update_all_CRN(dt, eff_erosion_rate, CRNp);
						break;
					case 2:
						(*part_iter).update_all_CRN_neutron_only(dt, eff_erosion_rate, CRNp);
						break;
					case 3:
						(*part_iter).update_10Be_conc(dt,eff_erosion_rate, CRNp);
						break;
					default:
						(*part_iter).update_all_CRN_neutron_only(dt, eff_erosion_rate, CRNp);

				}
			}
			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			// PARTICLES IN FLUFF ZONE NO DOWNSLOPE MOVEMENT
			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			// if the particle is in the fluff zone, it then gets'fluffed' as a
			// function of how elevated in the fluff zone it is.
			// these particles do not move downslope and are not subject to
			// random vertical motions
			else if (z_loc < eta_old_local && z_loc >= eta_new_local)
			{
				// bin number does not change
				post_move_bn = bn;

				zeta_new_local = ((zeta[ h_node_ds[bn] ] - zeta[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + zeta[ h_node_us[bn] ];

				fluff_local = ((fluff[ h_node_ds[bn] ] - fluff[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + fluff[ h_node_us[bn] ];

				// update the z location within the fluff zone
				psi = (z_loc-eta_new_local)/(eta_old_local-eta_new_local);
				z_loc = z_loc+psi*fluff_local;

				// now update the effective depth
				// note no rock density is needed because the particle is
				// now in the soil
				d = zeta_new_local-z_loc;

				// d could be less than zero if there is surface
				// erosion, so set to zero if this is the case
				if (d<0)
				{
					eff_d = 0;
				}
				else
				{
					eff_d = 0.1*(rho_s*d);
							// the 0.1 factor is to convert from
							// kg/m^2 to g/cm^2
				}

				// now get the erosion rate
				delta_eff_d = eff_d - old_eff_d;
				eff_erosion_rate = delta_eff_d/dt;

				// update the z location
				(*part_iter).update_zetaLoc(z_loc);

				// update the depth and effective depth of the particle
				(*part_iter).update_depths(d, eff_d);

				// update the CRN concentrations
				switch (CRN_switch)
				{
					case 0:
						break;
					case 1:
						(*part_iter).update_all_CRN(dt, eff_erosion_rate, CRNp);
						break;
					case 2:
						(*part_iter).update_all_CRN_neutron_only(dt, eff_erosion_rate, CRNp);
						break;
					case 3:
						(*part_iter).update_10Be_conc(dt,eff_erosion_rate, CRNp);
						break;
					default:
						(*part_iter).update_all_CRN_neutron_only(dt, eff_erosion_rate, CRNp);

				}

				// test if the age has been reset
				if ((*part_iter).getAge()<0)
					(*part_iter).SoilAgeExpose();


				// the update the age
				(*part_iter).incrementAge(dt);
			}
			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			// PARTICLES MOVING DOWNSLOPE: PARtICLES IN SOIL
			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			// now look at particles that are in the soil column
			// to calculate the downslope velocity of a particle
			// in the soil zone,
			else if (z_loc > eta_old_local)
			{
				int move_ds = 1;		// a switch, if this == 1 then
										// the particle moves downslope
				// see if the particle is active:
				// if the activity, Omega == 1,
				// all particles are active, whereas is it is less than
				// one only some particles are active
				if (Omega < 1)
				{
					// if the random number is less than omega,
					// the particle moves, if not it doesn't move
					if (ran3(&seed) > Omega)
					{
						move_ds = 0;
					}
				}

				// now we have two logic loops, if the particle moves downslope
				// and if it does not move downslope
				if (move_ds ==1)
				{
					// the aveage velocity of a particle is given by
					// v_avg = F/(Omega*rho_s*h*b)
					// these parameters are interpolated locally, using the
					// values at the beginning of the timestep
					h_old_local = ((old_h[ h_node_ds[bn] ] - old_h[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + old_h[ h_node_us[bn] ];
					b_local = (Slope_b[bn])*(s_loc- s_us_b[bn]) + b_us[bn];
					F_local = ((Mass_Flux[ b_node_ds[bn] ] - Mass_Flux[ b_node_us[bn] ])/
								dx_b[bn])*(s_loc- s_us_b[bn]) + Mass_Flux[ b_node_us[bn] ];

					// you also have to get the psi value in order to place the particle
					// in the new location to do this you need the
					// old zeta value
					zeta_old_local = ((old_zeta[ h_node_ds[bn] ] - old_zeta[ h_node_us[bn] ])/
									dx_h[bn])*(s_loc- s_us_h[bn]) + old_zeta[ h_node_us[bn] ];

					// get the fractional elevation of the particle in the soil
					// this will be preserved after downslope movement
					psi = (z_loc-eta_old_local)/(zeta_old_local-eta_old_local);

					// get the average downslope velocity
					v_avg = F_local/(Omega*b_local*h_old_local*rho_s);

					//
					// NOT IMPLEMENTED YET: VARY THE VELOCITY RANDOMLY DOWNSLOPE
					//
					// for now the velocity is just v_avg
					v_ds = v_avg;

					// now move the particle
					s_loc+= v_ds*dt;

					// check to see if it is in the same bin,
					// if not assign it to a new bin

					//cout << "LINE 601 upslope bn: " << bin_edge_loc[bn]
					//	 << " and downslope: " << bin_edge_loc[bn+1]
					//	 << " and s loc: " << s_loc << endl;


					if (s_loc >= bin_edge_loc[bn] && s_loc <= bin_edge_loc[bn+1])
					{
						post_move_bn = bn;
					}
					else if (s_loc < bin_edge_loc[bn])
					{
						post_move_bn = bn-1;
						//cout << "LINE 622, particle moved into upslope bin!\n";
					}
					else if (s_loc > bin_edge_loc[bn+1])
					{
						post_move_bn = bn+1;
					}

					// if the particle has moved above the upslope boundary, reflect
					// it back into the upslope bin
					if (post_move_bn < 0)
					{
						reflect_distance = bin_edge_loc[0]-s_loc;
						s_loc = bin_edge_loc[0]+reflect_distance;
						post_move_bn = 0;
					}

					// now update the particle properties if it has not left the
					// model domain
					if ( post_move_bn < n_bins)
					{
						// now get the new z elevation by interpolating the new distance
						intermediate_zeta_local = ((intermediate_zeta[ h_node_ds[post_move_bn] ]
												- intermediate_zeta[ h_node_us[post_move_bn] ])/
												dx_h[post_move_bn])*(s_loc- s_us_h[post_move_bn])
												+ intermediate_zeta[ h_node_us[post_move_bn] ];
						intermediate_eta_local = ((old_eta[ h_node_ds[post_move_bn] ]
											- old_eta[ h_node_us[post_move_bn] ])/
											dx_h[post_move_bn])*(s_loc- s_us_h[post_move_bn])
											+ old_eta[ h_node_us[post_move_bn] ];

						// reproduce the original psi value at the new location
						z_loc = psi*(intermediate_zeta_local-intermediate_eta_local)
									+ intermediate_eta_local;

						// now fluff the particle
						fluff_local = ((fluff[ h_node_ds[post_move_bn] ]
									- fluff[ h_node_us[post_move_bn] ])/
									dx_h[post_move_bn])*(s_loc- s_us_h[post_move_bn])
									+ fluff[ h_node_us[post_move_bn] ];
						z_loc += fluff_local;

						// now vertically mix particle
						// right now the particle moves at the same velocity
						// either up or down
						if (ran3(&seed) > 0.5)
						{
							z_loc += vert_vel_fluctuating*dt;
						}
						else
						{
							z_loc -= vert_vel_fluctuating*dt;
						}

						// now reflect the particle if it is outside the
						// soil
						if (z_loc > intermediate_zeta_local)
						{
							reflect_distance = z_loc - intermediate_zeta_local;
							z_loc = intermediate_zeta_local- reflect_distance;

							// update the OSL_age
							(*part_iter).OSLexpose();

						}
						else if (z_loc < intermediate_eta_local)
						{
							reflect_distance = intermediate_eta_local-z_loc;
							z_loc = intermediate_eta_local+ reflect_distance;
						}


						// now get the new zeta location (used for depth calucaltions)
						zeta_new_local = ((zeta[ h_node_ds[post_move_bn] ]
												- zeta[ h_node_us[post_move_bn] ])/
												dx_h[post_move_bn])*(s_loc- s_us_h[post_move_bn])
												+ zeta[ h_node_us[post_move_bn] ];


					}
					else
					{
						// logic if particle has left domain
						zeta_new_local = zeta[n_h_nodes-1];
					}
				} // !! end logic if particle moved
				// note at this stage the particle has not moved into the next bin of the particle
				// list, even if it has in fact moved



				// now the logic if the particle didn't move downslope
				else
				{
					// bin number does not change
					post_move_bn = bn;

					zeta_new_local = ((zeta[ h_node_ds[bn] ] - zeta[ h_node_us[bn] ])/
									dx_h[bn])*(s_loc- s_us_h[bn]) + zeta[ h_node_us[bn] ];

					fluff_local = ((fluff[ h_node_ds[bn] ] - fluff[ h_node_us[bn] ])/
									dx_h[bn])*(s_loc- s_us_h[bn]) + fluff[ h_node_us[bn] ];

					// update the z location within the fluff zone
					z_loc = z_loc+psi*fluff_local;

					// now vertically mix particle
					// right now the particle moves at the same velocity
					// either up or down
					if (ran3(&seed) > 0.5)
					{
						z_loc += vert_vel_fluctuating*dt;
					}
					else
					{
						z_loc -= vert_vel_fluctuating*dt;
					}
					// now reflect the particle if it is outside the
					// soil
					if (z_loc > intermediate_zeta_local)
					{
						reflect_distance = z_loc - intermediate_zeta_local;
						z_loc = intermediate_zeta_local- reflect_distance;

						// update the OSL_age
						(*part_iter).OSLexpose();

					}
					else if (z_loc < intermediate_eta_local)
					{
						reflect_distance = intermediate_eta_local-z_loc;
						z_loc = intermediate_eta_local+ reflect_distance;
					}

				}

				// now update the particle properties
				// now update the effective depth
				// note no rock density is needed because the particle is
				// now in the soil
				d = zeta_new_local-z_loc;

				// d could be less than zero if there is surface
				// erosion, so set to zero if this is the case
				if (d<0)
				{
					eff_d = 0;
				}
				else
				{
					eff_d = 0.1*(rho_s*d);
							// the 0.1 factor is to convert from
							// kg/m^2 to g/cm^2
				}

				// now get the erosion rate
				delta_eff_d = eff_d - old_eff_d;
				eff_erosion_rate = delta_eff_d/dt;

				// update the z location
				(*part_iter).update_zetaLoc(z_loc);

				// update s_loc
				(*part_iter).update_xLoc(s_loc);

				// update the depth and effective depth of the particle
				(*part_iter).update_depths(d, eff_d);

				// update the CRN concentrations
				switch (CRN_switch)
				{
					case 0:
						break;
					case 1:
						(*part_iter).update_all_CRN(dt, eff_erosion_rate, CRNp);
						break;
					case 2:
						(*part_iter).update_all_CRN_neutron_only(dt, eff_erosion_rate, CRNp);
						break;
					case 3:
						(*part_iter).update_10Be_conc(dt,eff_erosion_rate, CRNp);
						break;
					default:
						(*part_iter).update_all_CRN_neutron_only(dt, eff_erosion_rate, CRNp);

				}

				// the update the age
				(*part_iter).incrementAge(dt);


			}			// end logic for if particle is in the soil
			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
			//

			// now see if the particle has been eroded from the
			// surface layer or out of the particle domain
			if (s_loc < bin_edge_loc[n_bins])
			{
				if (post_move_bn >= n_bins)
				{
					cout << "error LINE 805 the updated bin is outside of the domain!";
					exit(1);
				}

				if (d < 0)
				{
					// the bin number of the eroded bins
					eroded_bins[bn].push_back(*part_iter);
					remove_iter = part_iter;
					part_iter++;
					particle_bins[bn].erase( remove_iter );
					n_eroded++;
				}
				else if (post_move_bn != bn)
				{
					//cout << endl << "LINE 819 particle moved into downslope bin!" << endl;
					//cout << "particle number is: " << (*part_iter).getType() << endl;
					//cout << "s_loc: " << s_loc << " bin number: "
					//     << bn << " new bn: " << post_move_bn << endl;
					//cout << "upslope bn: " << bin_edge_loc[post_move_bn]
					//     << " and downslope: " << bin_edge_loc[post_move_bn+1] << endl;
					//exit(1);

					// now remove the particle if it is no longer in the same bin and place
					// it in the move particle vector
					//cout << "moved bins: " << moved_bins.size() << endl;
					moved_bins[post_move_bn].push_back(*part_iter);
					remove_iter = part_iter;
					part_iter++;
					particle_bins[bn].erase( remove_iter );
				}
				else
				{
					part_iter++;
				}
			}
			// these particles have been eroded from the downslope boundary
			else
			{
				eroded_bins[n_bins].push_back(*part_iter);
				remove_iter = part_iter;
				part_iter++;
				particle_bins[bn].erase( remove_iter );
			}

		}				// end loop through list
	}					// end loop through bins


	// now loop through the moved bins
	// if there are any moved particles, take them from
	// the moved_bins list and place them in
	// the normal particle bins
	int n_parts_in_list;
	for (int bn = 0; bn < n_bins; bn++)
	{
		n_parts_in_list = moved_bins[bn].size();
		if (n_parts_in_list > 0)
		{
			part_iter = moved_bins[bn].begin();
			while (part_iter != moved_bins[bn].end())
			{
				particle_bins[bn].push_back(*part_iter);
				part_iter++;
			}
		}
	}


  //cout << "Number of particles eroded: " << n_eroded << endl;
	return eroded_bins;
}

void CRN_tParticle_bins::update_CRN_conc_const(double C_10Be, double C_26Al,
										double C_36Cl, double C_14C,
										 double C_21Ne, double C_3He)
{
	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// loop through all the bins
	for (int bn = 0; bn< n_bins; bn++)
	{
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			// update the cosmo concentrations
			(*part_iter).update_cosmo_conc_const(C_10Be, C_26Al, C_36Cl,
								  C_14C, C_21Ne, C_3He);
			part_iter++;
		}
	}

}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
// sampling functions

// this function returns the bin number given a x location
int CRN_tParticle_bins::get_bin_number_of_sloc(double s_loc)
{
	int bin_number;
	//cout << "n_bins is: " << n_bins << endl;
	//for (int i = 0; i<n_bins; i++)
	//{
	//	cout << "bin number " << i << " us edge = " << bin_edge_loc[i] << " ds edge: " << bin_edge_loc[i+1] << endl;
	//}

	if (s_loc < bin_edge_loc[0])
	{
		cout << "Error: you have selected a s location that is upslope of the first node!" << endl
		     << "Assigning bin number 0" << endl;
		bin_number = 0;
	}
	else if (s_loc > bin_edge_loc[n_bins])
	{
		cout << "Error, you have selected an s location that is downslope of the last node!" << endl
		     << "Assigning last bin" << endl;
		bin_number = n_bins+1;
	}
	else
	{
		//cout << endl << endl;
		int bn = 0;
		// loop through bins until you get the right bin location
		while (bn < n_bins)
		{
			// if the s location of the downslope bin edge is less
			// than the selected s location, then stop looping and take the previous
			// bin
			//cout << "bn: " << bn << " us bin: " << bin_edge_loc[bn+1] << " and s_loc " << s_loc << endl;
			if (bin_edge_loc[bn+1] > s_loc && bin_edge_loc[bn] <= s_loc)
			{
				//cout << "triggered catch; bin: " << bn << endl;
				bin_number = bn;
				bn = n_bins;
			}
			else
			{
				bn++;
			}
		}
	}
	return bin_number;
}





//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
// printing functions
void CRN_tParticle_bins::print_particle_stats(double t_ime, ofstream& particle_out)
{
	double s_loc;				// distance downslope of the particle
	double z_loc;				// elevation of the particle
	double d_loc;				// depth of the particle
	double pID;					// the particle identification number

	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// loop through all the bins
	for (int bn = 0; bn< n_bins; bn++)
	{
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			// get the z_location of the particle
			pID = (*part_iter).getType();
			z_loc = (*part_iter).get_zetaLoc();
			s_loc = (*part_iter).getxLoc();
			d_loc = (*part_iter).getdLoc();

			particle_out << t_ime << " " << bn << " " << pID << " " << z_loc
						 << " " << s_loc << " " << d_loc << endl;
			part_iter++;
		}
	}

}

void CRN_tParticle_bins::print_particle_stats(double t_ime, flowtube ft,
								ofstream& particle_out)
{
	double s_loc;				// distance downslope of the particle
	double z_loc;				// elevation of the particle
	double d_loc;				// depth of the particle
	double pID;					// the particle identification number
	double pAge;
	double pOSLage;
	double C10Be;
	double C14C;
	double C21Ne;

	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// hillslope properties
	vector<double> eta = ft.get_eta();
	double eta_new_local;
	// loop through all the bins
	for (int bn = 0; bn< n_bins; bn++)
	{
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			// get the z_location of the particle
			pID = (*part_iter).getType();
			z_loc = (*part_iter).get_zetaLoc();
			s_loc = (*part_iter).getxLoc();
			d_loc = (*part_iter).getdLoc();
			pAge = (*part_iter).getAge();
			pOSLage = (*part_iter).getOSLage();
			C10Be = (*part_iter).getConc_10Be();
			C14C = (*part_iter).getConc_14C();
			C21Ne = (*part_iter).getConc_21Ne();

			particle_out << t_ime << " " << bn << " " << pID << " " << z_loc
						 << " " << s_loc << " " << d_loc << " ";
			eta_new_local = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + eta[ h_node_us[bn] ];
			if (z_loc > eta_new_local)
			{
				particle_out << 1 << " ";
			}
			else
			{
				particle_out << 0 << " ";
			}
			particle_out << pAge << " " << pOSLage << " " << C10Be << " "
			             << C14C << " " << C21Ne << endl;
			part_iter++;
		}
	}

}


void CRN_tParticle_bins::print_particle_stats_soil(double t_ime, flowtube ft,
								ofstream& particle_out)
{
	double s_loc;				// distance downslope of the particle
	double z_loc;				// elevation of the particle
	double d_loc;				// depth of the particle
	double pID;					// the particle identification number
	double pAge;
	double pOSLage;
	double C10Be;
	double C14C;
	double C21Ne;

	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// hillslope properties
	vector<double> eta = ft.get_eta();
	double eta_new_local;
	// loop through all the bins
	for (int bn = 0; bn< n_bins; bn++)
	{
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			// get the z_location of the particle
			pID = (*part_iter).getType();
			z_loc = (*part_iter).get_zetaLoc();
			s_loc = (*part_iter).getxLoc();
			d_loc = (*part_iter).getdLoc();
			pAge = (*part_iter).getAge();
			pOSLage = (*part_iter).getOSLage();
			C10Be = (*part_iter).getConc_10Be();
			C14C = (*part_iter).getConc_14C();
			C21Ne = (*part_iter).getConc_21Ne();


			eta_new_local = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + eta[ h_node_us[bn] ];
			if (z_loc > eta_new_local)
			{
				particle_out << t_ime << " " << bn << " " << pID << " " << z_loc
						 << " " << s_loc << " " << d_loc << " ";
				particle_out << 1 << " ";
				particle_out << pAge << " " << pOSLage << " " << C10Be << " "
			             << C14C << " " << C21Ne << endl;
			}

			part_iter++;
		}
	}

}

void CRN_tParticle_bins::print_eroded_stats(double t_ime,
							vector< list<CRN_tParticle> > eroded_list_vec,
							flowtube ft,
								ofstream& particle_out)
{
	double s_loc;				// distance downslope of the particle
	double z_loc;				// elevation of the particle
	double d_loc;				// depth of the particle
	double pID;					// the particle identification number
	double pAge;
	double pOSLage;
	double C10Be;
	double C14C;
	double C21Ne;


	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// hillslope properties
	vector<double> eta = ft.get_eta();
	// loop through all the bins
	for (int bn = 0; bn<= n_bins; bn++)
	{

		// now loop through each particle in the bin
		part_iter = eroded_list_vec[bn].begin();
		while (part_iter != eroded_list_vec[bn].end())
		{
			// get the z_location of the particle
			pID = (*part_iter).getType();
			z_loc = (*part_iter).get_zetaLoc();
			s_loc = (*part_iter).getxLoc();
			d_loc = (*part_iter).getdLoc();
			pAge = (*part_iter).getAge();
			pOSLage = (*part_iter).getOSLage();
			C10Be = (*part_iter).getConc_10Be();
			C14C = (*part_iter).getConc_14C();
			C21Ne = (*part_iter).getConc_21Ne();

			particle_out << t_ime << " " << bn << " " << pID << " " << "-99"
						 << " " << s_loc << " " << "-99" << " ";

			particle_out << "-99 ";	// this line is a placeholder where in the normal
									// particle stats printing it tells you if
			particle_out << pAge << " " << pOSLage << " " << C10Be << " "
			             << C14C << " " << C21Ne << endl;
			part_iter++;
		}
	}

}


void CRN_tParticle_bins::print_particle_stats_vtk(double t_ime, flowtube ft,
								string vtk_fname)
{
	vector<double> s_loc;				// distance downslope of the particle
	vector<double> z_loc;				// elevation of the particle
	vector<double> d_loc;				// depth of the particle
	vector<int> pID;					// the particle identification number
	vector<double> pAge;
	vector<double> pOSLage;
	vector<double> C10Be;
	vector<double> C14C;
	vector<double> C21Ne;
	vector<int> soil_switch;

	list<CRN_tParticle>::iterator part_iter;	// list iterator
	string time_bit = itoa( int(t_ime+0.5) );
	string vtk_ext = ".vtk";
	string fname = vtk_fname+time_bit+vtk_ext;
	//cout << "LINE 1149 CRN_tparticle_ bins, vtk filename is: " << fname << " and time: " << t_ime << endl;
	ofstream vtk_out;
	vtk_out.open(fname.c_str());

	// hillslope properties
	vector<double> eta = ft.get_eta();
	double eta_new_local;
	// loop through all the bins
	for (int bn = 0; bn< n_bins-1; bn++)
	{
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		int counter = 0;
		while (part_iter != particle_bins[bn].end())
		{
			// get the data from each particle
			s_loc.push_back( (*part_iter).getxLoc() );
			z_loc.push_back( (*part_iter).get_zetaLoc() );
			d_loc.push_back( (*part_iter).getdLoc() );
			pID.push_back( (*part_iter).getType() );
			pAge.push_back( (*part_iter).getAge() );
			pOSLage.push_back( (*part_iter).getOSLage() );
			C10Be.push_back( (*part_iter).getConc_10Be() );
			C14C.push_back( (*part_iter).getConc_14C() );
			C21Ne.push_back( (*part_iter).getConc_21Ne() );

			eta_new_local = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*( (*part_iter).getxLoc()- s_us_h[bn]) + eta[ h_node_us[bn] ];
			if ( (*part_iter).get_zetaLoc() > eta_new_local)
			{
				soil_switch.push_back(1);
			}
			else
			{
				soil_switch.push_back(0);
			}

			part_iter++;
			counter ++;
		}
		//cout << "bin number " << bn << " and counter is: " << counter << endl;
		//cout << "n_parts_total: " << d_loc.size() << endl;
	}

	// find the number of particles
	int n_parts = d_loc.size();
	//cout << "d_loc size: " << n_parts << endl;

	vtk_out << "# vtk DataFile Version 2.0" << endl << "Unstructured Grid Ptrack"
	        << endl << "ASCII" << endl << endl << "DATASET UNSTRUCTURED_GRID" << endl
	        << "POINTS " << n_parts << " float" << endl;

	for (int i = 0; i< n_parts; i++)
	{
		vtk_out << s_loc[i] << " " << z_loc[i] << " 0.0" <<endl;
	}

	vtk_out << endl << "POINT_DATA "<<n_parts << endl << "SCALARS Particle_age float 1"
	        << endl << "LOOKUP_TABLE default" << endl;

	for (int i = 0; i< n_parts; i++)
	{
		vtk_out << pAge[i] <<endl;
	}

	vtk_out << "SCALARS Be10_conc float 1"
	        << endl << "LOOKUP_TABLE default" << endl;

	for (int i = 0; i< n_parts; i++)
	{
		vtk_out << C10Be[i] <<endl;
	}

	vtk_out.close();
}


void CRN_tParticle_bins::print_age_cdfpdf_bulk(double t_ime, double max_age, int n_spacings,
							double K_times_D, double D, double sigma,
						  flowtube& ft,
						  ofstream& cdf_out,
						  ofstream& pdf_out)
{
	// hillslope properties
	vector<double> eta = ft.get_eta();

	double eta_new_local;

	double age_bin_width = max_age/double(n_spacings);
  	vector<double> bins(n_spacings);
  	vector<double> age_pdf(n_spacings,0.0);
  	vector<double> age_cdf(n_spacings,0.0);
  	double total_mass = 0;
  	double total_mass_weighted_age = 0;
  	double mean_age;

  	for (int i = 0; i<n_spacings; i++)
  	{
  		bins[i] = double(i)*age_bin_width+age_bin_width/2;
	}

	//cout << "LINE 954 looped bins" << endl;

	double s_loc;				// distance downslope of the particle
	double z_loc;				// elevation of the particle
	double d_loc;				// depth of the particle
	double pID;					// the particle identification number
	double pAge;
	double pOSLage;
	double C10Be;
	double C14C;
	double C21Ne;
	double f_remaining;
	double K;
	double powT;
	int age_bin_num;

	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// loop through the bins collecting data
	int n_bins = particle_bins.size();
	//cout << "LINE 1109 n_bins is: " << n_bins << endl;
	int n_parts = 0;
	for (int bn = 0; bn< n_bins; bn++)
	{

		//cout << "LINE 964 starting bin number " << bn << endl;
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			// get the z_location of the particle
			pID = (*part_iter).getType();
			z_loc = (*part_iter).get_zetaLoc();
			s_loc = (*part_iter).getxLoc();
			d_loc = (*part_iter).getdLoc();
			pAge = (*part_iter).getAge();
			pOSLage = (*part_iter).getOSLage();
			C10Be = (*part_iter).getConc_10Be();
			C14C = (*part_iter).getConc_14C();
			C21Ne = (*part_iter).getConc_21Ne();

			//cout << "LINE 1130, pID: " << pID << endl;

			// get the age_bin_number for the particle, based on the age
			age_bin_num = int(double(pAge/age_bin_width));

			//if (age_bin_num < 0)
			//	cout << "LINE 985 uhoh <0" << endl;
			//if (age_bin_num >=n_spacings

			//cout << "Line 988, getting eta" << endl;
			eta_new_local = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + eta[ h_node_us[bn] ];
			//cout << "Line 991, got eta" << endl;

			// only record age stats if the particle is in the soil
			if (z_loc > eta_new_local)
			{
				n_parts++;

				// calcualte the fraction of mass remaining based on chemical weahtering
				// for now we assume no weathering
				powT = pow(pAge,sigma+1);
				K = K_times_D/D;
				f_remaining = exp(-K*powT/(sigma+1));



				total_mass += f_remaining;
				total_mass_weighted_age+=f_remaining*pAge;
				//cout << "LINE 996 in soil, pAge:" << pAge << " pID: " << pID
				//     << " total_mass: " << total_mass << endl;
				//cout << "LINE 1002, recording a particle" << endl;
				if (age_bin_num < n_spacings)
				{
					//cout << "LINE 1005 age_bin_num is: " << age_bin_num
					//     << " and n_spacings is: " << n_spacings << endl;
					age_pdf[age_bin_num]+= f_remaining;

					for (int abn = age_bin_num; abn<n_spacings; abn++)
					{
						age_cdf[abn]+=f_remaining;
					}
				}
				//cout << "LINE 1012, recorded a particle" << endl;
			}
			//else
			//{
				//cout << "LINE 1019 not in soil " << endl;
			//}
			//cout << "LINE 1173 going to next particle" << endl;
			part_iter++;


		}
		//cout << "LINE 1024 ending bin number " << bn << endl;
	}


    //cout << "LINE 1011 looped parts" << endl;

	//cout << "total_mass: " << total_mass << endl;


	for (int abn = 0; abn<n_spacings; abn++)
	{

		age_pdf[abn] = age_pdf[abn]/total_mass;
		age_cdf[abn] = age_cdf[abn]/total_mass;

	}
	mean_age = total_mass_weighted_age/total_mass;
	//cout << "mean_age: " << mean_age << endl;

	//cout << "LINE 1021 looped totals" << endl;

	cdf_out << n_spacings << " " << n_parts << endl;
	cdf_out << t_ime << " " << mean_age << endl;
	pdf_out << n_spacings << " " << n_parts << endl;
	pdf_out << t_ime << " " << mean_age << endl;
	for (int abn = 0; abn<n_spacings; abn++)
	{
		cdf_out << bins[abn] << " " << age_cdf[abn] << endl;
		pdf_out << bins[abn] << " " << age_pdf[abn] << endl;
	}

	//cout << "LINE 1162 looped printing" << endl;

}



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// these function give you the mean age, and cdf and pdf of partcles
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRN_tParticle_bins::print_age_cdfpdf_bins(double t_ime, double max_age, int n_spacings,
						  double K_times_D, double D, double sigma,
						  flowtube& ft,
						  ofstream& cdf_out,
						  ofstream& pdf_out)
{
	// hillslope properties
	vector<double> eta = ft.get_eta();


	double eta_new_local;

	double age_bin_width = max_age/double(n_spacings);
  	vector<double> bins(n_spacings);
  	vector<double> age_pdf;
  	vector<double> age_cdf;
  	vector<double> zero_vec(n_spacings,0.0);



  	for (int i = 0; i<n_spacings; i++)
  	{
  		bins[i] = double(i)*age_bin_width+age_bin_width/2;
	}

	//cout << "LINE 954 looped bins" << endl;

	double s_loc;				// distance downslope of the particle
	double z_loc;				// elevation of the particle
	double d_loc;				// depth of the particle
	double pID;					// the particle identification number
	double pAge;
	double pOSLage;
	double C10Be;
	double C14C;
	double C21Ne;
	double K;
	double powT;
	double f_remaining;

	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// loop through the bins collecting data
	int n_bins = particle_bins.size();
	vector< vector<double> > age_pdf_vec(n_bins);
	vector< vector<double> > age_cdf_vec(n_bins);
	vector<double> mean_age(n_bins,0.0);
  	vector<double> total_mass(n_bins,0.0);
  	vector<double> total_mass_weighted_age(n_bins,0.0);
  	vector<int> n_parts(n_bins,0);
	for (int bn = 0; bn< n_bins; bn++)
	{

		// reset the pdf and cfd vectors
		age_pdf = zero_vec;
		age_cdf = zero_vec;

		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			// get the z_location of the particle
			pID = (*part_iter).getType();
			z_loc = (*part_iter).get_zetaLoc();
			s_loc = (*part_iter).getxLoc();
			d_loc = (*part_iter).getdLoc();
			pAge = (*part_iter).getAge();
			pOSLage = (*part_iter).getOSLage();
			C10Be = (*part_iter).getConc_10Be();
			C14C = (*part_iter).getConc_14C();
			C21Ne = (*part_iter).getConc_21Ne();

			// get the age_bin_number for the particle, based on the age
			int age_bin_num = int(double(pAge/age_bin_width));

			//if (age_bin_num < 0)
			//	cout << "LINE 985 uhoh <0" << endl;
			//if (age_bin_num >=n_spacings

			//cout << "Line 988, getting eta" << endl;
			eta_new_local = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*(s_loc- s_us_h[bn]) + eta[ h_node_us[bn] ];
			//cout << "Line 991, got eta" << endl;

			// only record age stats if the particle is in the soil
			if (z_loc > eta_new_local)
			{
				n_parts[bn]++;
				//cout << "LINE 996 in soil" << endl;
				// calcualte the fraction of mass remaining based on chemical weahtering
				powT = pow(pAge,sigma+1);
				K = K_times_D/D;
				f_remaining = exp(-K*powT/(sigma+1));

				total_mass[bn] += f_remaining;
				total_mass_weighted_age[bn]+=f_remaining*pAge;
				//cout << "LINE 1002, recording a particle" << endl;
				if (age_bin_num < n_spacings)
				{
					//cout << "LINE 1005 age_bin_num is: " << age_bin_num
					//     << " and n_spacings is: " << n_spacings << endl;
					age_pdf[age_bin_num]+= f_remaining;

					for (int abn = age_bin_num; abn<n_spacings; abn++)
					{
						age_cdf[abn]+=f_remaining;
					}
				}
				//cout << "LINE 1012, recorded a particle" << endl;
			}
			//else
			//{
				//cout << "LINE 1019 not in soil " << endl;
			//}

			part_iter++;


		}
		age_pdf_vec[bn] = age_pdf;
		age_cdf_vec[bn] = age_cdf;


	}

    //cout << "LINE 1011 looped parts" << endl;
	for (int bn = 0; bn<n_bins; bn++)
	{
		age_pdf = age_pdf_vec[bn];
		age_cdf = age_cdf_vec[bn];
		for (int abn = 0; abn<n_spacings; abn++)
		{
			age_pdf[abn] = age_pdf[abn]/total_mass[bn];
			age_cdf[abn] = age_cdf[abn]/total_mass[bn];
		}
		mean_age[bn] = total_mass_weighted_age[bn]/total_mass[bn];
		age_pdf_vec[bn] = age_pdf;
		age_cdf_vec[bn] = age_cdf;
	}


	// now loop through bins printing to file
	cdf_out << n_bins << " " << t_ime;
	pdf_out << n_bins << " " << t_ime;
	for (int abn = 0; abn<n_spacings; abn++)
	{
		cdf_out << " -99";
		pdf_out << " -99";
	}
	cdf_out << endl;
	pdf_out << endl;

	pdf_out << "-99 "<< n_spacings;
	cdf_out << "-99 "<< n_spacings;
	for (int abn = 0; abn<n_spacings; abn++)
	{
		cdf_out << " " << bins[abn];
		pdf_out << " " << bins[abn];
	}
	cdf_out << endl;
	pdf_out << endl;

	for (int bn = 0; bn<n_bins; bn++)
	{
		age_cdf = age_cdf_vec[bn];
		age_pdf = age_pdf_vec[bn];


		cdf_out << n_parts[bn] << " " << mean_age[bn];
		pdf_out << n_parts[bn] << " " << mean_age[bn];
		for (int abn = 0; abn<n_spacings; abn++)
		{
			cdf_out << " " << age_cdf[abn];
			pdf_out << " " << age_pdf[abn];
		}
		cdf_out << endl;
		pdf_out << endl;
	}


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// these function give you the mean age, and cdf and pdf of partcles
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRN_tParticle_bins::print_age_cdfpdf_eroded_bulk(double t_ime, double max_age, int n_spacings,
						  vector< list<CRN_tParticle> >& eroded_particle_bins,
						  double K_times_D, double D, double sigma,
						  ofstream& cdf_out,
						  ofstream& pdf_out)
{
	double age_bin_width = max_age/double(n_spacings);
  	vector<double> bins(n_spacings);
  	vector<double> age_pdf(n_spacings,0.0);
  	vector<double> age_cdf(n_spacings,0.0);
  	double total_mass = 0;
  	double total_mass_weighted_age = 0;
  	double mean_age;

  	for (int i = 0; i<n_spacings; i++)
  	{
  		bins[i] = double(i)*age_bin_width+age_bin_width/2;
	}

	//cout << "LINE 954 looped bins" << endl;

	double s_loc;				// distance downslope of the particle
	double z_loc;				// elevation of the particle
	double d_loc;				// depth of the particle
	double pID;					// the particle identification number
	double pAge;
	double pOSLage;
	double C10Be;
	double C14C;
	double C21Ne;
	double K;
	double powT;
	double f_remaining;

	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// loop through the bins collecting data
	int n_bins = eroded_particle_bins.size();
	int n_parts = 0;
	for (int bn = 0; bn< n_bins; bn++)
	{
		//cout << "
		//cout << "LINE 964 starting bin number " << bn << endl;
		// now loop through each particle in the bin
		part_iter = eroded_particle_bins[bn].begin();
		while (part_iter != eroded_particle_bins[bn].end())
		{
			n_parts++;

			// get the z_location of the particle
			pID = (*part_iter).getType();
			z_loc = (*part_iter).get_zetaLoc();
			s_loc = (*part_iter).getxLoc();
			d_loc = (*part_iter).getdLoc();
			pAge = (*part_iter).getAge();
			pOSLage = (*part_iter).getOSLage();
			C10Be = (*part_iter).getConc_10Be();
			C14C = (*part_iter).getConc_14C();
			C21Ne = (*part_iter).getConc_21Ne();

			// get the age_bin_number for the particle, based on the age
			int age_bin_num = int(double(pAge/age_bin_width));

			//cout << "LINE 996 in soil" << endl;
			// calcualte the fraction of mass remaining based on chemical weahtering
			// for now we assume no weathering
			powT = pow(pAge,sigma+1);
			K = K_times_D/D;
			f_remaining = exp(-K*powT/(sigma+1));

			total_mass += f_remaining;
			total_mass_weighted_age+=f_remaining*pAge;
			//cout << "LINE 1002, recording a particle" << endl;
			if (age_bin_num < n_spacings)
			{
				//cout << "LINE 1005 age_bin_num is: " << age_bin_num
				//     << " and n_spacings is: " << n_spacings << endl;
				age_pdf[age_bin_num]+= f_remaining;

				for (int abn = age_bin_num; abn<n_spacings; abn++)
				{
					age_cdf[abn]+=f_remaining;
				}
			}
			part_iter++;
		}
	}

	for (int abn = 0; abn<n_spacings; abn++)
	{
		age_pdf[abn] = age_pdf[abn]/total_mass;
		age_cdf[abn] = age_cdf[abn]/total_mass;
	}
	mean_age = total_mass_weighted_age/total_mass;

	cdf_out << n_spacings << " " << n_parts <<  endl;
	cdf_out << t_ime << " " << mean_age << endl;
	pdf_out << n_spacings << " " << n_parts <<  endl;
	pdf_out << t_ime << " " << mean_age << endl;
	for (int abn = 0; abn<n_spacings; abn++)
	{
		cdf_out << bins[abn] << " " << age_cdf[abn] << endl;
		pdf_out << bins[abn] << " " << age_pdf[abn] << endl;
	}


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// these function give you the mean age, and cdf and pdf of partcles
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRN_tParticle_bins::print_age_cdfpdf_eroded_bins(double t_ime, double max_age, int n_spacings,
								double K_times_D, double D, double sigma,
						  		vector< list<CRN_tParticle> >& particle_bins,
						  		ofstream& cdf_out,
						  		ofstream& pdf_out)
{


	double age_bin_width = max_age/double(n_spacings);
  	vector<double> bins(n_spacings);
  	vector<double> age_pdf;
  	vector<double> age_cdf;
  	vector<double> zero_vec(n_spacings,0.0);



  	for (int i = 0; i<n_spacings; i++)
  	{
  		bins[i] = double(i)*age_bin_width+age_bin_width/2;
	}

	//cout << "LINE 954 looped bins" << endl;

	double s_loc;				// distance downslope of the particle
	double z_loc;				// elevation of the particle
	double d_loc;				// depth of the particle
	double pID;					// the particle identification number
	double pAge;
	double pOSLage;
	double C10Be;
	double C14C;
	double C21Ne;
	double K;
	double powT;
	double f_remaining;

	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// loop through the bins collecting data
	int n_bins = particle_bins.size();
	vector< vector<double> > age_pdf_vec(n_bins);
	vector< vector<double> > age_cdf_vec(n_bins);
	vector<double> mean_age(n_bins,0.0);
  	vector<double> total_mass(n_bins,0.0);
  	vector<double> total_mass_weighted_age(n_bins,0.0);
  	vector<int> n_parts(n_bins,0);
	for (int bn = 0; bn< n_bins; bn++)
	{

		// reset the pdf and cfd vectors
		age_pdf = zero_vec;
		age_cdf = zero_vec;

		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		while (part_iter != particle_bins[bn].end())
		{
			n_parts[bn]++;

			// get the z_location of the particle
			pID = (*part_iter).getType();
			z_loc = (*part_iter).get_zetaLoc();
			s_loc = (*part_iter).getxLoc();
			d_loc = (*part_iter).getdLoc();
			pAge = (*part_iter).getAge();
			pOSLage = (*part_iter).getOSLage();
			C10Be = (*part_iter).getConc_10Be();
			C14C = (*part_iter).getConc_14C();
			C21Ne = (*part_iter).getConc_21Ne();

			// get the age_bin_number for the particle, based on the age
			int age_bin_num = int(double(pAge/age_bin_width));


			//cout << "LINE 996 in soil" << endl;
			// calcualte the fraction of mass remaining based on chemical weahtering
			// for now we assume no weathering
			powT = pow(pAge,sigma+1);
			K = K_times_D/D;
			f_remaining = exp(-K*powT/(sigma+1));

			total_mass[bn] += f_remaining;
			total_mass_weighted_age[bn]+=f_remaining*pAge;
			//cout << "LINE 1002, recording a particle" << endl;
			if (age_bin_num < n_spacings)
			{
				//cout << "LINE 1005 age_bin_num is: " << age_bin_num
				//     << " and n_spacings is: " << n_spacings << endl;
				age_pdf[age_bin_num]+= f_remaining;

				for (int abn = age_bin_num; abn<n_spacings; abn++)
				{
					age_cdf[abn]+=f_remaining;
				}
			}
			//cout << "LINE 1012, recorded a particle" << endl;


			part_iter++;


		}
		age_pdf_vec[bn] = age_pdf;
		age_cdf_vec[bn] = age_cdf;


	}


    //cout << "LINE 1011 looped parts" << endl;
	for (int bn = 0; bn<n_bins; bn++)
	{
		age_pdf = age_pdf_vec[bn];
		age_cdf = age_cdf_vec[bn];
		if (total_mass[bn] != 0)
		{
			for (int abn = 0; abn<n_spacings; abn++)
			{

				age_pdf[abn] = age_pdf[abn]/total_mass[bn];
				age_cdf[abn] = age_cdf[abn]/total_mass[bn];
			}
			mean_age[bn] = total_mass_weighted_age[bn]/total_mass[bn];
		}
		else
		{
			mean_age[bn] = 0.0;
			age_pdf = zero_vec;
			age_cdf = zero_vec;
		}
		age_pdf_vec[bn] = age_pdf;
		age_cdf_vec[bn] = age_cdf;
	}


	// now loop through bins printing to file
	cdf_out << n_bins << " " << t_ime;
	pdf_out << n_bins << " " << t_ime;
	for (int abn = 0; abn<n_spacings; abn++)
	{
		cdf_out << " -99";
		pdf_out << " -99";
	}
	cdf_out << endl;
	pdf_out << endl;

	pdf_out << "-99 "<< n_spacings;
	cdf_out << "-99 "<< n_spacings;
	for (int abn = 0; abn<n_spacings; abn++)
	{
		cdf_out << " " << bins[abn];
		pdf_out << " " << bins[abn];
	}
	cdf_out << endl;
	pdf_out << endl;

	for (int bn = 0; bn<n_bins; bn++)
	{
		age_cdf = age_cdf_vec[bn];
		age_pdf = age_pdf_vec[bn];


		cdf_out << n_parts[bn] << " " << mean_age[bn];
		pdf_out << n_parts[bn] << " " << mean_age[bn];
		for (int abn = 0; abn<n_spacings; abn++)
		{
			cdf_out << " " << age_cdf[abn];
			pdf_out << " " << age_pdf[abn];
		}
		cdf_out << endl;
		pdf_out << endl;
	}


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function caluclates values within cells
//
// NUMBERING OF THE VERTICES AND CELLS
//
// the total number of cells is (n_depthintervals_soil+n_depthintervals_parent)*n_bins
// on each edge of a bin there are (n_depthintervals_soil+n_depthintervals_parent)+1 depth
// intervals
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRN_tParticle_bins::cell_printing_vtk(double t_ime, flowtube ft,
								string vtk_fname,
								int n_depthintervals_soil, int n_depthintervals_parent,
								double bottom_depth)
{

	int total_depth_intervals = n_depthintervals_soil+n_depthintervals_parent;

	// get the file printed
	string time_bit = itoa( int(t_ime+0.5) );
	string vtk_ext = ".vtk";
	string fname = vtk_fname+time_bit+vtk_ext;
	//cout << "LINE 1891 CRN_tparticle_ bins, vtk filename is: " << fname << " and time: " << t_ime << endl;
	ofstream vtk_out;
	vtk_out.open(fname.c_str());

	// vectors for holding the vertices
	vector<double> verts_s;
	vector<double> verts_z;

	// this loops through each bin, getting the
	double us_zeta, us_eta;

	// get zeta and eta
	vector<double> eta = ft.get_eta();
	vector<double> zeta = ft.get_zeta();
	vector<double> s = ft.get_s_h();

	int n_z_nodes = eta.size();
	int n_bin_edges = bin_edge_loc.size();

	for (int i = 0; i<n_z_nodes; i++)
	{
		cout << "s["<< i<<"]: " << s[i] << " zeta: " << zeta[i] << " eta: " << eta[i] << endl;
	}
	for (int i = 0; i<n_bin_edges; i++)
	{
		cout << "bin_edge_loc["<< i<<"]: " << bin_edge_loc[i] << endl;
	}

	// loop through the bins collecting data
	int n_bins = particle_bins.size();

	for (int i = 0; i<n_z_nodes; i++)
	{
		cout << "s["<< i<<"]: " << s[i] << " zeta: " << zeta[i] << " eta: " << eta[i] << endl;
	}
	for (int i = 0; i<n_bin_edges; i++)
	{
		cout << "bin_edge_loc["<< i<<"]: " << bin_edge_loc[i] << endl;
	}
	cout << "n_bins: " << n_bins << endl;

	// the vectors hold the indices to the vertices
	vector<int> cell_node1;
	vector<int> cell_node2;
	vector<int> cell_node3;
	vector<int> cell_node4;
	vector<int> cell_code;			// used to differentiate soil from parent material. The numbers are
									// arbitrary, but for now 0 == parent and 1 == soil
	int soil_code = 1;
	int parent_code = 0;

	double us_soil_depthinterval,
	       us_parent_depthinterval;

	double s_up;
	double s_h_up;

	double fuzzy_boundary = 0.001;
	// now loop through the bins, calcualting vertices
	for (int bn = 0; bn< n_bins; bn++)
	{
		s_h_up = s[ h_node_us[bn]];
		s_up = bin_edge_loc[bn];

		// if the upslope h location and the upslope bin locations are the same...
		if (s_h_up <= s_up+fuzzy_boundary &&  s_h_up >= s_up-fuzzy_boundary )
		{
			us_zeta = zeta[ h_node_us[bn] ];
			us_eta = eta[ h_node_us[bn] ];
		}
		// if not ...
		else
		{
			us_zeta = ( (zeta[ h_node_us[bn] ]-zeta[ h_node_ds[bn] ])/dx_h[bn] )
			           *(s_up-s_h_up)+ zeta[ h_node_ds[bn] ];
			us_eta = ( (eta[ h_node_us[bn] ]-eta[ h_node_ds[bn] ])/dx_h[bn] )
			           *(s_up-s_h_up)+ eta[ h_node_ds[bn] ];
		}

		if (bn == 0)
		{
			us_zeta = zeta[ h_node_us[bn] ];
			us_eta = eta[ h_node_us[bn] ];
		}

		// calcualte the depth of the bins in the soil
		us_soil_depthinterval = (us_zeta-us_eta)/double(n_depthintervals_soil);
		us_parent_depthinterval = (bottom_depth-(us_zeta-us_eta))/double(n_depthintervals_parent);

		// first get the soil nodes
		for (int sdi = 0; sdi<n_depthintervals_soil; sdi++)
		{
			// get the s-coordinate vertex
			verts_s.push_back(bin_edge_loc[bn]);

			// get the vertical index
			verts_z.push_back( us_zeta - (double(sdi))*us_soil_depthinterval);
		}

		// now the parent material nodes
		for (int pdi = 0; pdi<n_depthintervals_parent; pdi++)
		{
			// get the s-coordinate vertex
			verts_s.push_back(bin_edge_loc[bn]);

			// get the vertical index
			verts_z.push_back( us_eta- (double(pdi))*us_parent_depthinterval);
		}

		// and now the bottom node
		verts_s.push_back(bin_edge_loc[bn]);
		verts_z.push_back(us_zeta-bottom_depth);

	}

	// loop through cells plotting vertices
	//for (int bn = 0; bn<n_bins-1; bn++)
	for (int bn = 0; bn<n_bins-1; bn++)
	{
		for (int di = 0; di<total_depth_intervals; di++)
		{
			cell_node1.push_back( bn*(total_depth_intervals+1)+di );
			cell_node2.push_back( (bn+1)*(total_depth_intervals+1)+di );
			cell_node3.push_back( (bn+1)*(total_depth_intervals+1)+di+1 );;
			cell_node4.push_back( bn*(total_depth_intervals+1)+di+1 );

			if (di < n_depthintervals_soil)
			{
				cell_code.push_back(soil_code);
			}
			else
			{
				cell_code.push_back(parent_code);
			}
		}
	}

	// now print the vtk file
	// find the number of particles
	int n_verts = verts_z.size();
	//cout << "verts size: " << n_verts << endl;

	vtk_out << "# vtk DataFile Version 2.0" << endl << "Unstructured Grid Ptrack_cells"
	        << endl << "ASCII" << endl << endl << "DATASET UNSTRUCTURED_GRID" << endl
	        << "POINTS " << n_verts << " float" << endl;

	for (int i = 0; i< n_verts; i++)
	{
		vtk_out << verts_s[i] << " " << verts_z[i] << " 0.0" <<endl;
	}

	int n_cells = cell_code.size();
	vtk_out << endl << "CELLS "<<n_cells << " " << n_cells*5 << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_out << 4 << " " << cell_node1[i] << " " << cell_node2[i] << " "
		        << cell_node3[i] << " " << cell_node4[i] << endl;
	}

	vtk_out << endl << "CELL_TYPES " << n_cells << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_out << "9" << endl;
	}


	vtk_out << endl << "CELL_DATA " << n_cells << endl;
	vtk_out << "SCALARS SOIL_OR_PARENT int 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_out << cell_code[i] <<endl;
	}

	vtk_out << "SCALARS CELL_NUM int 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_out << i <<endl;
	}

	vtk_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function caluclates values within cells
//
// NUMBERING OF THE VERTICES AND CELLS
//
// the total number of cells is (n_depthintervals_soil+n_depthintervals_parent)*n_bins
// on each edge of a bin there are (n_depthintervals_soil+n_depthintervals_parent)+1 depth
// intervals
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRN_tParticle_bins::cell_and_particle_printing_vtk(double t_ime, flowtube ft,
								 string vtk_particle_fname, string vtk_cell_fname,
								int n_depthintervals_soil, int n_depthintervals_parent,
								double bottom_depth)
{

	int total_depth_intervals = n_depthintervals_soil+n_depthintervals_parent;


	// get the file printed
	string time_bit = itoa( int(t_ime+0.5) );
	string vtk_cell_ext = ".vtk";
	string fname = vtk_cell_fname+time_bit+vtk_cell_ext;
	ofstream vtk_cell_out;
	vtk_cell_out.open(fname.c_str());

	// vectors for holding the vertices
	vector<double> verts_s;
	vector<double> verts_z;

	// this loops through each bin, getting the
	double us_zeta, us_eta;

	// get zeta, eta and s
	vector<double> eta = ft.get_eta();
	vector<double> zeta = ft.get_zeta();
	vector<double> s = ft.get_s_h();

	// loop through the bins collecting data
	int n_bins = particle_bins.size();

	// the vectors hold the indices to the vertices
	vector<int> cell_node1;
	vector<int> cell_node2;
	vector<int> cell_node3;
	vector<int> cell_node4;
	vector<int> cell_code;			// used to differentiate soil from parent material. The numbers are
									// arbitrary, but for now 0 == parent and 1 == soil
	int soil_code = 1;
	int parent_code = 0;

	double us_soil_depthinterval,
	       us_parent_depthinterval;

	double s_up;
	double s_h_up;

	double fuzzy_boundary = 0.001;
	// now loop through the bins, calcualting vertices
	for (int bn = 0; bn< n_bins; bn++)
	{
		s_h_up = s[ h_node_us[bn]];
		s_up = bin_edge_loc[bn];

		// if the upslope h location and the upslope bin locations are the same...
		if (s_h_up <= s_up+fuzzy_boundary &&  s_h_up >= s_up-fuzzy_boundary )
		{
			us_zeta = zeta[ h_node_us[bn] ];
			us_eta = eta[ h_node_us[bn] ];
		}
		// if not ...
		else
		{
			us_zeta = ( (zeta[ h_node_us[bn] ]-zeta[ h_node_ds[bn] ])/dx_h[bn] )
			           *(s_up-s_h_up)+ zeta[ h_node_ds[bn] ];
			us_eta = ( (eta[ h_node_us[bn] ]-eta[ h_node_ds[bn] ])/dx_h[bn] )
			           *(s_up-s_h_up)+ eta[ h_node_ds[bn] ];
		}

		if (bn == 0)
		{
			us_zeta = zeta[ h_node_us[bn] ];
			us_eta = eta[ h_node_us[bn] ];
		}

		// calcualte the depth of the bins in the soil
		us_soil_depthinterval = (us_zeta-us_eta)/double(n_depthintervals_soil);
		us_parent_depthinterval = (bottom_depth-(us_zeta-us_eta))/double(n_depthintervals_parent);

		// first get the soil nodes
		for (int sdi = 0; sdi<n_depthintervals_soil; sdi++)
		{
			// get the s-coordinate vertex
			verts_s.push_back(bin_edge_loc[bn]);

			// get the vertical index
			verts_z.push_back( us_zeta - (double(sdi))*us_soil_depthinterval);
		}

		// now the parent material nodes
		for (int pdi = 0; pdi<n_depthintervals_parent; pdi++)
		{
			// get the s-coordinate vertex
			verts_s.push_back(bin_edge_loc[bn]);

			// get the vertical index
			verts_z.push_back( us_eta- (double(pdi))*us_parent_depthinterval);
		}

		// and now the bottom node
		verts_s.push_back(bin_edge_loc[bn]);
		verts_z.push_back(us_zeta-bottom_depth);

	}

	// loop through cells plotting vertices
	//for (int bn = 0; bn<n_bins-1; bn++)
	for (int bn = 0; bn<n_bins-1; bn++)
	{
		for (int di = 0; di<total_depth_intervals; di++)
		{
			cell_node1.push_back( bn*(total_depth_intervals+1)+di );
			cell_node2.push_back( (bn+1)*(total_depth_intervals+1)+di );
			cell_node3.push_back( (bn+1)*(total_depth_intervals+1)+di+1 );;
			cell_node4.push_back( bn*(total_depth_intervals+1)+di+1 );

			if (di < n_depthintervals_soil)
			{
				cell_code.push_back(soil_code);
			}
			else
			{
				cell_code.push_back(parent_code);
			}
		}
	}

	// now print the vtk file
	// find the number of particles
	int n_verts = verts_z.size();
	//cout << "verts size: " << n_verts << endl;

	vtk_cell_out << "# vtk DataFile Version 2.0" << endl << "Unstructured Grid Ptrack_cells"
	        << endl << "ASCII" << endl << endl << "DATASET UNSTRUCTURED_GRID" << endl
	        << "POINTS " << n_verts << " float" << endl;

	for (int i = 0; i< n_verts; i++)
	{
		vtk_cell_out << verts_s[i] << " " << verts_z[i] << " 0.0" <<endl;
	}

	int n_cells = cell_code.size();
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
	vtk_cell_out << "SCALARS SOIL_OR_PARENT int 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << cell_code[i] <<endl;
	}

	vtk_cell_out << "SCALARS CELL_NUM int 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << i <<endl;
	}

	// now to do the particles
	vector<double> s_loc;				// distance downslope of the particle
	vector<double> z_loc;				// elevation of the particle
	vector<double> d_loc;				// depth of the particle
	vector<int> pID;					// the particle identification number
	vector<double> pAge;
	vector<double> pOSLage;
	vector<double> C10Be;
	vector<double> C14C;
	vector<double> C21Ne;
	vector<int> soil_switch;

	list<CRN_tParticle>::iterator part_iter;	// list iterator
	string vtk_particle_ext = ".vtk";
	fname = vtk_particle_fname+time_bit+vtk_particle_ext;
	//cout << "LINE 1149 CRN_tparticle_ bins, vtk filename is: " << fname << " and time: " << t_ime << endl;
	ofstream vtk_particle_out;
	vtk_particle_out.open(fname.c_str());

	// hillslope properties

	double eta_new_local;
	// loop through all the bins
	for (int bn = 0; bn< n_bins-1; bn++)
	{
		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		int counter = 0;
		while (part_iter != particle_bins[bn].end())
		{
			// get the data from each particle
			s_loc.push_back( (*part_iter).getxLoc() );
			z_loc.push_back( (*part_iter).get_zetaLoc() );
			d_loc.push_back( (*part_iter).getdLoc() );
			pID.push_back( (*part_iter).getType() );
			pAge.push_back( (*part_iter).getAge() );
			pOSLage.push_back( (*part_iter).getOSLage() );
			C10Be.push_back( (*part_iter).getConc_10Be() );
			C14C.push_back( (*part_iter).getConc_14C() );
			C21Ne.push_back( (*part_iter).getConc_21Ne() );

			eta_new_local = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*( (*part_iter).getxLoc()- s_us_h[bn]) + eta[ h_node_us[bn] ];
			if ( (*part_iter).get_zetaLoc() > eta_new_local)
			{
				soil_switch.push_back(1);
			}
			else
			{
				soil_switch.push_back(0);
			}

			part_iter++;
			counter ++;
		}
		//cout << "bin number " << bn << " and counter is: " << counter << endl;
		//cout << "n_parts_total: " << d_loc.size() << endl;
	}

	// find the number of particles
	int n_parts = d_loc.size();
	//cout << "d_loc size: " << n_parts << endl;

	// now loop back through the bins, checking which particles are in which depth interval
	vector<int> particle_in_depth_interval(n_parts,-1);	// default is -1 for particles that don't get
														// picked up by a bin
	int n_parts_in_bin;
	int depth_interval_index;

	int p_start = 0;			// index of the starting particle in bin
	int p_end = 0;			// index of the finishing particle in bin
	int pip_test;			// this is 1 if particle is in cell and 0 otherwise

	vector<double> sverts_for_pip(4);
	vector<double> zverts_for_pip(4);
	for (int bn = 0; bn< n_bins-1; bn++)
	{
		n_parts_in_bin = particle_bins[bn].size();
		p_start = p_end;
		p_end = p_start+n_parts_in_bin;

		for (int di = 0; di< total_depth_intervals; di++)
		{
			// get the index for the cell
			depth_interval_index = bn*total_depth_intervals+di;

			// populate the vertices vectors
			sverts_for_pip[0] = verts_s[ cell_node1[depth_interval_index] ];
			sverts_for_pip[1] = verts_s[ cell_node2[depth_interval_index] ];
			sverts_for_pip[2] = verts_s[ cell_node3[depth_interval_index] ];
			sverts_for_pip[3] = verts_s[ cell_node4[depth_interval_index] ];

			zverts_for_pip[0] = verts_z[ cell_node1[depth_interval_index] ];
			zverts_for_pip[1] = verts_z[ cell_node2[depth_interval_index] ];
			zverts_for_pip[2] = verts_z[ cell_node3[depth_interval_index] ];
			zverts_for_pip[3] = verts_z[ cell_node4[depth_interval_index] ];

			// now loop through all the particles in this bin
			for (int part_i = p_start; part_i < p_end; part_i++)
			{
				pip_test = pnpoly(sverts_for_pip, zverts_for_pip, s_loc[part_i], z_loc[part_i]);
				if (pip_test == 1)
				{
					particle_in_depth_interval[part_i] = depth_interval_index;
				}
			}
		}
	}

	vtk_particle_out << "# vtk DataFile Version 2.0" << endl << "Unstructured Grid Ptrack"
	        << endl << "ASCII" << endl << endl << "DATASET UNSTRUCTURED_GRID" << endl
	        << "POINTS " << n_parts << " float" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << s_loc[i] << " " << z_loc[i] << " 0.0" <<endl;
	}

	vtk_particle_out << endl << "POINT_DATA "<<n_parts << endl 
            << "SCALARS Particle_age float 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << pAge[i] <<endl;
	}

	vtk_particle_out << "SCALARS Be10_conc float 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << C10Be[i] <<endl;
	}

	vtk_particle_out << "SCALARS Ptype int 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
			vtk_particle_out << pID[i] <<endl;
	}

	vtk_particle_out << "SCALARS depth_interval_number int 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << particle_in_depth_interval[i] <<endl;
	}

	// close the files
	vtk_cell_out.close();
	vtk_particle_out.close();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// 
// This function prints a simple version of the particles that just gives location, 
// type and the mass remaining
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRN_tParticle_bins::vtk_print_basic_volume_particles(double t_ime, 
								 string vtk_particle_fname, int reference_frame_switch)
{

 	// get the file printed
	string time_bit = itoa( int(t_ime+0.5) );
	string vtk_cell_ext = ".vtk";

	vector<double> s_loc;				// distance downslope of the particle
	vector<double> z_loc;				// elevation of the particle
	vector<double> d_loc;				// depth of the particle
	vector<int> pType;					// the particle identification number
	vector<double> Mass_remain;
	double Mass,StartingMass;

  // get the iterator and print the filename
	list<CRN_tParticle>::iterator part_iter;	// list iterator
	string vtk_particle_ext = ".vtk";
	string fname = vtk_particle_fname+time_bit+vtk_particle_ext;
	cout << "LINE 3177 CRN_tparticle_ bins, vtk filename is: " << fname << " and time: " << t_ime << endl;
	ofstream vtk_particle_out;
	vtk_particle_out.open(fname.c_str());

	// loop through all the bins
	for (int bn = 0; bn< n_bins; bn++)
	{
		//cout << "CRN_tparticle.cpp, LINE 2671 bin number is: " << bn << endl;

		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		int counter = 0;
		while (part_iter != particle_bins[bn].end())
		{
			// get the data from each particle
			s_loc.push_back( (*part_iter).getxLoc() );
			z_loc.push_back( (*part_iter).get_zetaLoc() );
			d_loc.push_back( (*part_iter).getdLoc() );
      pType.push_back( (*part_iter).getType() );
      Mass = (*part_iter).getMass();
      StartingMass = (*part_iter).getStartingMass();
      Mass_remain.push_back(Mass/StartingMass);

			part_iter++;
		}
	}

	// find the number of particles
	int n_parts = d_loc.size();
	//cout << "d_loc size: " << n_parts << endl;

	vtk_particle_out << "# vtk DataFile Version 2.0" << endl << "Unstructured Grid Ptrack"
	        << endl << "ASCII" << endl << endl << "DATASET UNSTRUCTURED_GRID" << endl
	        << "POINTS " << n_parts << " float" << endl;

	if (reference_frame_switch == 1)
	{
		for (int i = 0; i< n_parts; i++)
		{
			vtk_particle_out << s_loc[i] << " " << -d_loc[i] << " 0.0" <<endl;
		}
	}
	else
	{
		for (int i = 0; i< n_parts; i++)
		{
			vtk_particle_out << s_loc[i] << " " << z_loc[i] << " 0.0" <<endl;
		}
	}

	vtk_particle_out << endl << "POINT_DATA "<<n_parts << endl 
                   << "SCALARS Particle_type int 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << pType[i] <<endl;
	}

	vtk_particle_out << "SCALARS Fraction_of_mass_remaining float 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
			vtk_particle_out << Mass_remain[i] <<endl;
	}


	// close the files
	vtk_particle_out.close();

}
								
								


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function caluclates values within cells
//
// NUMBERING OF THE VERTICES AND CELLS
//
// the total number of cells is (n_depthintervals_soil+n_depthintervals_parent)*n_bins
// on each edge of a bin there are (n_depthintervals_soil+n_depthintervals_parent)+1 depth
// intervals
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRN_tParticle_bins::cell_and_particle_chemistry_printing_vtk(double t_ime, flowtube ft,
								 Particle_info& pi,
								 string vtk_particle_fname, string vtk_cell_fname,
								int n_depthintervals_soil, int n_depthintervals_parent,
								double bottom_depth, int reference_frame_switch)
{

	int total_depth_intervals = n_depthintervals_soil+n_depthintervals_parent;


	// get the file printed
	string time_bit = itoa( int(t_ime+0.5) );
	string vtk_cell_ext = ".vtk";
	string fname = vtk_cell_fname+time_bit+vtk_cell_ext;
	//cout << "LINE 1891 CRN_tparticle_ bins, vtk filename is: " << fname << " and time: " << t_ime << endl;
	ofstream vtk_cell_out;
	vtk_cell_out.open(fname.c_str());

	// vectors for holding the vertices
	vector<double> verts_s;
	vector<double> verts_z;
	vector<double> verts_d;

	// this loops through each bin, getting the
	double us_zeta, us_eta;

	// get zeta, eta and s
	vector<double> eta = ft.get_eta();
	vector<double> zeta = ft.get_zeta();
	vector<double> s = ft.get_s_h();

	// loop through the bins collecting data
	int n_bins = particle_bins.size();

	// the vectors hold the indices to the vertices
	vector<int> cell_node1;
	vector<int> cell_node2;
	vector<int> cell_node3;
	vector<int> cell_node4;
	vector<int> cell_code;			// used to differentiate soil from parent material. The numbers are
									// arbitrary, but for now 0 == parent and 1 == soil
	int soil_code = 1;
	int parent_code = 0;

	double us_soil_depthinterval,
	       us_parent_depthinterval;

	double s_up;
	double s_h_up;

	double fuzzy_boundary = 0.001;
	// now loop through the bins, calcualting vertices
	for (int bn = 0; bn< n_bins; bn++)
	{
		s_h_up = s[ h_node_us[bn]];
		s_up = bin_edge_loc[bn];

		// if the upslope h location and the upslope bin locations are the same...
		if (s_h_up <= s_up+fuzzy_boundary &&  s_h_up >= s_up-fuzzy_boundary )
		{
			us_zeta = zeta[ h_node_us[bn] ];
			us_eta = eta[ h_node_us[bn] ];
		}
		// if not ...
		else
		{
			us_zeta = ( (zeta[ h_node_us[bn] ]-zeta[ h_node_ds[bn] ])/dx_h[bn] )
			           *(s_up-s_h_up)+ zeta[ h_node_ds[bn] ];
			us_eta = ( (eta[ h_node_us[bn] ]-eta[ h_node_ds[bn] ])/dx_h[bn] )
			           *(s_up-s_h_up)+ eta[ h_node_ds[bn] ];
		}

		if (bn == 0)
		{
			us_zeta = zeta[ h_node_us[bn] ];
			us_eta = eta[ h_node_us[bn] ];
		}

		// calcualte the depth of the bins in the soil
		us_soil_depthinterval = (us_zeta-us_eta)/double(n_depthintervals_soil);
		us_parent_depthinterval = (bottom_depth-(us_zeta-us_eta))/
                            double(n_depthintervals_parent);

		// first get the soil nodes
		for (int sdi = 0; sdi<n_depthintervals_soil; sdi++)
		{
			// get the s-coordinate vertex
			verts_s.push_back(bin_edge_loc[bn]);

			// get the vertical index
			verts_z.push_back( us_zeta - (double(sdi))*us_soil_depthinterval);
			verts_d.push_back( double(sdi)*us_soil_depthinterval );
		}

		// now the parent material nodes
		for (int pdi = 0; pdi<n_depthintervals_parent; pdi++)
		{
			// get the s-coordinate vertex
			verts_s.push_back(bin_edge_loc[bn]);

			// get the vertical index
			verts_z.push_back( us_eta- (double(pdi))*us_parent_depthinterval);
			verts_d.push_back( double(pdi)*us_parent_depthinterval);
		}

		// and now the bottom node
		verts_s.push_back(bin_edge_loc[bn]);
		verts_z.push_back(us_zeta-bottom_depth);
		verts_d.push_back(bottom_depth);

	}


	// loop through cells plotting vertices
	//for (int bn = 0; bn<n_bins-1; bn++)
	for (int bn = 0; bn<n_bins-1; bn++)
	{
		for (int di = 0; di<total_depth_intervals; di++)
		{
			cell_node1.push_back( bn*(total_depth_intervals+1)+di );
			cell_node2.push_back( (bn+1)*(total_depth_intervals+1)+di );
			cell_node3.push_back( (bn+1)*(total_depth_intervals+1)+di+1 );;
			cell_node4.push_back( bn*(total_depth_intervals+1)+di+1 );

			if (di < n_depthintervals_soil)
			{
				cell_code.push_back(soil_code);
			}
			else
			{
				cell_code.push_back(parent_code);
			}
		}
	}

	// now print the vtk file
	// find the number of particles
	int n_verts = verts_z.size();
	//cout << "verts size: " << n_verts << endl;

	vtk_cell_out << "# vtk DataFile Version 2.0" << endl << "Unstructured Grid Ptrack_cells"
	        << endl << "ASCII" << endl << endl << "DATASET UNSTRUCTURED_GRID" << endl
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

	int n_cells = cell_code.size();
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
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << cell_code[i] <<endl;
	}

	vtk_cell_out << "SCALARS CELL_NUM int 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << i <<endl;
	}

	// now to do the particles

	vector<double> s_loc;				// distance downslope of the particle
	vector<double> z_loc;				// elevation of the particle
	vector<double> d_loc;				// depth of the particle
	vector<int> pID;					// the particle identification number
	vector<double> pAge;
	vector<double> pOSLage;
	vector<double> C10Be;
	vector<double> C14C;
	vector<double> C21Ne;
	vector<int> soil_switch;
	vector<double> mass_fraction;
	vector<double> rate_fraction;
	vector<double> clay_fraction;

	list<CRN_tParticle>::iterator part_iter;	// list iterator
	string vtk_particle_ext = ".vtk";
	fname = vtk_particle_fname+time_bit+vtk_particle_ext;
	//cout << "LINE 1149 CRN_tparticle_ bins, vtk filename is: " << fname << " and time: " << t_ime << endl;
	ofstream vtk_particle_out;
	vtk_particle_out.open(fname.c_str());

	// hillslope properties

	double eta_new_local;
	// loop through all the bins
	for (int bn = 0; bn< n_bins-1; bn++)
	{
		//cout << "CRN_tparticle.cpp, LINE 2671 bin number is: " << bn << endl;

		// now loop through each particle in the bin
		part_iter = particle_bins[bn].begin();
		int counter = 0;
		while (part_iter != particle_bins[bn].end())
		{


			// get the data from each particle
			s_loc.push_back( (*part_iter).getxLoc() );
			z_loc.push_back( (*part_iter).get_zetaLoc() );
			d_loc.push_back( (*part_iter).getdLoc() );
			pID.push_back( (*part_iter).getType() );
			pAge.push_back( (*part_iter).getAge() );
			pOSLage.push_back( (*part_iter).getOSLage() );
			C10Be.push_back( (*part_iter).getConc_10Be() );
			C14C.push_back( (*part_iter).getConc_14C() );
			C21Ne.push_back( (*part_iter).getConc_21Ne() );

			vector<double> part_characteristics = pi.calculate_weathering_on_demand(
				                                  (*part_iter).getAge() , (*part_iter).getType() );
			mass_fraction.push_back( part_characteristics[0] );
			clay_fraction.push_back( part_characteristics[1] );
			rate_fraction.push_back( part_characteristics[2] );

			eta_new_local = ((eta[ h_node_ds[bn] ] - eta[ h_node_us[bn] ])/
								dx_h[bn])*( (*part_iter).getxLoc()- s_us_h[bn]) + eta[ h_node_us[bn] ];
			if ( (*part_iter).get_zetaLoc() > eta_new_local)
			{
				soil_switch.push_back(1);
			}
			else
			{
				soil_switch.push_back(0);
			}

			part_iter++;
			counter ++;
		}
		//cout << "bin number " << bn << " and counter is: " << counter << endl;
		//cout << "n_parts_total: " << d_loc.size() << endl;
	}

	// find the number of particles
	int n_parts = d_loc.size();
	//cout << "d_loc size: " << n_parts << endl;

	// now loop back through the bons, checking which particles are in which depth interval
	vector<int> particle_in_depth_interval(n_parts,total_depth_intervals-1);	// default is -1 for particles that don't get
														// picked up by a bin
	int n_parts_in_bin;
	int depth_interval_index;

	int p_start = 0;			// index of the starting particle in bin
	int p_end = 0;			// index of the finishing particle in bin
	int pip_test;			// this is 1 if particle is in cell and 0 otherwise

	vector<double> sverts_for_pip(4);
	vector<double> zverts_for_pip(4);
	for (int bn = 0; bn< n_bins-1; bn++)
	{
		n_parts_in_bin = particle_bins[bn].size();
		p_start = p_end;
		p_end = p_start+n_parts_in_bin;

		for (int di = 0; di< total_depth_intervals; di++)
		{
			//cout << "depth_interval: " << di << " of " << total_depth_intervals << endl;

			// get the index for the cell
			depth_interval_index = bn*total_depth_intervals+di;

			// populate the vertices vectors
			sverts_for_pip[0] = verts_s[ cell_node1[depth_interval_index] ];
			sverts_for_pip[1] = verts_s[ cell_node2[depth_interval_index] ];
			sverts_for_pip[2] = verts_s[ cell_node3[depth_interval_index] ];
			sverts_for_pip[3] = verts_s[ cell_node4[depth_interval_index] ];

			zverts_for_pip[0] = verts_z[ cell_node1[depth_interval_index] ];
			zverts_for_pip[1] = verts_z[ cell_node2[depth_interval_index] ];
			zverts_for_pip[2] = verts_z[ cell_node3[depth_interval_index] ];
			zverts_for_pip[3] = verts_z[ cell_node4[depth_interval_index] ];

			// now loop through all the particles in this bin
			//cout << "p_start: " << p_start << " and p_end: " << p_end << endl;
			//cout << "number of parts: " << n_parts << endl;
			for (int part_i = p_start; part_i < p_end; part_i++)
			{
				pip_test = pnpoly(sverts_for_pip, zverts_for_pip, s_loc[part_i], z_loc[part_i]);
				if (pip_test == 1)
				{
					particle_in_depth_interval[part_i] = depth_interval_index;
				}
			}
		}
	}
	//cout << " done with getting verts " << endl;
	//cout << "ncells is: " << n_cells << endl;

	// now create particle based cell vectors
	vector<int> n_parts_in_cell(n_cells,0);
	vector< vector<double> > age_in_depth_intervals(n_cells);
	vector<double> age_in_depth_intervals_sum(n_cells,0.0);
	vector<double> mean_age_in_depth_intervals(n_cells,-99.0);

	vector< vector<double> > mfrac_in_depth_intervals(n_cells);
	vector<double> mfrac_in_depth_intervals_sum(n_cells,0.0);
	vector<double> mean_mfrac_in_depth_intervals(n_cells,-99.0);

	vector< vector<double> > cfrac_in_depth_intervals(n_cells);
	vector<double> cfrac_in_depth_intervals_sum(n_cells,0.0);
	vector<double> mean_cfrac_in_depth_intervals(n_cells,-99.0);

	vector< vector<double> > rfrac_in_depth_intervals(n_cells);
	vector<double> rfrac_in_depth_intervals_sum(n_cells,0.0);
	vector<double> mean_rfrac_in_depth_intervals(n_cells,-99.0);

	vector< vector<double> > C10Be_in_depth_intervals(n_cells);
	vector<double> C10Be_in_depth_intervals_sum(n_cells,0.0);
	vector<double> mean_C10Be_in_depth_intervals(n_cells,-99.0);

	// loop through the particles collecting information from each particle and dropping it into
	// cells
	//cout << " looping through parts " << endl;
	//cout << "n_parts: " << n_parts << endl;
	//cout << "pindi size: " << particle_in_depth_interval.size() << endl;
	//cout << "page size: " << pAge.size() << endl;
	//cout << "mf size: " << mass_fraction.size() << endl;
	//cout << "cf size: " << clay_fraction.size() << endl;
	//cout << "rf size: " << rate_fraction.size() << endl;
	//cout << "n_cells is: " << n_cells << endl;
	for (int part_i = 0; part_i<n_parts; part_i++)
	{
		//cout << "i: " << part_i << " depth interval: " << particle_in_depth_interval[part_i] << endl;
		age_in_depth_intervals[ particle_in_depth_interval[part_i] ].push_back(pAge[part_i]);
		mfrac_in_depth_intervals[ particle_in_depth_interval[part_i] ].push_back(mass_fraction[part_i]);
		cfrac_in_depth_intervals[ particle_in_depth_interval[part_i] ].push_back(clay_fraction[part_i]);
		rfrac_in_depth_intervals[ particle_in_depth_interval[part_i] ].push_back(rate_fraction[part_i]);
		C10Be_in_depth_intervals[ particle_in_depth_interval[part_i] ].push_back(C10Be[part_i]);
		n_parts_in_cell[ particle_in_depth_interval[part_i] ]++;
	}
	//cout << " done with getting n_parts in cell" << endl;

	//for (int i = 0; i< n_cells; i++)
	//{
	//	cout << "cell number: " << i << " and n_parts in cell: " 
  //         << n_parts_in_cell[i] << endl;
	//}

	vector<double>::iterator v_iter;
	// now get the mean age
	for (int cell_i=0; cell_i < n_cells; cell_i++)
	{
		v_iter = age_in_depth_intervals[ cell_i ].begin();
		while (v_iter != age_in_depth_intervals[ cell_i ].end())
		{
			age_in_depth_intervals_sum[cell_i] += (*v_iter);
			//cout << " age is: " << (*v_iter) << endl;

			v_iter++;
		}

		if (n_parts_in_cell[cell_i] >0)
		{
			mean_age_in_depth_intervals[cell_i] = age_in_depth_intervals_sum[cell_i] /
											double( n_parts_in_cell[cell_i] );
		}
		else
		{
			mean_age_in_depth_intervals[cell_i] = -100;
		}
	}

	// now get the mean mass frac
	for (int cell_i=0; cell_i < n_cells; cell_i++)
	{
		v_iter = mfrac_in_depth_intervals[ cell_i ].begin();
		while (v_iter != mfrac_in_depth_intervals[ cell_i ].end())
		{
			mfrac_in_depth_intervals_sum[cell_i] += (*v_iter);
			v_iter++;
		}
		if (n_parts_in_cell[cell_i] >0)
		{
			mean_mfrac_in_depth_intervals[cell_i] = mfrac_in_depth_intervals_sum[cell_i] /
											double( n_parts_in_cell[cell_i] );
		}
		else
		{
			mean_mfrac_in_depth_intervals[cell_i] = -100;
		}
	}

	// now get the mean clay frac
	for (int cell_i=0; cell_i < n_cells; cell_i++)
	{
		v_iter = cfrac_in_depth_intervals[ cell_i ].begin();
		while (v_iter != cfrac_in_depth_intervals[ cell_i ].end())
		{
			cfrac_in_depth_intervals_sum[cell_i] += (*v_iter);
			v_iter++;
		}

		if (n_parts_in_cell[cell_i] >0)
		{
			mean_cfrac_in_depth_intervals[cell_i] = cfrac_in_depth_intervals_sum[cell_i] /
											double( n_parts_in_cell[cell_i] );
		}
		else
		{
			mean_cfrac_in_depth_intervals[cell_i] = -100;
		}
	}

	// now get the mean rate frac
	for (int cell_i=0; cell_i < n_cells; cell_i++)
	{
		v_iter = rfrac_in_depth_intervals[ cell_i ].begin();
		while (v_iter != rfrac_in_depth_intervals[ cell_i ].end())
		{
			rfrac_in_depth_intervals_sum[cell_i] += (*v_iter);
			v_iter++;
		}

		if (n_parts_in_cell[cell_i] >0)
		{
			mean_rfrac_in_depth_intervals[cell_i] = rfrac_in_depth_intervals_sum[cell_i] /
											double( n_parts_in_cell[cell_i] );
		}
		else
		{
			mean_rfrac_in_depth_intervals[cell_i] = -100;
		}
	}

	// now get the mean 10Be Concentration
	// NOTE: this takes no account of mass loss due to chemical weathering
	for (int cell_i=0; cell_i < n_cells; cell_i++)
	{
		v_iter = C10Be_in_depth_intervals[ cell_i ].begin();
		while (v_iter != C10Be_in_depth_intervals[ cell_i ].end())
		{
			C10Be_in_depth_intervals_sum[cell_i] += (*v_iter);
			v_iter++;
		}

		if (n_parts_in_cell[cell_i] >0)
		{
			mean_C10Be_in_depth_intervals[cell_i] = C10Be_in_depth_intervals_sum[cell_i] /
											double( n_parts_in_cell[cell_i] );
		}
		else
		{
			mean_C10Be_in_depth_intervals[cell_i] = -100;
		}
	}


	vtk_cell_out << "SCALARS p_in_cell int 1" << endl 
               << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << n_parts_in_cell[i] <<endl;
	}

	vtk_cell_out << "SCALARS mean_age float 1" << endl 
              << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << mean_age_in_depth_intervals[i] <<endl;
	}
	vtk_cell_out << "SCALARS mean_mass_frac float 1" << endl 
               << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << mean_mfrac_in_depth_intervals[i] <<endl;
	}
	vtk_cell_out << "SCALARS mean_clay_frac float 1" << endl 
               << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << mean_cfrac_in_depth_intervals[i] <<endl;
	}
	vtk_cell_out << "SCALARS mean_rate_frac float 1" << endl 
               << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << mean_rfrac_in_depth_intervals[i] <<endl;
	}
	vtk_cell_out << "SCALARS Conc_10Be float 1" << endl 
               << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_cells; i++)
	{
		vtk_cell_out << mean_C10Be_in_depth_intervals[i] <<endl;
	}

	vtk_particle_out << "# vtk DataFile Version 2.0" << endl << "Unstructured Grid Ptrack"
	        << endl << "ASCII" << endl << endl << "DATASET UNSTRUCTURED_GRID" << endl
	        << "POINTS " << n_parts << " float" << endl;

	if (reference_frame_switch == 1)
	{
		for (int i = 0; i< n_parts; i++)
		{
			vtk_particle_out << s_loc[i] << " " << -d_loc[i] << " 0.0" <<endl;
		}
	}
	else
	{
		for (int i = 0; i< n_parts; i++)
		{
			vtk_particle_out << s_loc[i] << " " << z_loc[i] << " 0.0" <<endl;
		}
	}

	vtk_particle_out << endl << "POINT_DATA "<<n_parts << endl 
                   << "SCALARS Particle_age float 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << pAge[i] <<endl;
	}

	vtk_particle_out << "SCALARS Be10_conc float 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << C10Be[i] <<endl;
	}

	vtk_particle_out << "SCALARS mfrac float 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << mass_fraction[i] <<endl;
	}

	vtk_particle_out << "SCALARS cfrac float 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << clay_fraction[i] <<endl;
	}

	vtk_particle_out << "SCALARS rfrac float 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << rate_fraction[i] <<endl;
	}

	vtk_particle_out << "SCALARS Ptype int 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
			vtk_particle_out << pID[i] <<endl;
	}

	vtk_particle_out << "SCALARS depth_interval_number int 1"
	        << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i< n_parts; i++)
	{
		vtk_particle_out << particle_in_depth_interval[i] <<endl;
	}

	// close the files
	vtk_cell_out.close();
	vtk_particle_out.close();


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// this function caluclates values within cells
//
// NUMBERING OF THE VERTICES AND CELLS
//
// the total number of cells is (n_depthintervals_soil+n_depthintervals_parent)*n_bins
// on each edge of a bin there are (n_depthintervals_soil+n_depthintervals_parent)+1 depth
// intervals
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRN_tParticle_bins::update_particles_cell_index(flowtube ft,
								int n_depthintervals_soil, int n_depthintervals_parent,
								double bottom_depth,
								vector<double>& verts_s, vector<double>& verts_z,
								vector<double>& verts_d,
								vector<int>& cell_node1, vector<int>& cell_node2,
								vector<int>& cell_node3, vector<int>& cell_node4)
{

	int total_depth_intervals = n_depthintervals_soil+n_depthintervals_parent;

	// vectors for holding the vertices
	// of the individual cells
	vector<double> empty_d_vec;
	verts_s = empty_d_vec;
	verts_z = empty_d_vec;
	verts_d = empty_d_vec;

	// this loops through each bin, getting the
	double us_zeta, us_eta;

	// get zeta, eta and s from the flowtube
	vector<double> eta = ft.get_eta();
	vector<double> zeta = ft.get_zeta();
	vector<double> s = ft.get_s_h();
	
	//for(int i = 0; i<s.size();i++)
	//{
  //  cout << "s[" << i << "]: " << s[i] << endl;
  //}

	// loop through the bins collecting data
	int n_bins = particle_bins.size();

	// the vectors hold the indices to the vertices
	vector<int> empty_i_vec;
	cell_node1 = empty_i_vec;
	cell_node2 = empty_i_vec;
	cell_node3 = empty_i_vec;
	cell_node4 = empty_i_vec;
	vector<int> cell_code;			// used to differentiate soil from parent material. 
                  // The numbers are
									// arbitrary, but for now 0 == parent and 1 == soil
	int soil_code = 1;
	int parent_code = 0;

	double us_soil_depthinterval,
	       us_parent_depthinterval;

	double s_up;
	double s_h_up;

	double fuzzy_boundary = 0.001;
	// now loop through the bins, calculating vertices
	for (int bn = 0; bn< n_bins; bn++)
	{
		s_h_up = s[ h_node_us[bn]];
		s_up = bin_edge_loc[bn];

		// if the upslope h location and the upslope bin locations are the same...
		if (s_h_up <= s_up+fuzzy_boundary &&  s_h_up >= s_up-fuzzy_boundary )
		{
			us_zeta = zeta[ h_node_us[bn] ];
			us_eta = eta[ h_node_us[bn] ];
		}
		// if not interpolate between the nodes
		else
		{
			us_zeta = ( (zeta[ h_node_us[bn] ]-zeta[ h_node_ds[bn] ])/dx_h[bn] )
			           *(s_up-s_h_up)+ zeta[ h_node_ds[bn] ];
			us_eta = ( (eta[ h_node_us[bn] ]-eta[ h_node_ds[bn] ])/dx_h[bn] )
			           *(s_up-s_h_up)+ eta[ h_node_ds[bn] ];
		}

		if (bn == 0)
		{
			us_zeta = zeta[ h_node_us[bn] ];
			us_eta = eta[ h_node_us[bn] ];
		}

		// calculate the depth of the bins in the soil
		us_soil_depthinterval = (us_zeta-us_eta)/double(n_depthintervals_soil);
		us_parent_depthinterval = (bottom_depth-(us_zeta-us_eta))/
                               double(n_depthintervals_parent);

		// first get the soil nodes
		for (int sdi = 0; sdi<n_depthintervals_soil; sdi++)
		{
			// get the s-coordinate vertex
			verts_s.push_back(bin_edge_loc[bn]);

      //cout << "The bin edge location: " <<  bin_edge_loc[bn] << "for bin: " << bn << endl;

			// get the vertical index
			verts_z.push_back( us_zeta - (double(sdi))*us_soil_depthinterval);

			// get the depth index
			verts_d.push_back( (double(sdi))*us_soil_depthinterval);

		}

		// now the parent material nodes
		for (int pdi = 0; pdi<n_depthintervals_parent; pdi++)
		{
			// get the s-coordinate vertex
			verts_s.push_back(bin_edge_loc[bn]);

			// get the vertical index
			verts_z.push_back( us_eta - (double(pdi))*us_parent_depthinterval);

			// get the depth index
			verts_d.push_back( (us_zeta-us_eta)+(double(pdi))*us_parent_depthinterval);
		}

		// and now the bottom node
		verts_s.push_back(bin_edge_loc[bn]);
		verts_z.push_back(us_zeta-bottom_depth);
		verts_d.push_back(bottom_depth);
	}

	// now get the last wall of the bin
	int lbn = n_bins-1;
  double s_h_ds = s[ h_node_ds[lbn]];
	double s_ds = bin_edge_loc[lbn+1];
	double ds_zeta = zeta[ h_node_ds[lbn] ];
	double ds_eta = eta[ h_node_ds[lbn] ];

	// calculate the depth of the bins in the soil
	double ds_soil_depthinterval = (ds_zeta-ds_eta)/double(n_depthintervals_soil);
	double ds_parent_depthinterval = (bottom_depth-(ds_zeta-us_eta))/
                             double(n_depthintervals_parent);

	// first get the soil nodes
	for (int sdi = 0; sdi<n_depthintervals_soil; sdi++)
	{
		// get the s-coordinate vertex
		verts_s.push_back(s_ds);

		// get the vertical index
		verts_z.push_back( ds_zeta - (double(sdi))*ds_soil_depthinterval);

		// get the depth index
		verts_d.push_back( (double(sdi))*ds_soil_depthinterval);

	}

	// now the parent material nodes
	for (int pdi = 0; pdi<n_depthintervals_parent; pdi++)
	{
		// get the s-coordinate vertex
		verts_s.push_back(s_ds);

		// get the vertical index
		verts_z.push_back( ds_eta - (double(pdi))*ds_parent_depthinterval);

		// get the depth index
		verts_d.push_back( (ds_zeta-ds_eta)+(double(pdi))*ds_parent_depthinterval);
	}

	// and now the bottom node
	verts_s.push_back(s_ds);
	verts_z.push_back(ds_zeta-bottom_depth);
	verts_d.push_back(bottom_depth);




	// loop through cells plotting vertices
	//for (int bn = 0; bn<n_bins-1; bn++)
	for (int bn = 0; bn<n_bins; bn++)
	{
		for (int di = 0; di<total_depth_intervals; di++)
		{
			cell_node1.push_back( bn*(total_depth_intervals+1)+di );
			cell_node2.push_back( (bn+1)*(total_depth_intervals+1)+di );
			cell_node3.push_back( (bn+1)*(total_depth_intervals+1)+di+1 );;
			cell_node4.push_back( bn*(total_depth_intervals+1)+di+1 );

			if (di < n_depthintervals_soil)
			{
				cell_code.push_back(soil_code);
			}
			else
			{
				cell_code.push_back(parent_code);
			}
		}
	}

	// now to do the particles
	list<CRN_tParticle>::iterator part_iter;	// list iterator

	int n_parts_in_bin;
	int depth_interval_index;

	int p_start = 0;			// index of the starting particle in bin
	int p_end = 0;			// index of the finishing particle in bin
	int pip_test;			// this is 1 if particle is in cell and 0 otherwise

	vector<double> sverts_for_pip(4);
	vector<double> zverts_for_pip(4);
	for (int bn = 0; bn< n_bins; bn++)
	{
		n_parts_in_bin = particle_bins[bn].size();
		p_start = p_end;
		p_end = p_start+n_parts_in_bin;

		part_iter = particle_bins[bn].begin();
		int counter = 0;


		for (int di = 0; di< total_depth_intervals; di++)
		{
			//cout << "depth_interval: " << di << " of " << total_depth_intervals << endl;

			// get the index for the cell
			depth_interval_index = bn*total_depth_intervals+di;

			// populate the vertices vectors
			sverts_for_pip[0] = verts_s[ cell_node1[depth_interval_index] ];
			sverts_for_pip[1] = verts_s[ cell_node2[depth_interval_index] ];
			sverts_for_pip[2] = verts_s[ cell_node3[depth_interval_index] ];
			sverts_for_pip[3] = verts_s[ cell_node4[depth_interval_index] ];

			zverts_for_pip[0] = verts_z[ cell_node1[depth_interval_index] ];
			zverts_for_pip[1] = verts_z[ cell_node2[depth_interval_index] ];
			zverts_for_pip[2] = verts_z[ cell_node3[depth_interval_index] ];
			zverts_for_pip[3] = verts_z[ cell_node4[depth_interval_index] ];

			// now loop through all the particles in this bin
			//cout << "p_start: " << p_start << " and p_end: " << p_end << endl;
			//cout << "number of parts: " << n_parts << endl;
			part_iter = particle_bins[bn].begin();
			while (part_iter != particle_bins[bn].end())
			{
        // check to see if the particle is within the cell
				pip_test = pnpoly(sverts_for_pip, zverts_for_pip, (*part_iter).getxLoc(), 
                           (*part_iter).get_zetaLoc());
				
        // if not, change the cell
        if (pip_test == 1)
				{
					(*part_iter).setCellIndex(depth_interval_index);
				}
				part_iter++;
			}
		}
	}
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function is purely for checking the cell index
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRN_tParticle_bins::check_particles_in_cells(int bn, vector<double>& verts_s,
                vector<double>& verts_z,
								vector<double>& verts_d,
								vector<int>& cell_node1, vector<int>& cell_node2,
								vector<int>& cell_node3, vector<int>& cell_node4)
{
	// now to do the particles
	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// parameters from the cells
	int cell_index;
	double zLoc;
	double dLoc;
	double sLoc;

	part_iter = particle_bins[bn].begin();
	int n_particles_in_bin = particle_bins[bn].size();
	cout << endl << endl << endl << " n_parts in bin: " 
       << n_particles_in_bin << endl;
	while (part_iter != particle_bins[bn].end())
	{
		cell_index = (*part_iter).getCellIndex();
		zLoc = (*part_iter).get_zetaLoc();
		sLoc = (*part_iter).getxLoc();
		dLoc = (*part_iter).getdLoc();

		// now print out the particle information
		cout << endl << "Particle zLoc: " << zLoc << " sLoc: " << sLoc 
         << " and d loc: " << dLoc
		     << " cell index: " << cell_index << endl;
		cout << "s1: " << verts_s[ cell_node1[cell_index] ]
		     << " z1: " << verts_z[ cell_node1[cell_index] ]
		     << " d1: " << verts_d[ cell_node1[cell_index] ] << endl;
		cout << "s2: " << verts_s[ cell_node2[cell_index] ]
		     << " z2: " << verts_z[ cell_node2[cell_index] ]
		     << " d2: " << verts_d[ cell_node2[cell_index] ] << endl;
		cout << "s3: " << verts_s[ cell_node3[cell_index] ]
		     << " z3: " << verts_z[ cell_node3[cell_index] ]
		     << " d3: " << verts_d[ cell_node3[cell_index] ] << endl;
		cout << "s4: " << verts_s[ cell_node4[cell_index] ]
		     << " z4: " << verts_z[ cell_node4[cell_index] ]
		     << " d4: " << verts_d[ cell_node4[cell_index] ] << endl;

		part_iter++;
	}

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function is purely for checking the cell index
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRN_tParticle_bins::get_data_by_cell(int bn, int n_PDZ_intervals, int n_CAZ_intervals,
									double bottom_depth, vector<double> verts_s, vector<double> verts_z,
								vector<double>& verts_d,
								vector<int>& cell_node1, vector<int>& cell_node2,
								vector<int>& cell_node3, vector<int>& cell_node4)

{
	int n_cells = (n_PDZ_intervals+n_CAZ_intervals)*n_bins;
	int starting_cell = bn*(n_PDZ_intervals+n_CAZ_intervals);
	int ending_cell = (bn+1)*(n_PDZ_intervals+n_CAZ_intervals);
	vector<int> particles_in_cells(n_cells,0);

	// now to do the particles
	list<CRN_tParticle>::iterator part_iter;	// list iterator

	// parameters from the cells
	int cell_index;
	double zLoc;
	double dLoc;
	double sLoc;

	part_iter = particle_bins[bn].begin();
	int n_particles_in_bin = particle_bins[bn].size();
	cout << endl << endl << endl << " n_parts in bin: " << n_particles_in_bin << endl;
	while (part_iter != particle_bins[bn].end())
	{
		cell_index = (*part_iter).getCellIndex();
		zLoc = (*part_iter).get_zetaLoc();
		sLoc = (*part_iter).getxLoc();
		dLoc = (*part_iter).getdLoc();

		particles_in_cells[cell_index]++;

		part_iter++;
	}

	for(int ci = starting_cell; ci< ending_cell; ci++)
	{
		cout << "cell number: " << ci << endl;
		cout << "n parts in cell: " << particles_in_cells[ci] << endl;
		cout << "s1: " << verts_s[ cell_node1[ci] ]
		     << " z1: " << verts_z[ cell_node1[ci] ]
		     << " d1: " << verts_d[ cell_node1[ci] ] << endl;
		cout << "s2: " << verts_s[ cell_node2[ci] ]
		     << " z2: " << verts_z[ cell_node2[ci] ]
		     << " d2: " << verts_d[ cell_node2[ci] ] << endl;
		cout << "s3: " << verts_s[ cell_node3[ci] ]
		     << " z3: " << verts_z[ cell_node3[ci] ]
		     << " d3: " << verts_d[ cell_node3[ci] ] << endl;
		cout << "s4: " << verts_s[ cell_node4[ci] ]
		     << " z4: " << verts_z[ cell_node4[ci] ]
		     << " d4: " << verts_d[ cell_node4[ci] ] << endl;

	}
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function loops through the bin getting mass fractions and depletions
// It assumes the cell indices and verts have been calcualted already
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRN_tParticle_bins::get_mineral_mass_loss_and_mfracs_volumetric(int bn, int n_PDZ_intervals, 
                  int n_CAZ_intervals,
									VolumeParticleInfo& vpi,
									list< vector<double> >& mineral_mfracs,
									list< vector<double> >& mineral_depletion)
{
	// reset the list vectors
	list< vector<double> > empty_list_vec;
	mineral_mfracs = empty_list_vec;
	mineral_depletion = empty_list_vec;

	// get some cell information
	int n_cells = (n_PDZ_intervals+n_CAZ_intervals)*n_bins;
	int starting_cell = bn*(n_PDZ_intervals+n_CAZ_intervals);
	int ending_cell = (bn+1)*(n_PDZ_intervals+n_CAZ_intervals);

	// initiate data storage
	vector<double> mass_in_cell(n_cells,0.0);
	vector<double> empty_vec(n_cells,0.0);
	list< vector<double> > mass_of_types_in_cells;
	list< vector<double> > starting_mass_of_types_in_cells;

	// set up empty list vecs
	int n_types = vpi.get_n_types();
	for (int type = 0; type<n_types; type++)
	{
		mass_of_types_in_cells.push_back(empty_vec);
		starting_mass_of_types_in_cells.push_back(empty_vec);
		mineral_depletion.push_back(empty_vec);
		mineral_mfracs.push_back(empty_vec);
	}

	// now to do the particles
	list<CRN_tParticle>::iterator part_iter;	// list iterator
	list< vector<double> >::iterator mass_iter;
	list< vector<double> >::iterator starting_mass_iter;
	list< vector<double> >::iterator mfrac_iter;
	list< vector<double> >::iterator depletion_iter;	

	// parameters from the cells
	int cell_index;
	int particleType;
	double particleMass;
	double particleStartingMass;
	double thick_upslope,thick_downslope;
	double cell_volume;

	// loop through bins collecting particles in different cells
	part_iter = particle_bins[bn].begin();
	int n_particles_in_bin = particle_bins[bn].size();
	while (part_iter != particle_bins[bn].end())
	{
	  // get the properties of this particle
		cell_index = (*part_iter).getCellIndex();
		particleType = (*part_iter).getType();
		particleMass = (*part_iter).getMass();
		particleStartingMass = (*part_iter).getStartingMass();

		// collect data elements based on individual particles
		if (cell_index >= 0)
		{
			mass_in_cell[cell_index]+=particleMass;

			mass_iter = mass_of_types_in_cells.begin();
			starting_mass_iter = starting_mass_of_types_in_cells.begin();

			for (int count = 0; count < particleType; count++)
			{
				starting_mass_iter++;
				mass_iter++;
			}

			// add to running total of starting mass and
      // mass for different particle types.
			(*mass_iter)[cell_index]+=particleMass;
			(*starting_mass_iter)[cell_index]+=particleStartingMass;
		}

		part_iter++;
	}

  // now loop through the cells, aggregating the information
	for(int ci = starting_cell; ci< ending_cell; ci++)
	{
	  // get the iterators. These refer to the element on the list
	  // corresponding to the type of interest. 
		mass_iter = mass_of_types_in_cells.begin();
    starting_mass_iter = starting_mass_of_types_in_cells.begin();
    mfrac_iter = mineral_mfracs.begin();
    depletion_iter = mineral_depletion.begin();
    
    // loop through the types. After collecting cell information
    // of all the cells it moves on to the next type by incrementing
    // the iterators
		for (int i = 0; i<n_types; i++)
		{
			if ( (*mass_iter)[ci] == 0)
			{
				(*mfrac_iter)[ci] = 0;
				(*depletion_iter)[ci] = 1;
			}
			else
			{
				(*mfrac_iter)[ci] = (*mass_iter)[ci]/mass_in_cell[ci];
				(*depletion_iter)[ci] = ((*mass_iter)[ci]/( (*starting_mass_iter)[ci] ))-1;
			}
			mass_iter++;
			starting_mass_iter++;
			mfrac_iter++;
			depletion_iter++;
		}
	}
	
  // note: at this stage we do not need to update the listvecs that were
  // passed to the function since they have already been updated in 
  // the loops above
}
									

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function is used to collect data about the volume fractions and 
// specific surface areas in
// cells that are then fed into CRUNCHflow
//
// it creates some list vecs where the vector elements are for the cells (i.e. the boxes in which the
// particles reside) and the list elements are for the different types.
//
// !!!NOTE!!! This function resets the list vecs and returns new ones.
// It is an extrmeley stupid way to go about this since the vectors are storing
// lots of zeros which are not ever used, but I don't have the time to 
// fix this. SMM 09/08/2014
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRN_tParticle_bins::get_data_by_cell_volumetric_for_CRUNCH(int bn, int n_PDZ_intervals, 
                  int n_CAZ_intervals,
									double bottom_depth, vector<double> verts_s, vector<double> verts_z,
									vector<double>& verts_d,
									vector<int>& cell_node1, vector<int>& cell_node2,
									vector<int>& cell_node3, vector<int>& cell_node4,
									VolumeParticleInfo vpi,
									list< vector<double> >& mineral_vpercents,
									list< vector<double> >& mineral_ssa,
									list< vector<double> >& mineral_surface_area,
									list< vector<double> >& mineral_mass)

{

	// reset the list vectors
	list< vector<double> > empty_list_vec;
	mineral_vpercents = empty_list_vec;
	mineral_ssa = empty_list_vec;

	// get some cell information
	int n_cells = (n_PDZ_intervals+n_CAZ_intervals)*n_bins;
	int starting_cell = bn*(n_PDZ_intervals+n_CAZ_intervals);
	int ending_cell = (bn+1)*(n_PDZ_intervals+n_CAZ_intervals);
	vector<int> particles_in_cells(n_cells,0);

	// initiate data storage
	vector<double> mass_in_cell(n_cells,0.0);
	vector<double> volume_of_particles_in_cell(n_cells,0.0);
	vector<double> porosity_in_cell(n_cells,0.0);
	vector<double> empty_vec(n_cells,0.0);
	list< vector<double> > surface_areas_of_types_in_cells;
	list< vector<double> > mass_of_types_in_cells;
	list< vector<double> > volumes_of_types_in_cells;

	list< vector<double> > ssa_of_types_in_cells;
	list< vector<double> > volume_percents_of_types_in_cells;

	// set up empty list vecs
	int n_types = vpi.get_n_types();
	for (int type = 0; type<n_types; type++)
	{
		surface_areas_of_types_in_cells.push_back(empty_vec);
		volumes_of_types_in_cells.push_back(empty_vec);
		mass_of_types_in_cells.push_back(empty_vec);
		ssa_of_types_in_cells.push_back(empty_vec);
		volume_percents_of_types_in_cells.push_back(empty_vec);
	}

	// now to do the particles
	list<CRN_tParticle>::iterator part_iter;	// list iterator
	list< vector<double> >::iterator surf_area_iter;
	list< vector<double> >::iterator volume_iter;
	list< vector<double> >::iterator mass_iter;
	list< vector<double> >::iterator ssa_iter;
	list< vector<double> >::iterator volume_percents_iter;

	// parameters from the cells
	int cell_index;
	int particleType;
	double zLoc;
	double dLoc;
	double sLoc;
	double particleVolume;
	double particleMass;
	double particleSurfaceArea;
	double thick_upslope,thick_downslope;
	double cell_volume;

	// loop through bins collecting particles in different cells
	part_iter = particle_bins[bn].begin();
	int n_particles_in_bin = particle_bins[bn].size();
	//cout << endl << endl << endl << " n_parts in bin: " << n_particles_in_bin << endl;
	//cout << "LINE 4202 First part, zl: " << (*part_iter).get_zetaLoc() << " sl: " << (*part_iter).getxLoc()
	//     << " dL: " << (*part_iter).getdLoc() << endl;
	while (part_iter != particle_bins[bn].end())
	{
		cell_index = (*part_iter).getCellIndex();
		zLoc = (*part_iter).get_zetaLoc();
		sLoc = (*part_iter).getxLoc();
		dLoc = (*part_iter).getdLoc();

		particleType = (*part_iter).getType();
		particleVolume = (*part_iter).update_surface_area_and_get_volume(vpi);
		particleMass = (*part_iter).getMass();
		particleSurfaceArea =  (*part_iter).getSurfaceArea();

		//cout << "LINE 4091: type: " << particleType << " cell index: " << cell_index << endl
		//     <<" Vol:" << particleVolume << " Mass: " << particleMass << " SA: " << particleSurfaceArea << endl;
		//cout << "zloc: " << zLoc << " sLoc: " << sLoc << " and dLoc: " << dLoc << endl;

		// collect data elements based on individual particles
		if (cell_index >= 0)
		{
			particles_in_cells[cell_index]++;
			mass_in_cell[cell_index]+=particleMass;
			volume_of_particles_in_cell[cell_index]+=particleVolume;

			surf_area_iter = surface_areas_of_types_in_cells.begin();
			volume_iter = volumes_of_types_in_cells.begin();
			mass_iter = mass_of_types_in_cells.begin();

			for (int count = 0; count < particleType; count++)
			{
				surf_area_iter++;
				volume_iter++;
				mass_iter++;
			}

			// add to running total of surface area, volume and 
      // mass for different particle types.
			(*surf_area_iter)[cell_index]+=particleSurfaceArea;
			(*volume_iter)[cell_index]+=particleVolume;
			(*mass_iter)[cell_index]+=particleMass;

		}

		part_iter++;
	}

	//cout << "starting cell: " << starting_cell << " ending cell: " << ending_cell << endl;
	//cout << "cell_node size: " << cell_node1.size() << " " << cell_node2.size() << " "
	//     << cell_node3.size() << " " << cell_node4.size() << endl;

	for(int ci = starting_cell; ci< ending_cell; ci++)
	{
		// calculate the thicknesses of the upslope and downslope sections of the cells
		thick_upslope = verts_d[ cell_node4[ci] ] - verts_d[ cell_node1[ci] ];
		thick_downslope = verts_d[ cell_node3[ci] ] - verts_d[ cell_node2[ci] ];


		// now calcualte the volume of the cell
		if (thick_upslope > thick_downslope)
		{
			cell_volume = A_bins[bn]*(thick_downslope + 0.5*(thick_upslope-thick_downslope));
		}
		else if (thick_upslope < thick_downslope)
		{
			cell_volume = A_bins[bn]*(thick_upslope + 0.5*(thick_downslope-thick_upslope));
		}
		else
		{
			cell_volume = A_bins[bn]*thick_upslope;
		}

		porosity_in_cell[ci] = volume_of_particles_in_cell[ci]/cell_volume;


		//cout << endl << "cell number: " << ci << endl;
		//cout << "thick upslope: " << thick_upslope << " thick_downslope: "
         // << thick_downslope << " area: " << A_bins[bn] << endl;
		//cout << "n parts in cell: " << particles_in_cells[ci] << " " 
    //     << "volume of cell: " << cell_volume << endl;
		//cout << "volume of parts: " << volume_of_particles_in_cell[ci] 
    //     << " mass of particles: " << mass_in_cell[ci] << endl;
		//cout << "porosity: " << porosity_in_cell[ci] << " and bulk density: " 
    //    << mass_in_cell[ci]/cell_volume << endl;
		//cout << "s1: " << verts_s[ cell_node1[ci] ]
		//     << " z1: " << verts_z[ cell_node1[ci] ]
		//     << " d1: " << verts_d[ cell_node1[ci] ] << endl;
		//cout << "s2: " << verts_s[ cell_node2[ci] ]
		//     << " z2: " << verts_z[ cell_node2[ci] ]
		//     << " d2: " << verts_d[ cell_node2[ci] ] << endl;
		//cout << "s3: " << verts_s[ cell_node3[ci] ]
		//     << " z3: " << verts_z[ cell_node3[ci] ]
		//     << " d3: " << verts_d[ cell_node3[ci] ] << endl;
		//cout << "s4: " << verts_s[ cell_node4[ci] ]
		//     << " z4: " << verts_z[ cell_node4[ci] ]
		//     << " d4: " << verts_d[ cell_node4[ci] ] << endl;


		surf_area_iter = surface_areas_of_types_in_cells.begin();
		volume_iter = volumes_of_types_in_cells.begin();
		mass_iter = mass_of_types_in_cells.begin();
		/*
    cout << "mass fractions: " << endl;
		for (int i = 0; i<n_types; i++)
		{
			cout << " " << (*mass_iter)[ci]/mass_in_cell[ci];
			mass_iter++;
		}
		cout << endl;

		cout << "volume fractions: " << endl;
		for (int i = 0; i<n_types; i++)
		{
			cout << " " << (*volume_iter)[ci]/cell_volume;
			volume_iter++;
		}
		cout << endl;

		cout << "surface area: " << endl;
		for (int i = 0; i<n_types; i++)
		{
			cout << " " << (*surf_area_iter)[ci];
			surf_area_iter++;
		}
		cout << endl;
    */
    
		ssa_iter = ssa_of_types_in_cells.begin();
		mass_iter = mass_of_types_in_cells.begin();
		surf_area_iter = surface_areas_of_types_in_cells.begin();
		//cout << "ssa: " << endl;
		for (int i = 0; i<n_types; i++)
		{
			if ( (*mass_iter)[ci] == 0)
			{
				//cout << " 0";
				(*ssa_iter)[ci] = 0;
			}
			else
			{
				// the factor of 0.001 is because crunch takes ssa in m^2/g but
				// mass is stored in kg
				//cout << " " << 0.001*(*surf_area_iter)[ci]/(*mass_iter)[ci];
				(*ssa_iter)[ci] = 0.001*(*surf_area_iter)[ci]/(*mass_iter)[ci];
			}
			mass_iter++;
			surf_area_iter++;
			ssa_iter++;
		}
		//cout << endl;

		volume_percents_iter = volume_percents_of_types_in_cells.begin();
		volume_iter = volumes_of_types_in_cells.begin();
		for (int i = 0; i<n_types; i++)
		{
			// the 100 factor is because it is in a percent
			(*volume_percents_iter)[ci] = 100*(*volume_iter)[ci]/cell_volume;
			volume_iter++;
			volume_percents_iter++;
		}

	}
	
	

	// update the listvecs
	mineral_vpercents = volume_percents_of_types_in_cells;
	mineral_ssa = ssa_of_types_in_cells;
	mineral_surface_area = surface_areas_of_types_in_cells;
	mineral_mass = mass_of_types_in_cells;
	//cout << "LINE 4138 finished getting particle data" << endl;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-




//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function weathers the particles.
// it takes two list_vecs: both of which are the volume fractions of the particles
// the two files are volum fraction before and after. The volume fraction 
// is per volume porous media
//
// For each particle type, then, this function calculates a loss or gain 
// of that type of particle.
// The algorithm then needs to calcualte the corresponding change 
// in mass as well as the total surface area
// of the particle.
//
// the mass loss for each particle is then divided amongst 
// the particles of different size classes:
// a mass loss per surafce area is calucalted followed and 
// thus the mass loss is calcualted by
// multiplying the surface area of a given particle type by its total surface area
//
// ***WARNING*** CRUNCH files are only output with precision of 6, 
// so there are rounding errors.
// these in fact can artificially increase or decrease weathering, 
// and for slowly weathering
// minerals over 1 year the error is the same magnitude as the weathering rate
//
// FIX: the weathering rate is calculated based _only_ 
// upon the precision returned by CRUNCH
//
// NOTE The list_vecs passed to this function store vectors that have indices to 
// ALL the cells, not just those cells in the bin. This is actually extremely stupid 
// since the vectors are replaced each timestep, but it will take too long at the moment
// to fix this!
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void CRN_tParticle_bins::weather_particles_from_CRUNCH(int bn, int n_PDZ_intervals, 
                  int n_CAZ_intervals,
									double bottom_depth, vector<double> verts_s, vector<double> verts_z,
									vector<double>& verts_d,
									vector<int>& cell_node1, vector<int>& cell_node2,
									vector<int>& cell_node3, vector<int>& cell_node4,
									VolumeParticleInfo vpi,
									list< vector<double> >& mineral_vpercents_old,
                  list< vector<double> >& mineral_vpercents_new,
                  list< vector<double> >& surface_area_in_cell_old,
                  list< vector<double> >& mass_in_cell_old)
{
	// set up some variables
	double thick_upslope;
	double thick_downslope;
	double cell_volume;
	double vf_in_CRUNCH_precision_old, vf_in_CRUNCH_precision_new;
	double vf_change_in_CRUNCH_precision;
	double mass_change;

	// set up iterators
	list< vector<double> >::iterator mvp_old_iter;
	list< vector<double> >::iterator mvp_new_iter;
	list< vector<double> >::iterator surface_area_old_iter;
	list< vector<double> >::iterator loss_per_surface_area_iter;
	list<CRN_tParticle>::iterator    part_iter;	// list iterator

	// get some information about the cells so that the list vecs can be built  
	int n_cells = (n_PDZ_intervals+n_CAZ_intervals)*n_bins;
	int starting_cell = bn*(n_PDZ_intervals+n_CAZ_intervals);
	int ending_cell = (bn+1)*(n_PDZ_intervals+n_CAZ_intervals);

	vector<double> empty_vec(n_cells,0.0);
	list< vector<double> > loss_per_surface_area;

	// set up empty list vecs
	int n_types = vpi.get_n_types();
	for (int type = 0; type<n_types; type++)
	{
		loss_per_surface_area.push_back(empty_vec);
	}

	// now go through the bin finding the change in volume and
	// the change in mass within the cells
	// cell nodes are the coordinates of the cell corners: there are 4
	for(int ci = starting_cell; ci< ending_cell; ci++)
	{
		// calculate the thicknesses of the upslope and downslope sections of the cells
		thick_upslope = verts_d[ cell_node4[ci] ] - verts_d[ cell_node1[ci] ];
		thick_downslope = verts_d[ cell_node3[ci] ] - verts_d[ cell_node2[ci] ];

		// now calcualte the volume of the cell
		if (thick_upslope > thick_downslope)
		{
			cell_volume = A_bins[bn]*(thick_downslope + 0.5*(thick_upslope-thick_downslope));
		}
		else if (thick_upslope < thick_downslope)
		{
			cell_volume = A_bins[bn]*(thick_upslope + 0.5*(thick_downslope-thick_upslope));
		}
		else
		{
			cell_volume = A_bins[bn]*thick_upslope;
		}

		// this loops through each of the types in each cell and
		// returns the mass lost per surface area
		mvp_old_iter = mineral_vpercents_old.begin();
		mvp_new_iter = mineral_vpercents_new.begin();
		surface_area_old_iter = surface_area_in_cell_old.begin();
		loss_per_surface_area_iter = loss_per_surface_area.begin();
		cout.precision(12);
		//cout << endl << "cell: " << ci << " volume percents; type | old | new" << endl;
		//cout << "volume in cell: " << cell_volume << endl;
		for (int i = 0; i<n_types; i++)
		{
			//cout << "density: " << vpi.get_type_density(i) << endl;

			// the mass in a cell may be calculated quite exactly, 
      // but CRUNCHflow only takes a volume fraction data element
			// with a precision of 7, and only reports a volume fraction
      // with a precision of 7. So here I have to convert
			// the volume percent to something that is a volume fraction 
      // with a precision of 7 and calcualte
			// the rate based on this mass
			// there also has to be a very annoying series of if else 
      // statements to acocunt for the fact the CRUNCHflow
			// does things in scientific notation so if the volume percent is, 
      // say 3, the logiuc has to be different
			// than for a volume percent of 10.
			if ((*mvp_old_iter)[ci] == 100)
			{
				vf_in_CRUNCH_precision_old = floor((*mvp_old_iter)[ci]*10000+0.5) * 0.000001;
			}
			else if ((*mvp_old_iter)[ci] < 100 && (*mvp_old_iter)[ci] >= 10)
			{
				vf_in_CRUNCH_precision_old = floor((*mvp_old_iter)[ci]*100000+0.5) * 0.0000001;
			}
			else if ((*mvp_old_iter)[ci] < 10 && (*mvp_old_iter)[ci] >= 1)
			{
				vf_in_CRUNCH_precision_old = floor((*mvp_old_iter)[ci]*1000000+0.5) * 0.00000001;
			}
			else if ((*mvp_old_iter)[ci] < 1 && (*mvp_old_iter)[ci] >= 0.1)
			{
				vf_in_CRUNCH_precision_old = floor((*mvp_old_iter)[ci]*10000000+0.5) * 0.000000001;
			}
			else if ( (*mvp_old_iter)[ci] < 0.1 )
			{
				vf_in_CRUNCH_precision_old = 0;
			}

			if ((*mvp_new_iter)[ci] == 100)
			{
				vf_in_CRUNCH_precision_new = floor((*mvp_new_iter)[ci]*10000+0.5) * 0.000001;
			}
			else if ((*mvp_new_iter)[ci] < 100 && (*mvp_new_iter)[ci] >= 10)
			{
				vf_in_CRUNCH_precision_new = floor((*mvp_new_iter)[ci]*100000+0.5) * 0.0000001;
			}
			else if ((*mvp_new_iter)[ci] < 10 && (*mvp_new_iter)[ci] >= 1)
			{
				vf_in_CRUNCH_precision_new = floor((*mvp_new_iter)[ci]*1000000+0.5) * 0.00000001;
			}
			else if ((*mvp_new_iter)[ci] < 1 && (*mvp_new_iter)[ci] >= 0.1)
			{
				vf_in_CRUNCH_precision_new = floor((*mvp_new_iter)[ci]*10000000+0.5) * 0.000000001;
			}
			else if ( (*mvp_new_iter)[ci] < 0.1 )
			{
				vf_in_CRUNCH_precision_new = 0;
			}
			vf_change_in_CRUNCH_precision = vf_in_CRUNCH_precision_new-vf_in_CRUNCH_precision_old;
			mass_change = vf_change_in_CRUNCH_precision*cell_volume*vpi.get_type_density(i);

			//cout << "vf Cp old: " << vf_in_CRUNCH_precision_old << " and vf Cp new: " 
      //     << vf_in_CRUNCH_precision_new
			//     << " and vf_change: " << vf_change_in_CRUNCH_precision << endl;

			if ( (*surface_area_old_iter)[ci] <=0)
			{
				(*loss_per_surface_area_iter)[ci] = 0;
			}
			else
			{
				(*loss_per_surface_area_iter)[ci] = mass_change/(*surface_area_old_iter)[ci];
			}

			// << "mass loss: " << mass_change
			//     << " and per surf area loss: " << (*loss_per_surface_area_iter)[ci] 
      //     << " and surf area: " << (*surface_area_old_iter)[ci] << endl;

			//cout << "mass loss in CRUNCH precision: " << mass_in_CRUNCH_precision - 
      // (*mass_old_iter)[ci] <<endl;

			mvp_old_iter++;
			mvp_new_iter++;
			surface_area_old_iter++;
			loss_per_surface_area_iter++;

		}
		//cout << endl;
	}

	// now loop through all the particles, removing mass based on their surface areas
	double massLoss;
	part_iter = particle_bins[bn].begin();
	int n_particles_in_bin = particle_bins[bn].size();
	//cout << endl << endl << endl << " n_parts in bin: " << n_particles_in_bin << endl;
	while (part_iter != particle_bins[bn].end())
	{
		massLoss = (*part_iter).weather_particle(vpi,loss_per_surface_area);
		part_iter++;
	}
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function retrieves copies of particles that are in a sampling interval
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
list<CRN_tParticle> CRN_tParticle_bins::get_particles_in_depth_interval(int bn, 
                double d_top, double d_bottom)
{
	// list that evenutally gets returns
	list<CRN_tParticle> sampled_list;

	// particle iterator
	list<CRN_tParticle>::iterator part_iter;

	double particle_d;
	// loop through the bin bn collecting particles within the depth interval
	part_iter = particle_bins[bn].begin();
	while (part_iter != particle_bins[bn].end())
	{
		particle_d = (*part_iter).getdLoc();
		if (particle_d > d_top && particle_d <= d_bottom)
		{
			CRN_tParticle sampled_part = (*part_iter);
			sampled_list.push_back(sampled_part);
		}
		part_iter++;
	}
	return sampled_list;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this loops through vectors holding the top and bottom locations of the samples
// and calcualtes the averages of their parameter values
void CRN_tParticle_bins::calculate_sample_averages(int bn, 
                   vector<double>& d_top_locs, vector<double>& d_bottom_locs,
									 vector<double>& sample_mean_age, vector<double>& sample_mean_C10Be,
									 vector<double>& sample_mean_Cf10Be)
{
	// get the number of samples
	int n_samples = d_top_locs.size();

	// initiate some variables for calculating averages
	list<CRN_tParticle> sample_list;
	double d_top,d_bottom;
	list<CRN_tParticle>::iterator list_iter;

	// an empty vector for resetting mean vectors
	vector<double> empty;
	sample_mean_age = empty;
	sample_mean_C10Be = empty;
	sample_mean_Cf10Be = empty;

	int n_particles;
	double tot_age, tot_C10Be, tot_Cf10Be;

	// loop through the samples
	for (int i = 0; i< n_samples; i++)
	{
		// reset the totals
		tot_age = 0;
		tot_C10Be = 0;
		tot_Cf10Be = 0;

		// get the sample locations
		d_top = d_top_locs[i];
		d_bottom = d_bottom_locs[i];
		sample_list = get_particles_in_depth_interval(bn, d_top, d_bottom);

		// get the number of particles in the sample
		n_particles = sample_list.size();

		// now collect data from each sample
		list_iter = sample_list.begin();
		while (list_iter != sample_list.end())
		{
			tot_age += (*list_iter).getAge();
			tot_C10Be += (*list_iter).getConc_10Be();
			tot_Cf10Be += (*list_iter).getConc_f10Be();

			list_iter++;
		}

		sample_mean_age.push_back( tot_age/double(n_particles) );
		sample_mean_C10Be.push_back( tot_C10Be/double(n_particles) );
		sample_mean_Cf10Be.push_back( tot_Cf10Be/double(n_particles) );
	}
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
void CRN_tParticle_bins::calculate_sample_averages(vector<double>& s_locs,
								vector<double>& d_top_locs, vector<double>& d_bottom_locs,
								Particle_info& pi,
								vector<double>& sample_mean_age, vector<double>& sample_mean_C10Be,
								vector<double>& sample_enrich, vector<double>& sample_frac_clay,
								vector<double>& sample_frac_t1, vector<double>& sample_frac_t2,
								vector<double>& sample_frac_t3, vector<double>& sample_frac_t4,
								vector<double>& sample_frac_t5)
{
	// get the number of samples
	int n_samples = d_top_locs.size();
	int bn;
	int type;
	vector<double> weath_info;

	// initiate some variables for calculating averages
	list<CRN_tParticle> sample_list;
	double d_top,d_bottom;
	list<CRN_tParticle>::iterator list_iter;

	// an empty vector for resetting mean vectors
	vector<double> empty;
	sample_mean_age = empty;
	sample_mean_C10Be = empty;
	sample_enrich = empty;
	sample_frac_clay = empty;
	sample_frac_t1 = empty;
	sample_frac_t2 = empty;
	sample_frac_t3 = empty;
	sample_frac_t4 = empty;
	sample_frac_t5 = empty;

	int n_particles;
	double tot_age, tot_C10Be, tot_mr1, tot_mr2, tot_mr3, tot_mr4, tot_mr5,
			tot_pf1, tot_pf2, tot_pf3, tot_pf4, tot_pf5,
			tot_cf1, tot_cf2, tot_cf3, tot_cf4, tot_cf5,
			tot_enrich;
	int nt_1, nt_2, nt_3, nt_4, nt_5;

	// loop through the samples
	for (int i = 0; i< n_samples; i++)
	{
		bn = get_bin_number_of_sloc(s_locs[i]);

		// reset the totals
		tot_age = 0;
		tot_C10Be = 0;
		tot_mr1 = 0;
		tot_mr2 = 0;
		tot_mr3 = 0;
		tot_mr4 = 0;
		tot_mr5 = 0;
		tot_pf1 = 0;
		tot_pf2 = 0;
		tot_pf3 = 0;
		tot_pf4 = 0;
		tot_pf5 = 0;
		tot_cf1 = 0;
		tot_cf2 = 0;
		tot_cf3 = 0;
		tot_cf4 = 0;
		tot_cf5 = 0;
		tot_enrich = 0;
		nt_1 = 0;
		nt_2 = 0;
		nt_3 = 0;
		nt_4 = 0;
		nt_5 = 0;

		// get the sample locations
		d_top = d_top_locs[i];
		d_bottom = d_bottom_locs[i];
		sample_list = get_particles_in_depth_interval(bn, d_top, d_bottom);



		// get the number of particles in the sample
		n_particles = sample_list.size();

		//cout << "s: " << s_locs[i] << " bin number: " << bn << " d_top: " << d_top
		//     << " d bottom: " << d_bottom << " n_parts: " << n_particles << endl << endl;

		// now collect data from each sample
		list_iter = sample_list.begin();
		while (list_iter != sample_list.end())
		{
			tot_age += (*list_iter).getAge();
			tot_C10Be += (*list_iter).getConc_10Be();
			type = (*list_iter).getType();
			weath_info = pi.calculate_weathering_on_demand((*list_iter).getAge(),type);

			if (type == 1)
			{
				nt_1++;
				tot_mr1 += weath_info[3];
				tot_pf1 += weath_info[0];
				tot_cf1 += weath_info[1];
			}
			else if (type == 2)
			{
				nt_2++;
				tot_mr2 += weath_info[3];
				tot_pf2 += weath_info[0];
				tot_cf2 += weath_info[1];
			}
			else if (type == 3)
			{
				nt_3++;
				tot_mr3 += weath_info[3];
				tot_pf3 += weath_info[0];
				tot_cf3 += weath_info[1];
			}
			else if (type == 4)
			{
				nt_4++;
				tot_mr4 += weath_info[3];
				tot_pf4 += weath_info[0];
				tot_cf4 += weath_info[1];
			}
			else if (type == 5)
			{
				nt_5++;
				tot_mr5 += weath_info[3];
				tot_pf5 += weath_info[0];
				tot_cf5 += weath_info[1];
			}

			list_iter++;
		}

		double total_p_all = tot_pf1+tot_pf2+tot_pf3+tot_pf4+tot_pf5;
		double total_c_all = tot_cf1+tot_cf2+tot_cf3+tot_cf4+tot_cf5;
		double total_m_rem = tot_mr1+tot_mr2+tot_mr3+tot_mr4+tot_mr5;
		sample_mean_age.push_back( tot_age/double(n_particles) );
		sample_mean_C10Be.push_back( tot_C10Be/double(n_particles) );
		sample_enrich.push_back( double(n_particles)/total_m_rem );
		sample_frac_clay.push_back(total_c_all/total_m_rem );
		sample_frac_t1.push_back(tot_pf1/total_m_rem );
		sample_frac_t2.push_back(tot_pf2/total_m_rem );
		sample_frac_t3.push_back(tot_pf3/total_m_rem );
		sample_frac_t4.push_back(tot_pf4/total_m_rem );
		sample_frac_t5.push_back(tot_pf5/total_m_rem );


	}

	//pi.details_to_screen();
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// this function takes PDZ and CAZ intervals and creates vectors of top and
// bottom depths for particle bins.
// Note there will be some discrepancy along the soil-saprolite boundary because the
// thickness of the soil changes in the downslope direction but the depth interval
// is only calcualted at the 'pit point'. However, by choosing a small dx this
// error can be minimized.
void CRN_tParticle_bins::partition_bins_into_cells(int bn, flowtube ft,
							int n_PDZ_intervals, int n_CAZ_intervals,
									double bottom_depth,
								   vector<double>& d_top_locs,
								   vector<double>& d_bottom_locs)
{

	//int total_depth_intervals = n_depthintervals_soil+n_depthintervals_parent;

	// get zeta, eta and s from the flowtube
	vector<double> eta = ft.get_eta();
	vector<double> zeta = ft.get_zeta();
	vector<double> s = ft.get_s_h();

	double fuzzy_boundary = 0.001;
	// now loop through the bins, calcualting vertices
	double s_h_up = s[ h_node_us[bn]];
	double s_up = bin_edge_loc[bn];
	double us_zeta;
	double us_eta;

	vector<double> empty_vec;
	d_top_locs = empty_vec;
	d_bottom_locs = empty_vec;

	// if the upslope h location and the upslope bin locations are the same...
	if (s_h_up <= s_up+fuzzy_boundary &&  s_h_up >= s_up-fuzzy_boundary )
	{
		us_zeta = zeta[ h_node_us[bn] ];
		us_eta = eta[ h_node_us[bn] ];
	}
	// if not interpolate between the nodes
	else
	{
		us_zeta = ( (zeta[ h_node_us[bn] ]-zeta[ h_node_ds[bn] ])/dx_h[bn] )
				   *(s_up-s_h_up)+ zeta[ h_node_ds[bn] ];
		us_eta = ( (eta[ h_node_us[bn] ]-eta[ h_node_ds[bn] ])/dx_h[bn] )
				   *(s_up-s_h_up)+ eta[ h_node_ds[bn] ];
	}

	if (bn == 0)
	{
		us_zeta = zeta[ h_node_us[bn] ];
		us_eta = eta[ h_node_us[bn] ];
	}

	//cout << "LINE 4037 h is: " << us_zeta-us_eta << endl;

	// calculate the depth of the bins in the soil
	double us_soil_depthinterval = (us_zeta-us_eta)/double(n_PDZ_intervals);
	double us_parent_depthinterval = (bottom_depth-(us_zeta-us_eta))/double(n_CAZ_intervals);

  //cout << "LINE 4892, soil thick: " << (us_zeta-us_eta) << " and parent thick: " 
  //    << (bottom_depth-(us_zeta-us_eta)) << " ncaz: " << n_CAZ_intervals << endl;

  //cout << "LINE 4892, getting depth intervals, soil: " << us_soil_depthinterval
  //     << " and parent: " << us_parent_depthinterval << endl;

	// first get the soil nodes
	d_top_locs.push_back(0);
	for (int sdi = 1; sdi<n_PDZ_intervals; sdi++)
	{
		// get the depth index
		d_top_locs.push_back( (double(sdi))*us_soil_depthinterval);
		d_bottom_locs.push_back( (double(sdi))*us_soil_depthinterval);
	}
	
	

	// now the parent material nodes
	for (int pdi = 0; pdi<n_CAZ_intervals; pdi++)
	{
		d_top_locs.push_back( (double(pdi))*us_parent_depthinterval+(us_zeta-us_eta));
		d_bottom_locs.push_back( (double(pdi))*us_parent_depthinterval+(us_zeta-us_eta));
	}

	// and now the bottom node
	d_bottom_locs.push_back(bottom_depth);

	//cout << "h is: " << us_zeta-us_eta << " and n PDZ: " << n_PDZ_intervals << " n CAZ: " << n_CAZ_intervals << endl;
	//int sz = d_bottom_locs.size();
	//for(int i = 0; i<sz; i++)
	//{
	//	cout << "i: " << i << " " << d_top_locs[i] << " " << d_bottom_locs[i] << endl;
	//}


}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#endif
