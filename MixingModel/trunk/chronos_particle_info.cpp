
#include <iostream>
#include <fstream>
#include <math.h>
#include "chronos_particle_info.hpp"
using namespace std;


void Particle_info::create()
 {
  cout << "error, you need to have a filename\n";
 }

// this function loads the parameters from the file with the name fname
void Particle_info::create( const char* fname )
 {
  ifstream in;
  in.open(fname);

  vector<int> t_type_index;
  vector<string> t_type_name;
  vector<double> t_a_vec;
  vector<double> t_alpha_vec;
  vector<double> t_b_vec;
  vector<double> t_beta_vec;
  vector<double> t_rho_vec;
  vector<double> t_rho_c_vec;
  vector<double> t_w_p_vec;
  vector<double> t_w_c_vec;
  vector<double> t_stoic_p_vec;
  vector<double> t_stoic_c_vec;
  vector<double> t_D_vec;
  vector<double> t_SiO2_stoic;
  vector<double> t_Na_stoic;
  vector<double> t_Ca_stoic;
  vector<double> t_K_stoic;
  vector<double> t_Mg_stoic;

  string t_ptype;
  int temp_ti;
  int n_t;
  double t_a,t_alpha,t_b,t_beta,t_rho_p,t_rho_c,t_w_p,t_w_c,t_stoic_p,t_stoic_c,t_D;
  double t_SiO2,t_Na,t_Ca,t_K,t_Mg;

  in >> n_t;

  n_types = n_t;
  for (int i = 1; i<=n_t; i++)
   {
    in >> temp_ti >> t_ptype >> t_a >> t_alpha >> t_b >> t_beta >> t_rho_p >> t_rho_c
       >> t_w_p >> t_w_c >> t_stoic_p >> t_stoic_c >> t_D
       >> t_SiO2 >> t_Na >> t_Ca >> t_K >> t_Mg;
    t_type_index.push_back(temp_ti);
    t_type_name.push_back(t_ptype);
    if(t_a == 0)
     t_a_vec.push_back(0);
    else
     t_a_vec.push_back(pow(10,t_a)*31536000);
    t_alpha_vec.push_back(t_alpha);
    t_b_vec.push_back(t_b);
    t_beta_vec.push_back(t_beta);
    t_rho_vec.push_back(t_rho_p);
    t_rho_c_vec.push_back(t_rho_c);
    t_w_p_vec.push_back(t_w_p);
    t_w_c_vec.push_back(t_w_c);
    t_stoic_p_vec.push_back(t_stoic_p);
    t_stoic_c_vec.push_back( t_stoic_c);
    t_D_vec.push_back(t_D);
    t_SiO2_stoic.push_back(t_SiO2);
    t_Na_stoic.push_back(t_Na);
    t_Ca_stoic.push_back(t_Ca);
    t_K_stoic.push_back(t_K);
    t_Mg_stoic.push_back(t_Mg);
   }
  type_index = t_type_index;
  type_name = t_type_name;
  a_vec = t_a_vec;
  alpha_vec = t_alpha_vec;
  b_vec = t_b_vec;
  beta_vec = t_beta_vec;
  rho_p_vec = t_rho_vec;
  rho_c_vec = t_rho_c_vec;
  w_p_vec = t_w_p_vec;
  w_c_vec = t_w_c_vec;
  stoic_p_vec = t_stoic_p_vec;
  stoic_c_vec = t_stoic_c_vec;
  D_vec = t_D_vec;
  SiO2_stoic = t_SiO2_stoic;
  Na_stoic = t_Na_stoic;
  Ca_stoic = t_Ca_stoic;
  K_stoic = t_K_stoic;
  Mg_stoic = t_Mg_stoic;

  in.close();
 }

// print the parameters of each particle type to file
void Particle_info::details_to_screen()
 {
  cout << "type index: ";
  for (int i = 0; i<n_types; i++)
   cout << type_index[i] << " ";
  cout << endl;
  cout << "type name: ";
  for (int i = 0; i<n_types; i++)
   cout << type_name[i] << " ";
  cout << endl;
  cout << "a: ";
  for (int i = 0; i<n_types; i++)
   cout << a_vec[i] << " ";
  cout << endl;
  cout << "alpha: ";
  for (int i = 0; i<n_types; i++)
   cout << alpha_vec[i] << " ";
  cout << endl;
  cout << "b: ";
  for (int i = 0; i<n_types; i++)
   cout << b_vec[i] << " ";
  cout << endl;
  cout << "beta: ";
  for (int i = 0; i<n_types; i++)
   cout << beta_vec[i] << " ";
  cout << endl;
  cout << "rho_p: ";
  for (int i = 0; i<n_types; i++)
   cout << rho_p_vec[i] << " ";
  cout << endl;
  cout << "rho_c: ";
  for (int i = 0; i<n_types; i++)
   cout << rho_c_vec[i] << " ";
  cout << endl;
  cout << "weight_p: ";
  for (int i = 0; i<n_types; i++)
   cout << w_p_vec[i] << " ";
  cout << endl;
  cout << "weight_c: ";
  for (int i = 0; i<n_types; i++)
   cout << w_c_vec[i] << " ";
  cout << endl;
  cout << "stoic_p: ";
  for (int i = 0; i<n_types; i++)
   cout << stoic_p_vec[i] << " ";
  cout << endl;
  cout << "stoic_c: ";
  for (int i = 0; i<n_types; i++)
   cout << stoic_c_vec[i] << " ";
  cout << endl;
  cout << "D: ";
  for (int i = 0; i<n_types; i++)
   cout << D_vec[i] << " ";
  cout << endl;
  cout << "Si_stoich: ";
  for (int i = 0; i<n_types; i++)
   cout << SiO2_stoic[i] << " ";
  cout << endl;
  cout << "Na_stoich: ";
  for (int i = 0; i<n_types; i++)
   cout << Na_stoic[i] << " ";
  cout << endl;
  cout << "Ca_stoich: ";
  for (int i = 0; i<n_types; i++)
   cout << Ca_stoic[i] << " ";
  cout << endl;
  cout << "K_stoich: ";
  for (int i = 0; i<n_types; i++)
   cout << K_stoic[i] << " ";
  cout << endl;
  cout << "Mg_stoich: ";
  for (int i = 0; i<n_types; i++)
   cout << Mg_stoic[i] << " ";
  cout << endl;
 }


// this function solves the fractional and rate equations through
// time, loading the data into vectors. Then, during the computations,
// the vectors are referenced, thus saving computational time. All the computation
// mass and height fractions is done at the front end within the particle info
// object
// dt is the time spacing, n_timesteps are the number of timesteps in the model run
void Particle_info::generate_time_vecs(double deltat,int nt)
 {
  dt = deltat;
  n_timesteps = nt;

  double exp_multiple; 			// multiples for the
  double r_exp_multiple;
  double t_pow;				// exponential and power laws
  double r_t_pow;
  double chem_frac_ratio;		// ratio of the chemical fraction
  double density_ratio;			// ratio of mineral and clay density
  double mpf;				// fractional mass of primary (temporary)
  double mcf;				// fractional mass of clay (temporary)
  double rf;				// rate fraction
  // each particle type spawns a vector
  vector<double> temp_vec;
  for (int i = 0; i<n_types; i++)
   {
    m_p_frac_vecs.push_back(temp_vec);
    m_c_frac_vecs.push_back(temp_vec);
    h_frac_vecs.push_back(temp_vec);
    m_rate_frac_vecs.push_back(temp_vec);
   }

  // loop through the particle types
  for (int i = 0; i<n_types; i++)
   {
    // get some constants for the given particle type
    // these calculations are used to save computational time
    if (a_vec[i] != 0)
     {
      r_t_pow = alpha_vec[i]+beta_vec[i];
      t_pow = 1+alpha_vec[i]+beta_vec[i];
      r_exp_multiple = -(6*a_vec[i]*b_vec[i]*w_p_vec[i])/(D_vec[i]*rho_p_vec[i]);
      exp_multiple = r_exp_multiple/t_pow;
      chem_frac_ratio = stoic_c_vec[i]*w_c_vec[i]/(stoic_p_vec[i]*w_p_vec[i]);
      density_ratio = rho_p_vec[i]/rho_c_vec[i];
     }

    // loop through time
    for (int tt = 0; tt<n_timesteps; tt++)
     {
      if (tt == 0)
       {
        m_p_frac_vecs[i].push_back(1.0);
        m_c_frac_vecs[i].push_back(0.0);
        h_frac_vecs[i].push_back(1.0);
        m_rate_frac_vecs[i].push_back(0.0);
       }
      else							// !! else 1
       {
        // note that the rate laws for white and brantley are in seconds, but the units of the
        // model are years, so we have to multiply time (in years) by 31536000 sec/year
        if (a_vec[i] != 0)
         {
          mpf = exp( pow(dt*tt,t_pow) * exp_multiple);
          mcf = (1-mpf)*chem_frac_ratio;
          rf  = r_exp_multiple*pow(dt*tt,r_t_pow);
          //cout << "time is: " << tt*dt << " mpf is: " << mpf << " and mcf is: " << mcf << endl;
          m_p_frac_vecs[i].push_back( mpf );
          m_c_frac_vecs[i].push_back(mcf);
          h_frac_vecs[i].push_back(mpf+density_ratio*mcf);
          m_rate_frac_vecs[i].push_back(rf);
         }
        else							// !! else 2
         {
          m_p_frac_vecs[i].push_back( 1.0 );
          m_c_frac_vecs[i].push_back(0.0);
          h_frac_vecs[i].push_back(1.0);
          m_rate_frac_vecs[i].push_back(0.0);
         }							// !! end else 2
       }							// !! end else 1
     }								// !! end time loop
   }								// !! end p-type loop
 }								// !! end function

// this function prints to file all of the fractions and weathering rates
// the first line format is
// n_ts n_types -9999 -9999 ...
// then the format goes:
// time
void Particle_info::print_time_vecs(const char* fname)
 {
  // set up iterateors for accessing the time series data
  vector< vector<double> >::iterator mpf_iter;
  vector< vector<double> >::iterator mcf_iter;
  vector< vector<double> >::iterator hf_iter;
  vector< vector<double> >::iterator rf_iter;
  mpf_iter = m_p_frac_vecs.begin();
  mcf_iter = m_c_frac_vecs.begin();
  hf_iter =  h_frac_vecs.begin();
  rf_iter =  m_rate_frac_vecs.begin();

  ofstream outf;
  outf.open(fname);

  outf << n_timesteps << " " << n_types;
  for (int i = 0; i<(4*n_types)-1; i++)
   outf << " -9999";
  outf << endl;

  for (int i = 0; i< n_timesteps; i++)
   {
    outf << dt*i;
    for (int nt = 0; nt<n_types; nt++)
     outf << " " << (*(mpf_iter+nt))[i] << " " << (*(mcf_iter+nt))[i]
          << " " << (*(hf_iter+nt))[i]  << " " << (*(rf_iter+nt))[i];
    outf << endl;
   }

 }

// this function gets the particle fractions based on mass fractions
vector<double> Particle_info::get_particle_fractions(vector<double> chi_frac)
 {
  int sz = chi_frac.size();
  if (sz != n_types)
   cout << "the mass fractions do not equal the number of types\n";
  double temp_theta;
  double denom;
  vector<double> theta;
  if (sz == 1)
   theta.push_back(1.0);
  else if (sz == 2)
   {
    double r1 = rho_p_vec[0];
    double r2 = rho_p_vec[1];
    double chi1 = chi_frac[0];
    temp_theta = -r2*chi1/(-r1+r1*chi1-r2*chi1);
    theta.push_back(temp_theta);
    temp_theta = -(r1-r1*chi1)/(-r1+r1*chi1-r2*chi1);
    theta.push_back(temp_theta);
   }
  else if (sz == 3)
   {
    double r1 = rho_p_vec[0];
    double r2 = rho_p_vec[1];
    double r3 = rho_p_vec[2];
    double chi1 = chi_frac[0];
    double chi2 = chi_frac[1];
    denom = r1*r2 - r1*r2*chi1 + r2*r3*chi1 - r1*r2*chi2 + r1*r3*chi2;
    temp_theta = r2*r3*chi1/denom;
    theta.push_back(temp_theta);
    temp_theta = r1*r3*chi2/denom;
    theta.push_back(temp_theta);
    temp_theta = (r1*r2 - r1*r2*chi1-r1*r2*chi2)/denom;
    theta.push_back(temp_theta);
   }
  else if (sz == 4)
   {
    double r1 = rho_p_vec[0];
    double r2 = rho_p_vec[1];
    double r3 = rho_p_vec[2];
    double r4 = rho_p_vec[3];
    double chi1 = chi_frac[0];
    double chi2 = chi_frac[1];
    double chi3 = chi_frac[2];
    denom = r1*r2*r3 - r1*r2*r3*chi1 + r2*r3*r4*chi1 - r1*r2*r3*chi2 +
            r1*r3*r4*chi2 - r1*r2*r3*chi3 + r1*r2*r4*chi3;
    temp_theta = r2*r3*r4*chi1/denom;
    theta.push_back(temp_theta);
    temp_theta = r1*r3*r4*chi2/denom;
    theta.push_back(temp_theta);
    temp_theta = r1*r2*r4*chi2/denom;
    theta.push_back(temp_theta);
    temp_theta = (r1*r2*r3      - r1*r2*r3*chi1-
                  r1*r2*r3*chi2 + r1*r2*r3*chi3)/denom;
    theta.push_back(temp_theta);
   }

   return theta;
  }

// this function returns the fraction of the height
double Particle_info::retrieve_h_frac(int Timestep, int type)
 {
  double h_frac;
  //cout << "timestep is: " << Timestep << " and type is: " << type << endl;
  vector< vector<double> >::iterator hf_iter;
  hf_iter = h_frac_vecs.begin();
  h_frac = (*(hf_iter+type-1))[Timestep];
  //cout << "number of types is: " << h_frac_vecs.size() << endl;
  //cout << "number of timesteps is: " << (*(hf_iter+type-1)).size() << endl;
  return h_frac;
 }

// this function returns the fractions and reaction rate for all the particles
vector<double> Particle_info::retrieve_mass_fracs(int Timestep, int type)
 {
  vector<double> m_fracs(3);
  //cout << "timestep is: " << Timestep << " and type is: " << type << endl;
  vector< vector<double> >::iterator mpf_iter;
  vector< vector<double> >::iterator mcf_iter;
  vector< vector<double> >::iterator rf_iter;
  mpf_iter = m_p_frac_vecs.begin();
  mcf_iter = m_c_frac_vecs.begin();
  rf_iter =  m_rate_frac_vecs.begin();
  int nt = type-1;
  //cout << "primary_frac: " << (*(mpf_iter+nt))[Timestep] << endl;

  // get the mass fraction of the primary mineral
  m_fracs[0] = (*(mpf_iter+nt))[Timestep];
  //cout << "primary fraction is: " << m_fracs[0] << endl;

  // get the mass fraction of the clay mineral

  m_fracs[1] = (*(mcf_iter+type-1))[Timestep];
  //cout << "clay fraction is: " << m_fracs[0] << endl;

  // get the rate mass fraction
  m_fracs[2] = (*(rf_iter+type-1))[Timestep];
  //cout << "rate fraction is: " << m_fracs[0] << endl;

  return m_fracs;
 }

vector<double> Particle_info::calculate_m_0_vec(double dx, double h_0, double phi)
 {
  vector<double> m_0_vec(n_types);
  for (int i = 0; i< n_types; i++)
   m_0_vec[i] = dx*h_0*(1-phi)*rho_p_vec[i];
  return m_0_vec;
 }

vector<double> Particle_info::calculate_h_congruent(tParticle& tpart, double deltat, double phi)
 {
  vector<double> return_data(4,0.0);
  int Type = tpart.getType();
  double age = tpart.getAge();
  //double m_0 = tpart.getm_0();
  double m_0 = 0.0;


  //cout << "\n\n\nLINE 419 age: " << age << " and m_0: " << m_0 << endl;
  int TimeStep = int(double(age/deltat));
  //cout << "LINE 421 timestep: " << TimeStep << endl;

  double m_frac;
  double r_frac;

  vector< vector<double> >::iterator mpf_iter;
  vector< vector<double> >::iterator rf_iter;
  mpf_iter = m_p_frac_vecs.begin();
  rf_iter =  m_rate_frac_vecs.begin();


  // get the mass fraction of the primary mineral
  m_frac = (*(mpf_iter+Type))[TimeStep];
  r_frac = (*(rf_iter +Type))[TimeStep];
  //cout << "LINE 433 primary fraction is: " << m_frac << " and mass is: " << m_frac*m_0 << endl;
  double m = m_frac*m_0;

  //cout << "LINE 436 TYPE: " << Type << endl;

  double h = m_frac*m_0/(rho_p_vec[Type]*(1-phi));
  return_data[0] = m_frac;
  return_data[1] = m;
  return_data[2] = r_frac*m;
  return_data[3] = h;

  return return_data;
 }


// this calcualtes weathering on a per particle basis
vector<double> Particle_info::calculate_weathering_on_demand(double t_ime,int type)
{
	vector<double> particle_characteristics(4,0.0);

	double exp_multiple; 			// multiples for the
    double r_exp_multiple;
    double t_pow;				// exponential and power laws
    double r_t_pow;
    double chem_frac_ratio;		// ratio of the chemical fraction
    double density_ratio;		// ratio of mineral and clay density
    double mpf;				// fractional mass of primary (temporary)
    double mcf;				// fractional mass of clay (temporary)
    double rf;				// rate fraction
    double tmf;				// total mass fraction (the fraction of the initial mass
    						// remaining

	int ti = type-1;
    // get some constants for the given particle type
    // these calculations are used to save computational time
    //
    // TO DO DEC 2010
    // PLACE THIS LOOP IN THE OBJECT DEFINITION TO SAVE TIME
    //
    if (alpha_vec[ti] != 0)
    {
		r_t_pow = alpha_vec[ti]+beta_vec[ti];
        t_pow = 1+alpha_vec[ti]+beta_vec[ti];
        r_exp_multiple = -(6*a_vec[ti]*b_vec[ti]*w_p_vec[ti])/(D_vec[ti]*rho_p_vec[ti]);
        exp_multiple = r_exp_multiple/t_pow;
        chem_frac_ratio = stoic_c_vec[ti]*w_c_vec[ti]/(stoic_p_vec[ti]*w_p_vec[ti]);
        density_ratio = rho_p_vec[ti]/rho_c_vec[ti];
    }

    // solve for the time
    if (t_ime <= 0)
    {
		mpf = 1.0;
		mcf = 0.0;
        rf = 0.0;
        tmf = 1.0;
    }
    else							// !! else 1
    {
        // note that the rate laws for white and brantley are in seconds, but the units of the
        // model are years, so we have to multiply time (in years) by 31536000 sec/year
        if (alpha_vec[ti] != 0)
        {
			mpf = exp( pow(t_ime,t_pow) * exp_multiple);
            mcf = (1-mpf)*chem_frac_ratio;
            rf  = r_exp_multiple*pow(t_ime,r_t_pow);
            tmf = mpf*(1-chem_frac_ratio)+chem_frac_ratio;
            //cout << "time is: " << t_ime << " type is: " << type << " alpha is: " << alpha_vec[ti] << " mpf is: " << mpf << " and mcf is: " << mcf << endl;
        }
        else							// !! else 2
        {
			mpf = 1.0;
            mcf = 0.0;
            rf = 0.0;
            tmf = 1.0;
        }							// !! end else 2
    }							// !! end else 1


    particle_characteristics[0] = mpf;
    particle_characteristics[1] = mcf;
    particle_characteristics[2] = rf;
    particle_characteristics[3] = tmf;

    return particle_characteristics;
 }								// !! end function



