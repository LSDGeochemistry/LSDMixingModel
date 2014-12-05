//tParticle.h
// header file for the discrete particle object
// this version of the dicrete particle keeps track
// of two integers, the time and the position

#include <iostream>
#include <vector>
#include <list>
#include "CRN_parameters.hpp"
#include "VolumeParticleInfo.hpp"
using namespace std;

#ifndef tParticle_H
#define tParticle_H

// empty class definition so that we can use friend functions
class CRN_parameters;

class tParticle
{
 public:
    tParticle()					{ create(); }
    tParticle( int StartType)			{ create(StartType); }
    tParticle( int StartType, double StartAge)	{ create(StartType, StartAge); }
    tParticle( int StartType, int StartCI, double StartAge, double StartOSLage, double StartxLoc, double StartdLoc)
    						{ create(StartType, StartCI, StartAge, StartOSLage, StartxLoc, StartdLoc); }
    tParticle( int StartType, double StartxLoc, double StartdLoc)
    						{ create(StartType, StartxLoc, StartdLoc); }

    int    getType() const			{ return Type; }
    int    getCellIndex() const		{ return CellIndex; }
    double getAge() const			{ return Age; }
    double getOSLage() const		{ return OSLage; }
    double getxLoc() const			{ return xLoc; }
    double getdLoc() const			{ return dLoc; }

    tParticle(const tParticle& tP)
    	{ create(tP.getType(),tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc()); }
    tParticle(tParticle& tP)
    	{ create(tP.getType(),tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc()); }

    tParticle& operator=(const tParticle& tP);

    std::ostream& operator<< (std::ostream&);

    void incrementAge(double dt);
    void setCellIndex(int CI);
    void OSLexpose();
    void SoilAgeExpose();
	void update_xLoc(double new_xLoc);
    void displaceReflect(double dx,double dd,double h,double dt);
    int  test_domain(double lambda);

 protected:
    int Type;				// note for volumetric calcualtions this
    						// is an INDEX into the vector that holds the
    						// type names.
    int CellIndex;
    double Age;
    double OSLage;
    double xLoc;
    double dLoc;

 private:
    void create();
    void create(int);
    void create(int, double);
    void create(int, int, double, double, double, double);
    void create(int, double, double);
};

// the CRN tracer particle object
// this object tracks the CRN concnentration in a number
// of particles
class CRN_tParticle: public tParticle
{
	public:
	CRN_tParticle()			{ create(); }
	CRN_tParticle(int startType, double startxLoc,double startdLoc,
					double start_effdloc, double startzloc)
							{ create(startType,startxLoc,startdLoc,
							  start_effdloc,startzloc); }
	CRN_tParticle(int startType, double startxLoc,double startzeta_Loc)
							{ create(startType,startxLoc,startzeta_Loc); }
	CRN_tParticle(int startType, int startGSDType, double startxLoc,
		            double startdLoc, double start_effdloc,
		            double startzLoc, double startMass,
		            double startSurfaceArea)
							{ create(startType, startGSDType, startxLoc, startdLoc, start_effdloc,
		            startzLoc, startMass, startSurfaceArea); }
	CRN_tParticle(int startType, double startxLoc,double startzeta_Loc,
					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C)
							{ create(startType,startxLoc,startzeta_Loc,
								startdLoc, start_effdLoc,
								start_C10Be, start_C14C); }
	CRN_tParticle(int startType, double startxLoc,double startzeta_Loc,
					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C, double start_21Ne)
							{ create(startType,startxLoc,startzeta_Loc,
								startdLoc, start_effdLoc,
								start_C10Be, start_C14C, start_21Ne); }
	CRN_tParticle(int startType, int startCellIndex, double startAge, double startOSLAge,
					double startxLoc,double startdLoc, double startefdLoc,
					double startzLoc, double start_C10Be, double start_C26Al,
					double start_C36Cl, double start_C14C,
					double start_C21Ne, double start_C3He,
					double start_Cf7Be, double start_Cf10Be,
					double start_Cf210Pb, double start_Cf137Cs,
					double start_Mass, double start_StartingMass,
					double start_SurfaceArea,
					int start_GSDType)
							{ create(startType, startCellIndex, startAge, startOSLAge,
								startxLoc, startdLoc, startefdLoc,
								startzLoc, start_C10Be, start_C26Al,
								start_C36Cl, start_C14C,
								start_C21Ne, start_C3He,
								start_Cf7Be, start_Cf10Be,
								start_Cf210Pb, start_Cf137Cs,
								start_Mass, start_StartingMass,
								start_SurfaceArea,
								start_GSDType); }


    CRN_tParticle(const CRN_tParticle& tP)
    	{ create(tP.getType(), tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc(),
    	         tP.geteffective_dLoc(), tP.get_zetaLoc(),tP.getConc_10Be(),
    	         tP.getConc_26Al(), tP.getConc_36Cl(), tP.getConc_14C(),
    	         tP.getConc_21Ne(), tP.getConc_3He(),
    	         tP.getConc_f7Be(), tP.getConc_f10Be(),
    	         tP.getConc_f210Pb(), tP.getConc_f137Cs(),
    	         tP.getMass(), tP.getStartingMass(),
				 tP.getSurfaceArea(),
    	         tP.getGSDType()); }
    CRN_tParticle(CRN_tParticle& tP)
    	{ create(tP.getType(), tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc(),
    	         tP.geteffective_dLoc(),tP.get_zetaLoc(),tP.getConc_10Be(),
    	         tP.getConc_26Al(), tP.getConc_36Cl(), tP.getConc_14C(),
    	         tP.getConc_21Ne(), tP.getConc_3He(),
    	         tP.getConc_f7Be(), tP.getConc_f10Be(),
    	         tP.getConc_f210Pb(), tP.getConc_f137Cs(),
    	         tP.getMass(), tP.getStartingMass(),
				 tP.getSurfaceArea(),
    	         tP.getGSDType());   }

	CRN_tParticle& operator=(const CRN_tParticle& tP);
	CRN_tParticle& operator=(CRN_tParticle& tP);



	// accessing function to get at the data elements
	double getConc_10Be() const				{ return Conc_10Be; }
	double getConc_26Al() const				{ return Conc_26Al; }
	double getConc_36Cl() const				{ return Conc_36Cl; }
	double getConc_14C() const				{ return Conc_14C; }
	double getConc_21Ne() const				{ return Conc_21Ne; }
	double getConc_3He() const				{ return Conc_3He; }
	double getConc_f7Be() const				{ return Conc_f7Be; }
	double getConc_f10Be() const			{ return Conc_f10Be; }
	double getConc_f210Pb() const			{ return Conc_f210Pb; }
	double getConc_f137Cs() const			{ return Conc_f137Cs; }
	double geteffective_dLoc() const		{ return effective_dLoc; }
	double get_zetaLoc() const				{ return zetaLoc; }

	// accessing data elements for volumetric data
	double getMass() const  				{ return Mass; }
	double getStartingMass() const			{ return StartingMass; }
	double getSurfaceArea()	const			{ return SurfaceArea; }
	int getGSDType()	const				{ return GSDType; }

	// update nuclide concentrations in the event of constant erosion
	// erosion is in g/cm^2/yr
	void update_10Be_conc(double dt,double erosion, CRN_parameters& CRNp);
	void update_26Al_conc(double dt,double erosion, CRN_parameters& CRNp);
	void update_14C_conc(double dt,double erosion, CRN_parameters& CRNp);
	void update_36Cl_conc(double dt,double erosion, CRN_parameters& CRNp);
	void update_21Ne_conc(double dt,double erosion, CRN_parameters& CRNp);
	void update_3He_conc(double dt,double erosion, CRN_parameters& CRNp);

	// update nuclide concentrations if erosion is increasing (or decreasing) linearly
	// alpha is the change in erosion rate--do not set alpha to zero
	void update_10Be_conc_linear_increase(double dt,double erosion_rate, double alpha, CRN_parameters& CRNp);

	// update nuclide concentrations using only spallation
	void update_10Be_conc_neutron_only(double dt,double erosion, CRN_parameters& CRNp);
	void update_26Al_conc_neutron_only(double dt,double erosion, CRN_parameters& CRNp);
	void update_14C_conc_neutron_only(double dt,double erosion, CRN_parameters& CRNp);
	void update_36Cl_conc_neutron_only(double dt,double erosion, CRN_parameters& CRNp);

	// update nuclide concentrations with constant values
	void update_cosmo_conc_const(double C_10Be, double C_26Al, double C_36Cl,
								 double C_14C, double C_21Ne, double C_3He);

	// update the concentration of nuclides if there is a steady state
	// profile
	void update_10Be_SSfull(double erosion, CRN_parameters& CRNp);
	void update_26Al_SSfull(double erosion, CRN_parameters& CRNp);
	void update_14C_SSfull(double erosion, CRN_parameters& CRNp);
	void update_36Cl_SSfull(double erosion, CRN_parameters& CRNp);
	void update_21Ne_SSfull(double erosion, CRN_parameters& CRNp);
	void update_3He_SSfull(double erosion, CRN_parameters& CRNp);

	// functions collecting the updating functions: these update
	// all nuclides at once
	void update_all_CRN(double dt, double erosion, CRN_parameters& CRNp);
	void update_all_CRN_neutron_only(double dt, double erosion, CRN_parameters& CRNp);
	void update_all_CRN_SSfull(double erosion_rate, CRN_parameters& CRNp);

	// caluclate the apparent erosion from nuclide concentrations
	double apparent_erosion_10Be_neutron_only(double rho, CRN_parameters& CRNp);
	double apparent_erosion_26Al_neutron_only(double rho, CRN_parameters& CRNp);
	double apparent_erosion_14C_neutron_only(double rho, CRN_parameters& CRNp);
	double apparent_erosion_36Cl_neutron_only(double rho, CRN_parameters& CRNp);
	double apparent_erosion_21Ne(double rho, CRN_parameters& CRNp);
	double apparent_erosion_3He(double rho, CRN_parameters& CRNp);

	// functions for dealing with fallout numclides
	void update_fallout10Be_simple_density(double dt, double M_supply_surface,
					double rho_s, double k_f10Be, double deltad, CRN_parameters& CRNp);
	void update_fallout10Be_simple_density_2exp(double dt, double M_supply_surface,
					double rho_skg, double k1_f10Be, double k2_f10Be, double chi_f10Be,
					double deltad_m, CRN_parameters& CRNp);

	// functions for managing the shielding depth, depth, and elevation of
	// particles
	void update_depths(double delta_d, double delta_ed);
	void update_zetaLoc(double new_zeta);
	void update_zetaLoc_with_new_surface(double new_zeta);
	void erode_mass_only(double dt, double mass_erosion_rate);
	void erode_mass_only_linear_increase(double dt, double mass_erosion_rate, double alpha);

	// fuctions for dealing with volumes and surface areas
	double update_surface_area_and_get_volume(VolumeParticleInfo vpi);
	double weather_particle(VolumeParticleInfo vpi, list< vector<double> >& loss_per_surface_area);

    protected:
    double effective_dLoc;		// the effective depth: in g/cm^2
	double zetaLoc;				// the elevation relative to arbitrary
								// datum (m)
	double Conc_10Be;			// concnetration of 10Be in atoms/g
	double Conc_26Al;			// a/g
	double Conc_36Cl;			// a/g
	double Conc_14C;			// a/g
	double Conc_21Ne;			// a/g
	double Conc_3He;			// a/g
	double Conc_f7Be;			// fallout units tba
	double Conc_f10Be; 			// fallout, units a/g
	double Conc_f210Pb;			// fallout, units tba
	double Conc_f137Cs;			// fallout, units tba

	double Mass;				// in kg
	double StartingMass;		// in kg
	double SurfaceArea;
								// in m^2
	int GSDType;				// an integer used to denote the GSD type
								// this is an INDEX into a vector in the VolumeParticleInfo
								// object


	private:
	// functions for creating particles
	void create();
	void create(int startType, double startxLoc,
	            double startzLoc);
	void create(int startType, double startxLoc,
	            double startdLoc, double start_effdloc,
	            double startzloc);
	void create(int startType, int startGSDType,
				double startxLoc,
	            double startdLoc, double start_effdloc,
	            double startzLoc, double startMass,
	            double startSurfaceArea);
	void create(int startType, double startxLoc,double startzeta_Loc,
					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C);
	void create(int startType, double startxLoc,double startzeta_Loc,
					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C, double start_21Ne);
	void create(int startType, int startCellIndex, double startAge, double startOSLAge,
	            double startxLoc,double startdLoc, double startefdLoc,
				double startzLoc, double start_C10Be, double start_C26Al,
				double start_C36Cl, double start_C14C,
				double start_C21Ne, double start_C3He,
				double start_Cf7Be, double start_Cf10Be,
				double start_Cf210Pb, double start_Cf137Cs,
				double start_Mass, double start_StartingMass,
				double start_SurfaceArea,
				int start_GSDType);
	void create(int startType, int startCellIndex, double startAge, double startOSLAge,
	            double startxLoc,double startdLoc, double startefdLoc,
				double startzLoc, double start_C10Be, double start_C26Al,
				double start_C36Cl, double start_C14C,
				double start_C21Ne, double start_C3He,
				double start_Cf7Be, double start_Cf10Be,
				double start_Cf210Pb, double start_Cf137Cs);

};


#endif