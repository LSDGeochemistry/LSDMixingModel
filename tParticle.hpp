//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// tParticle
// A particle object that is used for tracking particle properties
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
#include <vector>
#include <list>
#include "CRN_parameters.hpp"
#include "VolumeParticleInfo.hpp"
using namespace std;

#ifndef tParticle_H
#define tParticle_H

// empty class definition so that we can use friend functions
class CRN_parameters;

/// The particle object. It has CRN and other properties attached to each particle
/// These can be moved between locations during hillslope simulations
class tParticle
{
 public:
    /// @brief default constructor. Has defualt type and age (0) of a particle
    /// @author SMM
    /// @date 01/01/2007
    tParticle()           { create(); }
                
    /// @brief Constructor that simply assigns a type to the particle
    /// @param StartType The type you want the particle to have
    /// @author SMM
    /// @date 01/01/2011
    tParticle( int StartType)        { create(StartType); }

    /// @brief Constructor that simply assigns a type and age to the particle
    /// @param StartType The type you want the particle to have
    /// @param StartAge The starting age of the particle (years). 
    /// @author SMM
    /// @date 01/01/2011
    tParticle( int StartType, double StartAge)   { create(StartType, StartAge); }
    
    /// @brief Constructor that has loads of stuff assigned
    /// @param StartType The type you want the particle to have
    /// @param StartCI the starting cell index
    /// @param StartAge The starting age of the particle (years). 
    /// @param StartOSLage The starting OSL age (years)
    /// @param StartxLoc The starting x location (or s location in a flowtube) (metres)
    /// @param StartdLoc The starting depth of the particle (metres)
    /// @author SMM
    /// @date 01/01/2011
    tParticle( int StartType, int StartCI, double StartAge, double StartOSLage, double StartxLoc, double StartdLoc)
            { create(StartType, StartCI, StartAge, StartOSLage, StartxLoc, StartdLoc); }
            
    /// @brief Constructor that has a type and location
    /// @param StartType The type you want the particle to have
    /// @param StartxLoc The starting x location (or s location in a flowtube) (metres)
    /// @param StartdLoc The starting depth of the particle (metres)
    /// @author SMM
    /// @date 01/01/2011
    tParticle( int StartType, double StartxLoc, double StartdLoc)
            { create(StartType, StartxLoc, StartdLoc); }

    /// @brief This returns the "type", which is an integer key to the type of particle
    ///  the user can define types
    /// @return type and integer key to the type
    int    getType() const			{ return Type; }
    
    /// @brief This returns the "cell index, whioch is sometimes tracked by the particle
    ///  So that it knows where it is. 
    /// @return type and integer key to the type
    int    getCellIndex() const		{ return CellIndex; }
    
    // more getter functions
    double getAge() const			{ return Age; }
    double getOSLage() const		{ return OSLage; }
    double getxLoc() const			{ return xLoc; }
    double getdLoc() const			{ return dLoc; }

    // copy functions
    tParticle(const tParticle& tP)
       { create(tP.getType(),tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc()); }
    tParticle(tParticle& tP)
       { create(tP.getType(),tP.getCellIndex(), tP.getAge(),tP.getOSLage(),tP.getxLoc(),tP.getdLoc()); }

    // operators
    tParticle& operator=(const tParticle& tP);
    std::ostream& operator<< (std::ostream&);

    /// @brief Increase the age of the particle
    /// @param dt The time increment (years)
    /// @author SMM
    /// @date 01/01/2007
    void incrementAge(double dt);
    
    /// @brief Updates the cell index of the particle
    /// @param CI The new cell index
    /// @author SMM
    /// @date 01/01/2007
    void setCellIndex(int CI);
    
    /// @brief Resets OSL age to zero
    /// @author SMM
    /// @date 01/01/2007
    void OSLexpose();
    
    /// @brief Need to check what this does
    /// @author SMM
    /// @date 01/01/2007
    void SoilAgeExpose();

    /// @brief Updates the downslope location
    /// @param new_xLoc The update x location (metres)
    /// @author SMM
    /// @date 01/01/2007
    void update_xLoc(double new_xLoc);
    
    /// @brief This is for combined vertical and horizontal motion
    ///  THe particles have a vertical displacement but the "reflect" off
    ///  the soil-saprolite and soil-surface boundaries
    /// @param dx Change in the x location (metres)
    /// @param dd CHenge in the depth (metres)
    /// @param h soil thickness (metres)
    /// @param dt time interval of displacement (years)
    /// @author SMM
    /// @date 01/01/2007
    void displaceReflect(double dx,double dd,double h,double dt);


    /// @brief this function test to see if the particle is within a hillslope, 
    /// and if not sets the data to no data values
    /// @param lambda the length of the hillslope
    /// @author SMM
    /// @date 01/01/2010   
    int  test_domain(double lambda);

 protected:
 
    /// The type: used to identify, generally, what mineral the particle is
    int Type;
    
    /// The cell index, used for referencing the location of a particle in 
    /// a grid 
    int CellIndex;
    
    /// Age the particle has spent in the soil
    double Age;
    
    /// Optically stimulated luminescence age
    double OSLage;
    
    /// Horizontal location of the particle
    double xLoc;
    
    /// Other coordinate, used in 3D simulations
    double yLoc;
    
    /// Depth of the particle
    double dLoc; 


 private:
 
    /// @brief creates a default particle of type 0, 
    /// cell index of -1, age of 0 and OSLage of -9999
    /// @author SMM
    /// @date 01/01/2008
    void create();
    
    /// @brief creates a default particle 
    /// cell index of -1, age of 0 and OSLage of -9999
    /// @param StartType the starting type of the particle
    /// @author SMM
    /// @date 01/01/2008
    void create(int);
    
    /// @brief creates a default particle of 
    /// cell index of -1, and OSLage of -9999
    /// @param StartType the starting type of the particle
    /// @param StartAge the starting age of the particle
    /// @author SMM
    /// @date 01/01/2008    
    void create(int, double);
    
    /// @brief create function where user assigns all data members excepth
    /// yLoc (yLoc is only used in 3D simulations)
    void create(int, int, double, double, double, double);
    
    /// @brief create function where user assigns type, dLoc and xLoc
    void create(int, double, double);
};

/// the CRN tracer particle object
/// this object tracks the CRN concenentration in a number
/// of particles
class CRN_tParticle: public tParticle
{
  public:
  
    /// @brief default constructor. Has defualt type and age (0) of a particle
    /// @author SMM
    /// @date 01/01/2007
    CRN_tParticle()			{ create(); }
    
  /// @ brief constructor with starting x location, depth and z location
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

    /// @brief The copy constructor
    /// @param tP an LSDCRNParticle object
    /// @author SMM
    /// @date 01/01/2010
    CRN_tParticle& operator=(const CRN_tParticle& tP);
    
    /// @brief The copy constructor
    /// @param tP a constant LSDCRNParticle object
    /// @author SMM
    /// @date 01/01/2010
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

    /// @brief update the 10Be concentration based on a constant erosion rate
    /// using the full production range, including muons. 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr  POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_10Be_conc(double dt,double erosion, CRN_parameters& CRNp);

    /// @brief update the 26Al concentration based on a constant erosion rate
    /// using the full production range, including muons. 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr  POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_26Al_conc(double dt,double erosion, CRN_parameters& CRNp);
    
    /// @brief update the 14C concentration based on a constant erosion rate
    /// using the full production range, including muons. 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr  POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_14C_conc(double dt,double erosion, CRN_parameters& CRNp);

    /// @brief update the 36Cl concentration based on a constant erosion rate
    /// using the full production range, including muons. 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_36Cl_conc(double dt,double erosion, CRN_parameters& CRNp);

    /// @brief update the 21Ne concentration based on a constant erosion rate
    /// using the full production range, including muons. 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_21Ne_conc(double dt,double erosion, CRN_parameters& CRNp);

    /// @brief update the 3He concentration based on a constant erosion rate
    /// using the full production range, including muons. 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_3He_conc(double dt,double erosion, CRN_parameters& CRNp);

    /// @brief This updates 10Be nuclide concentration if erosion is increasing (or decreasing) linearly
    /// alpha is the change in erosion rate--do not set alpha to zero!
    /// @details This solves an analytical solution for cosmo concertration with a linear
    /// increase in cosmo concentration
    /// @param dt time step to calculate next cosmo concetration
    /// @param erosion_rate the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
    /// @param alpha the increase in erosion rate (in g/cm^2/yr^2)
    /// @param CRNp a CRN parameters object
    /// @author SMM
    /// @date 01/01/2010
    void update_10Be_conc_linear_increase(double dt,double erosion_rate, double alpha, CRN_parameters& CRNp);

    /// @brief update the 10Be concentration based on a constant erosion rate
    /// using the ONLY NEUTRON porduction. It is less computationally expensive
    /// than calcualting the full production and is a good approximation of the
    /// erosion rates are slow 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr      POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_10Be_conc_neutron_only(double dt,double erosion, CRN_parameters& CRNp);

    /// @brief update the 26Al concentration based on a constant erosion rate
    /// using the ONLY NEUTRON porduction. It is less computationally expensive
    /// than calcualting the full production and is a good approximation of the
    /// erosion rates are slow 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr       POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_26Al_conc_neutron_only(double dt,double erosion, CRN_parameters& CRNp);

    /// @brief update the 14C concentration based on a constant erosion rate
    /// using the ONLY NEUTRON porduction. It is less computationally expensive
    /// than calcualting the full production and is a good approximation of the
    /// erosion rates are slow 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_14C_conc_neutron_only(double dt,double erosion, CRN_parameters& CRNp);

    /// @brief update the 36Cl concentration based on a constant erosion rate
    /// using the ONLY NEUTRON porduction. It is less computationally expensive
    /// than calcualting the full production and is a good approximation of the
    /// erosion rates are slow 
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr) over a duration of dt
    /// @param dt the timestep over which erosion occurs
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_36Cl_conc_neutron_only(double dt,double erosion, CRN_parameters& CRNp);

    /// @brief This assigns nuclide concentrations with constant values
    /// @param C_10Be the 10Be concentration
    /// @param C_26Al the 26Al concentration
    /// @param C_36Cl the 36Cl concentration
    /// @param C_14C the 14C concentration
    /// @param C_21Ne the 21Ne concentration
    /// @param C_3He the 3He concentration
    void update_cosmo_conc_const(double C_10Be, double C_26Al, double C_36Cl,
                                 double C_14C, double C_21Ne, double C_3He);


    /// @brief Bring the 10Be concentration to steady state based
    /// on a constant erosion rate using full muogenic production.  
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr). It is an analytical
    /// steady state solution
    /// @param erosion the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_10Be_SSfull(double erosion, CRN_parameters& CRNp);
    
    /// @brief Bring the 26Al concentration to steady state based
    /// on a constant erosion rate using full muogenic production.  
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr). It is an analytical
    /// steady state solution
    /// @param erosion the erosion rate in g/cm^2/yr  POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_26Al_SSfull(double erosion, CRN_parameters& CRNp);
    
    /// @brief Bring the 14C concentration to steady state based
    /// on a constant erosion rate using full muogenic production.  
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr). It is an analytical
    /// steady state solution
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_14C_SSfull(double erosion, CRN_parameters& CRNp);
    
    /// @brief Bring the 36Cl concentration to steady state based
    /// on a constant erosion rate using full muogenic production.  
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr). It is an analytical
    /// steady state solution
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_36Cl_SSfull(double erosion, CRN_parameters& CRNp);
    
    /// @brief Bring the 21Ne concentration to steady state based
    /// on a constant erosion rate using full muogenic production.  
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr). It is an analytical
    /// steady state solution
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_21Ne_SSfull(double erosion, CRN_parameters& CRNp);

    /// @brief Bring the 3He concentration to steady state based
    /// on a constant erosion rate using full muogenic production.  
    /// @details This function solves for the updated concentration assuming
    /// a constant erosion rate (in g/cm^2/yr). It is an analytical
    /// steady state solution
    /// @param erosion the erosion rate in g/cm^2/yr   POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// to approximate production from the different production mechanisms
    /// @author SMM
    /// @date 01/01/2010
    void update_3He_SSfull(double erosion, CRN_parameters& CRNp);

    /// @brief A wrapper function to update all the nuclide concentrations
    /// in one go. It uses full muogenic production
    /// @param dt the timestep over which the concetrations are updated
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// @author SMM
    /// @date 01/01/2014
    void update_all_CRN(double dt, double erosion, CRN_parameters& CRNp);
  
    /// @brief A wrapper function to update all the nuclide concentrations
    /// in one go. It uses NEUTRON PRODUCTION ONLY to save computational
    /// expense. This is a reasonable approximation in slowly eroding landscapes
    /// @param dt the timestep over which the concetrations are updated
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// @author SMM
    /// @date 01/01/2014
    void update_all_CRN_neutron_only(double dt, double erosion, CRN_parameters& CRNp);
    
    /// @brief A wrapper function to update all the nuclide concentrations
    /// to steady state using full muon production
    /// @param erosion the erosion rate in g/cm^2/yr    POSITIVE FOR EROSION
    /// @param CRNp a CRN parameters object that stores the coefficients
    /// @author SMM
    /// @date 01/01/2014
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
