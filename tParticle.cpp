//tParticle.cpp
// implementation of the dParticle object
#include <fstream>
#include <math.h>
#include <iostream>
#include <list>
#include "mathutil.hpp"
#include "tParticle.hpp"
#include "VolumeParticleInfo.hpp"
using namespace std;

const double one_min_exp_neg_2 = 1-exp(-2);
const double one_min_exp_neg_5 = 1-exp(-5);

/**************************************************************************\
**
**  tParticle::tParticle:  Constructor for particles.
**			   The intial age of the particle is zero
**			   The location is assigned
\**************************************************************************/
void tParticle::create()
 {
  Type = 0;
  CellIndex = -1;
  Age = 0;
  OSLage = -9999;
 }


void tParticle::create( int StartType )
 {
  Type = StartType;
  CellIndex = -1;
  Age = 0;
  OSLage = -9999;
 }

void tParticle::create( int StartType, double StartAge )
 {
  Type = StartType;
  CellIndex = -1;
  Age = StartAge;
 }

void tParticle::create( int StartType, int StartCI,  double StartAge, double StartOSLage, double StartxLoc, double StartdLoc)
 {
  Type = StartType;
  CellIndex = StartCI;
  Age = StartAge;
  OSLage = StartOSLage;
  xLoc = StartxLoc;
  dLoc = StartdLoc;
 }

void tParticle::create(int StartType, double StartxLoc, double StartdLoc)
 {
  Type = StartType;
  CellIndex = -1;
  Age = 0;
  OSLage = -9999;
  xLoc = StartxLoc;
  dLoc = StartdLoc;
 }

tParticle& tParticle::operator=(const tParticle& rhs)
 {
  if (&rhs != this)
   {
    create(rhs.getType(),rhs.getCellIndex(),rhs.getAge(),rhs.getOSLage(), rhs.getxLoc(),rhs.getdLoc());
   }
  return *this;
 }

std::ostream& operator<<(std::ostream& os, const tParticle& tP)
 {
  os <<   tP.getType() << " " << tP.getCellIndex() << " " << tP.getAge() << " "
       << tP.getOSLage() << " " << tP.getxLoc() << " " << tP.getdLoc();
  return os;
 }

void tParticle::incrementAge(double dt)
 {
  if (Age < 0)
   Age = dt;
  else
   Age += dt;
  if (OSLage > 0)
   OSLage += dt;
 }

void tParticle::setCellIndex(int tempCI)
{
	CellIndex = tempCI;
}

void tParticle::OSLexpose()
 {
	 OSLage = 0;
 }

 void tParticle::SoilAgeExpose()
  {
 	 Age = 0;
 }

void tParticle::update_xLoc(double new_xLoc)
{
	xLoc = new_xLoc;
}

void tParticle::displaceReflect(double dx,double dd,double h, double dt)
 {
  xLoc += dx;

  //std::cout << "tPart.cpp LINE 77 dx  is: " << dd << std::endl;
  double dOld = dLoc;
  double dNew = dLoc+dd;
  if (dNew > h)
   {
    //std::cout << "tPart.cpp LINE 84 dNew  is: " << dNew << " and dd is: " << dd
    //          << " and dOld is: "<< dOld << std::endl;
    dLoc = 2*h - dd - dOld;
    //std::cout << "tPart.cpp LINE 86 dLoc  is: " << dLoc << std::endl;
   }
  else if (dNew <= 0)
   {
    dLoc = -(dd+dOld);
    OSLage = (dd+dOld)*dt/dd;
   }
  else
   dLoc = dNew;

  //if (dLoc > h)
  // std::cout << "tPart.cpp LINE 86 dLoc  is: " << dLoc << std::endl;
 }

int  tParticle::test_domain(double lambda)
 {
  int td;
  if (xLoc >= 0 && xLoc <= lambda)
   td = 1;
  else
   {
    td = -1;
    OSLage = -9999;
    Age = -9999;
    xLoc = -9999;
    dLoc = -9999;
   }
  return td;
 }


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// CRN_tParticle object
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// functions for the CRN loaded tracer particle
void CRN_tParticle::create()
 {
	double rho_r = 2650;		// in kg/m^3

	Type = 0;
	CellIndex = -1;
	Age = 0;
	OSLage = 0;
	xLoc = 0;			// in meters
	dLoc = 5;			// in meters
	zetaLoc = 100;
    Conc_10Be = 0.0;
	Conc_26Al = 0.0;
	Conc_36Cl = 0.0;
	Conc_14C = 0.0;
	Conc_21Ne = 0.0;
	Conc_3He = 0.0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
	effective_dLoc = rho_r*dLoc*0.1;	// the 0.1
	 									// converts between kg/m^2
	 									// to g/cm^2
 }

void CRN_tParticle::create(int startType, double startxLoc,
	            double startzLoc)
{
	double rho_r = 2650;		// in kg/m^3

	Type = startType;
	CellIndex = -1;
	Age = 0;
	OSLage = 0;
	xLoc = startxLoc;			// in meters
	dLoc = 0;			// in meters
	zetaLoc = startzLoc;
    Conc_10Be = 0.0;
	Conc_26Al = 0.0;
	Conc_36Cl = 0.0;
	Conc_14C = 0.0;
	Conc_21Ne = 0.0;
	Conc_3He = 0.0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
	effective_dLoc = rho_r*dLoc*0.1;	// the 0.1
	 									// converts between kg/m^2
	 									// to g/cm^2
 }

void CRN_tParticle::create(int startType, double startxLoc,
	            double startdLoc, double start_effdloc,
	            double startzLoc)
{
	Type = startType;
	CellIndex = -1;
	Age = -99;
	OSLage = -99;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	zetaLoc = startzLoc;
    Conc_10Be = 0.0;
	Conc_26Al = 0.0;
	Conc_36Cl = 0.0;
	Conc_14C = 0.0;
	Conc_21Ne = 0.0;
	Conc_3He = 0.0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
	effective_dLoc = start_effdloc;
}

void CRN_tParticle::create(int startType, int startGSDType,
				double startxLoc,
	            double startdLoc, double start_effdloc,
	            double startzLoc, double startMass,
	            double startSurfaceArea)
{
	Mass = startMass;					// in kg
	StartingMass = startMass;			// in kg
	SurfaceArea = startSurfaceArea;
										// in m^2/kg

	Type = startType;
	GSDType = startGSDType;
	CellIndex = -1;
	Age = -99;
	OSLage = -99;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	zetaLoc = startzLoc;
	Conc_26Al = 0.0;
	Conc_36Cl = 0.0;
	Conc_14C = 0.0;
	Conc_21Ne = 0.0;
	Conc_3He = 0.0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
	effective_dLoc = start_effdloc;
}

void CRN_tParticle::create(int startType, double startxLoc,double startzeta_Loc,
 					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C)
{

	Type = startType;
	CellIndex = -1;
	Age = -99;
	OSLage = -99;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	effective_dLoc = start_effdLoc;
	zetaLoc = startzeta_Loc;
    Conc_10Be = start_C10Be;
	Conc_26Al = 0;
	Conc_36Cl = 0;
	Conc_14C = start_C14C;
	Conc_21Ne = 0;
	Conc_3He = 0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
}

void CRN_tParticle::create(int startType, double startxLoc,double startzeta_Loc,
 					double startdLoc, double start_effdLoc,
					double start_C10Be,double start_C14C, double start_21Ne)
{
	Type = startType;
	CellIndex = -1;
	Age = -99;
	OSLage = -99;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	effective_dLoc = start_effdLoc;
	zetaLoc = startzeta_Loc;
    Conc_10Be = start_C10Be;
	Conc_26Al = 0;
	Conc_36Cl = 0;
	Conc_14C = start_C14C;
	Conc_21Ne = start_21Ne;
	Conc_3He = 0;
	Conc_f7Be = 0.0;
	Conc_f10Be = 0.0;
	Conc_f210Pb = 0.0;
	Conc_f137Cs = 0.0;
}

void CRN_tParticle::create(int startType, int startCellIndex, double startAge, double startOSLAge,
	            double startxLoc,double startdLoc, double startefdLoc,
				double startzLoc, double start_C10Be, double start_C26Al,
				double start_C36Cl, double start_C14C,
				double start_C21Ne, double start_C3He,
				double start_Cf7Be, double start_Cf10Be,
				double start_Cf210Pb, double start_Cf137Cs)
{
	Type = startType;
	CellIndex = startCellIndex;
	Age = startAge;
	OSLage = startOSLAge;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	effective_dLoc = startefdLoc;
	zetaLoc = startzLoc;
    Conc_10Be = start_C10Be;
	Conc_26Al = start_C26Al;
	Conc_36Cl = start_C36Cl;
	Conc_14C = start_C14C;
	Conc_21Ne = start_C21Ne;
	Conc_3He = start_C3He;
	Conc_f7Be = start_Cf7Be;
	Conc_f10Be = start_Cf10Be;
	Conc_f210Pb = start_Cf210Pb;
	Conc_f137Cs = start_Cf137Cs;
}

void CRN_tParticle::create(int startType, int startCellIndex, double startAge, double startOSLAge,
	            double startxLoc,double startdLoc, double startefdLoc,
				double startzLoc, double start_C10Be, double start_C26Al,
				double start_C36Cl, double start_C14C,
				double start_C21Ne, double start_C3He,
				double start_Cf7Be, double start_Cf10Be,
				double start_Cf210Pb, double start_Cf137Cs,
				double start_Mass, double start_StartingMass,
				double start_SurfaceArea,
				int start_GSDType)
{
	Type = startType;
	CellIndex = startCellIndex;
	Age = startAge;
	OSLage = startOSLAge;
	xLoc = startxLoc;			// in meters
	dLoc = startdLoc;			// in meters
	effective_dLoc = startefdLoc;
	zetaLoc = startzLoc;
    Conc_10Be = start_C10Be;
	Conc_26Al = start_C26Al;
	Conc_36Cl = start_C36Cl;
	Conc_14C = start_C14C;
	Conc_21Ne = start_C21Ne;
	Conc_3He = start_C3He;
	Conc_f7Be = start_Cf7Be;
	Conc_f10Be = start_Cf10Be;
	Conc_f210Pb = start_Cf210Pb;
	Conc_f137Cs = start_Cf137Cs;
	Mass = start_Mass;
	StartingMass = start_StartingMass;
	SurfaceArea =  start_SurfaceArea;
	GSDType = start_GSDType;
}


CRN_tParticle& CRN_tParticle::operator=(const CRN_tParticle& rhs)
{
	if (&rhs != this)
    {
		create(rhs.getType(), rhs.getCellIndex(), rhs.getAge(),rhs.getOSLage(),rhs.getxLoc(),rhs.getdLoc(),
    	         rhs.geteffective_dLoc(),rhs.get_zetaLoc(),rhs.getConc_10Be(),
    	         rhs.getConc_26Al(), rhs.getConc_36Cl(), rhs.getConc_14C(),
    	         rhs.getConc_21Ne(), rhs.getConc_3He(),
    	         rhs.getConc_f7Be(), rhs.getConc_f10Be(),
    	         rhs.getConc_f210Pb(), rhs.getConc_f137Cs(),
    	         rhs.getMass(), rhs.getStartingMass(),
    	         rhs.getSurfaceArea(),
    	         rhs.getGSDType());
    }
    return *this;
}

CRN_tParticle& CRN_tParticle::operator=(CRN_tParticle& rhs)
{
	if (&rhs != this)
    {
		create(rhs.getType(), rhs.getCellIndex(), rhs.getAge(),rhs.getOSLage(),rhs.getxLoc(),rhs.getdLoc(),
     	         rhs.geteffective_dLoc(),rhs.get_zetaLoc(),rhs.getConc_10Be(),
     	         rhs.getConc_26Al(), rhs.getConc_36Cl(), rhs.getConc_14C(),
    	         rhs.getConc_21Ne(), rhs.getConc_3He(),
    	         rhs.getConc_f7Be(), rhs.getConc_f10Be(),
    	         rhs.getConc_f210Pb(), rhs.getConc_f137Cs(),
    	         rhs.getMass(), rhs.getStartingMass(),
    	         rhs.getSurfaceArea(),
    	         rhs.getGSDType());
    }
    return *this;
}

// this function just sets the concentration of cosmogenic in situ nuclides
void CRN_tParticle::update_cosmo_conc_const(double C_10Be, double C_f10Be, double C_26Al, double C_36Cl,
								 double C_14C, double C_21Ne, double C_3He)
{
	Conc_10Be = C_10Be;
    Conc_f10Be = C_f10Be;
	Conc_26Al = C_26Al;
	Conc_36Cl = C_36Cl;
	Conc_14C  = C_14C;
	Conc_21Ne = C_21Ne;
	Conc_3He  = C_3He;
}


// this function updates the concntration of CRNs in a particle
// the model assumes that during the timestep the change in the
// 'depth' of the particle occurs ofver a constant rate.
// The depth in this case is an equvalent depth...it is linearly
// proportional to depth if overlying density is constant, but
// really represetnt the mass per area above a point in the subsurface
// thus the erosion rate represent a mass removal rate.
void CRN_tParticle::update_10Be_conc(double dt,double erosion_rate, CRN_parameters& CRNp)
{
	// fist berillium
	double Be_exp = exp(-dt*CRNp.lambda_10Be);
	//double Be_depth_exp = exp(-d_0/Gamma_Be1);
    //cout << "F_10Be =" << CRNp.F_10Be[0] << endl;
    //cout << "LINE 236 tParticle.cpp " << endl;
    //cout << lambda_10Be << " " << Be_exp << endl;
    //cout << "starting conc: " << Conc_10Be << endl;
    //cout << "starting depth: " << effective_dLoc << endl;
    //cout << "erosion_rate: "<< erosion_rate << endl;
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_10Be[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
		           (exp(dt*erosion_rate/CRNp.Gamma[i])-
		            exp(-dt*CRNp.lambda_10Be))/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_10Be);
		//cout << "F_10Be["<<i<<"]: " << CRNp.F_10Be[i] << " and sum_term: " << sum_term <<endl;
	}

    //cout << "and sum term is: " << sum_term << endl;

	Conc_10Be = Conc_10Be*Be_exp +  CRNp.S_t*CRNp.P0_10Be*sum_term;
	//cout << "and ending 10Be conc is: " << Conc_10Be <<endl;
}

// this function updates the 10Be concentration if there is a linear
// increase (or decrease) in erosion rate.
void CRN_tParticle::update_10Be_conc_linear_increase(double dt,double erosion_rate, double alpha, CRN_parameters& CRNp)
{
	// fist berillium
	double Be_exp = exp(-dt*CRNp.lambda_10Be);
	//double Be_depth_exp = exp(-d_0/Gamma_Be1);

    //cout << "LINE 236 tParticle.cpp " << endl;
    //cout << lambda_10Be << " " << Be_exp << endl;
    //cout << "starting conc: " << Conc_10Be << endl;
    //cout << "starting depth: " << effective_dLoc << endl;
    //cout << "erosion_rate: "<< erosion_rate << endl;
	double sum_term = 0;
	double L_j, M_j, A_j, B_j, D_j;
	double erfi_arg1, erfi_arg2;
	for (int j = 0; j<4; j++)
	{
		A_j = sqrt(0.5*CRNp.Gamma[j]/alpha);
		B_j = sqrt(2*alpha*CRNp.Gamma[j]);
		D_j = sqrt(0.5*alpha/CRNp.Gamma[j]);
		L_j = exp(-( (erosion_rate+CRNp.Gamma[j]*CRNp.lambda_10Be)*(erosion_rate+CRNp.Gamma[j]*CRNp.lambda_10Be)+
					2*alpha*(effective_dLoc+dt*CRNp.Gamma[j]*CRNp.lambda_10Be) )/(B_j*B_j) );
		erfi_arg1 = erosion_rate/B_j + A_j*CRNp.lambda_10Be;
		erfi_arg2 = D_j*dt + erfi_arg1;
		M_j = erfi(erfi_arg2) - erfi(erfi_arg1);

		//std::cout << "j is: " << j << endl;
		//std::cout << "erfi_arg1 is: " << erfi_arg1 << endl;
		//std::cout << "erfi_arg2 is: " << erfi_arg2 << endl;
		//std::cout << "M_j term is: " << M_j << endl;
		//std::cout << "L_j term is: " << L_j << endl;


		sum_term+= (CRNp.F_10Be[j]*A_j*L_j*M_j);
	}

    //cout << "and sum term is: " << sum_term << endl;

	Conc_10Be = Conc_10Be*Be_exp +  CRNp.S_t*CRNp.P0_10Be*sum_term*sqrt(M_PI);
	//cout << "and ending 10Be conc is: " << Conc_10Be <<endl;
}

void CRN_tParticle::update_10Be_conc_neutron_only(double dt,double erosion_rate, CRN_parameters& CRNp)
{

	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double Be_exp = exp(-dt*CRNp.lambda_10Be);

	double sum_term = 0;
	sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
		           (exp(dt*erosion_rate/Gamma_neutron)-
		            exp(-dt*CRNp.lambda_10Be))/
		           (erosion_rate+Gamma_neutron*CRNp.lambda_10Be);

	Conc_10Be = Conc_10Be*Be_exp +  CRNp.S_t*Be_exp*CRNp.P0_10Be*sum_term;
}

void CRN_tParticle::update_10Be_SSfull(double erosion_rate, CRN_parameters& CRNp)
{
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_10Be[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_10Be);
	}

    //cout << "and sum term is: " << sum_term << endl;

	Conc_10Be = CRNp.S_t*CRNp.P0_10Be*sum_term;
	//cout << "and ending 10Be conc is: " << Conc_10Be <<endl;
}

double CRN_tParticle::apparent_erosion_10Be_neutron_only(double rho, CRN_parameters& CRNp)
{

	double rho_r = rho/1000;
	// a few constants, all computed from Vermeesh 2007
	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double app_eros;

	app_eros = Gamma_neutron*(CRNp.S_t*CRNp.P0_10Be/Conc_10Be-CRNp.lambda_10Be)/rho_r;
	return -app_eros/100.0;
}

void CRN_tParticle::update_26Al_conc(double dt,double erosion_rate, CRN_parameters& CRNp)
{
	double Al_exp = exp(-dt*CRNp.lambda_26Al);

	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_26Al[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
		           (exp(dt*erosion_rate/CRNp.Gamma[i])-
		            exp(-dt*CRNp.lambda_26Al))/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_26Al);
	}

	Conc_26Al = Conc_26Al*Al_exp +  CRNp.S_t*Al_exp*CRNp.P0_26Al*sum_term;
}

void CRN_tParticle::update_26Al_conc_neutron_only(double dt,double erosion_rate, CRN_parameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2
	double Al_exp = exp(-dt*CRNp.lambda_26Al);

	double sum_term = 0;
	sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
		           (exp(dt*erosion_rate/Gamma_neutron)-
		            exp(-dt*CRNp.lambda_26Al))/
		           (erosion_rate+Gamma_neutron*CRNp.lambda_26Al);

	Conc_26Al = Conc_26Al*Al_exp +  CRNp.S_t*Al_exp*CRNp.P0_26Al*sum_term;
}

double CRN_tParticle::apparent_erosion_26Al_neutron_only(double rho, CRN_parameters& CRNp)
{
	double rho_r = rho/1000;
	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double app_eros;

	app_eros = Gamma_neutron*(CRNp.S_t*CRNp.P0_26Al/Conc_26Al-CRNp.lambda_26Al)/rho_r;
	return -app_eros/100.0;
}

void CRN_tParticle::update_26Al_SSfull(double erosion_rate, CRN_parameters& CRNp)
{
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_26Al[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_26Al);
	}

	Conc_26Al = CRNp.S_t*CRNp.P0_26Al*sum_term;
}

void CRN_tParticle::update_14C_conc(double dt,double erosion_rate, CRN_parameters& CRNp)
{
	double C_exp = exp(-dt*CRNp.lambda_14C);

	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_14C[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
		           (exp(dt*erosion_rate/CRNp.Gamma[i])-
		            exp(-dt*CRNp.lambda_14C))/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_14C);
	}

	Conc_14C = Conc_14C*C_exp +  CRNp.S_t*C_exp*CRNp.P0_14C*sum_term;
}

void CRN_tParticle::update_14C_conc_neutron_only(double dt,double erosion_rate, CRN_parameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2

	double C_exp = exp(-dt*CRNp.lambda_14C);

	double sum_term = 0;
	sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
		           (exp(dt*erosion_rate/Gamma_neutron)-
		            exp(-dt*CRNp.lambda_14C))/
		           (erosion_rate+Gamma_neutron*CRNp.lambda_14C);

	Conc_14C = Conc_14C*C_exp +  CRNp.S_t*C_exp*CRNp.P0_14C*sum_term;
}

double CRN_tParticle::apparent_erosion_14C_neutron_only(double rho, CRN_parameters& CRNp)
{
	double rho_r = rho/1000;			// in a/g/yr
	double Gamma_neutron=CRNp.Gamma[0];			// in g/cm^2
	double app_eros;

	//cout << "LINE 468 Conc 14C: " << Conc_14C << endl;
	//cout << "P0: " << P0_14C << " lambda: " << lambda_14C << endl;
	app_eros = Gamma_neutron*(CRNp.S_t*CRNp.P0_14C/Conc_14C-CRNp.lambda_14C)/rho_r;
	return -app_eros/100.0;
}

void CRN_tParticle::update_14C_SSfull(double erosion_rate, CRN_parameters& CRNp)
{
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_14C[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_14C);
	}

	Conc_14C = CRNp.S_t*CRNp.P0_14C*sum_term;
}

void CRN_tParticle::update_36Cl_conc(double dt,double erosion_rate, CRN_parameters& CRNp)
{
	double Cl_exp = exp(-dt*CRNp.lambda_36Cl);

	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_36Cl[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])*
		           (exp(dt*erosion_rate/CRNp.Gamma[i])-
		            exp(-dt*CRNp.lambda_36Cl))/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_36Cl);
	}

	Conc_36Cl = Conc_36Cl*Cl_exp +  CRNp.S_t*Cl_exp*CRNp.P0_36Cl*sum_term;
}

void CRN_tParticle::update_36Cl_conc_neutron_only(double dt,double erosion_rate, CRN_parameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2

	double Cl_exp = exp(-dt*CRNp.lambda_36Cl);

	double sum_term = 0;
	sum_term+= (exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron)*
		           (exp(dt*erosion_rate/Gamma_neutron)-
		            exp(-dt*CRNp.lambda_36Cl))/
		           (erosion_rate+Gamma_neutron*CRNp.lambda_36Cl);


	Conc_36Cl = Conc_36Cl*Cl_exp +  CRNp.S_t*Cl_exp*CRNp.P0_36Cl*sum_term;
}

double CRN_tParticle::apparent_erosion_36Cl_neutron_only(double rho, CRN_parameters& CRNp)
{
	double rho_r = rho/1000;
	double Gamma_neutron = CRNp.Gamma[0];			// in g/cm^2
	double app_eros;

	app_eros = Gamma_neutron*(CRNp.S_t*CRNp.P0_36Cl/Conc_36Cl-CRNp.lambda_36Cl)/rho_r;
	return -app_eros/100.0;
}

void CRN_tParticle::update_36Cl_SSfull(double erosion_rate, CRN_parameters& CRNp)
{
	double sum_term = 0;
	for (int i = 0; i<4; i++)
	{
		sum_term+= (CRNp.F_36Cl[i]*exp(-effective_dLoc/CRNp.Gamma[i])*CRNp.Gamma[i])/
		           (erosion_rate+CRNp.Gamma[i]*CRNp.lambda_36Cl);
	}

	Conc_36Cl = CRNp.S_t*CRNp.P0_36Cl*sum_term;
}

void CRN_tParticle::update_21Ne_conc(double dt,double erosion_rate, CRN_parameters& CRNp)
{
	double Gamma_neutron= CRNp.Gamma[0];					// in g/cm^2

	Conc_21Ne = Conc_21Ne +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_21Ne*
	             (exp(dt*erosion_rate/Gamma_neutron) - 1)/erosion_rate;
}

double CRN_tParticle::apparent_erosion_21Ne(double rho, CRN_parameters& CRNp)
{
	double rho_r = rho/1000;
	double Gamma_neutron = CRNp.Gamma[0];			// in g/cm^2
	double app_eros;

	app_eros = Gamma_neutron*(CRNp.S_t*CRNp.P0_21Ne/Conc_21Ne)/rho_r;
	return -app_eros/100.0;
}

void CRN_tParticle::update_21Ne_SSfull(double erosion_rate, CRN_parameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];					// in g/cm^2

	//cout << endl << endl <<"BEFORE Conc_21Ne: " << Conc_21Ne << endl;
	if (erosion_rate*erosion_rate < 0.0000000001)
	{
		Conc_21Ne = 0;
	}
	else
	{
		Conc_21Ne = CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_21Ne/erosion_rate;
	}
	//cout << "AFTER Conc_21Ne: " << Conc_21Ne << endl;
}

void CRN_tParticle::update_3He_conc(double dt,double erosion_rate, CRN_parameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];	// in g/cm^2
	Conc_3He = Conc_3He +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_3He*
	             (exp(dt*erosion_rate/Gamma_neutron) - 1)/erosion_rate;
}

double CRN_tParticle::apparent_erosion_3He(double rho, CRN_parameters& CRNp)
{
	double rho_r = rho/1000;
	double Gamma_neutron = CRNp.Gamma[0];	// in g/cm^2
	double app_eros;

	app_eros = Gamma_neutron*(CRNp.S_t*CRNp.P0_3He/Conc_3He)/rho_r;
	return -app_eros/100.0;
}

void CRN_tParticle::update_3He_SSfull(double erosion_rate, CRN_parameters& CRNp)
{
	double Gamma_neutron = CRNp.Gamma[0];	// in g/cm^2


	if (erosion_rate*erosion_rate < 0.0000000001)
	{
		Conc_3He = 0;
	}
	else
	{
		Conc_3He = Conc_3He +  CRNp.S_t*exp(-effective_dLoc/Gamma_neutron)*Gamma_neutron*CRNp.P0_3He/erosion_rate;
	}
}

void CRN_tParticle::update_all_CRN(double dt, double erosion_rate, CRN_parameters& CRNp)
{
	//cout << "LINE 445 CRN_tParticle.cpp updating CRN" << endl;
	update_10Be_conc(dt,erosion_rate, CRNp);
	update_26Al_conc(dt,erosion_rate, CRNp);
	update_36Cl_conc(dt,erosion_rate, CRNp);
	update_14C_conc(dt,erosion_rate, CRNp);
	update_21Ne_conc(dt,erosion_rate, CRNp);
	update_3He_conc(dt,erosion_rate, CRNp);
}

void CRN_tParticle::update_all_CRN_SSfull(double erosion_rate, CRN_parameters& CRNp)
{
	//cout << "LINE 445 CRN_tParticle.cpp updating CRN" << endl;
	update_10Be_SSfull(erosion_rate, CRNp);
	update_26Al_SSfull(erosion_rate, CRNp);
	update_36Cl_SSfull(erosion_rate, CRNp);
	update_14C_SSfull(erosion_rate, CRNp);
	update_21Ne_SSfull(erosion_rate, CRNp);
	update_3He_SSfull(erosion_rate, CRNp);
}

void CRN_tParticle::update_all_CRN_neutron_only(double dt, double erosion_rate, CRN_parameters& CRNp)
{
	update_10Be_conc_neutron_only(dt,erosion_rate, CRNp);
	update_26Al_conc_neutron_only(dt,erosion_rate, CRNp);
	update_36Cl_conc_neutron_only(dt,erosion_rate, CRNp);
	update_14C_conc_neutron_only(dt,erosion_rate, CRNp);
	update_21Ne_conc(dt,erosion_rate, CRNp);
	update_3He_conc(dt,erosion_rate, CRNp);
}

// this updates fallout radionuclides
// NOTE!!!!
// the untis of M_supply are in atoms/cm^2
// and conc is in atoms per gram
// we convert rho into g/cm^3
// adn deltad into cm
// withing the CRN_tParticle object
// the units of k_f10Be here are cm^2/g
void CRN_tParticle::update_fallout10Be_simple_density(double dt, double M_supply_surface,
					double rho_skg, double k_f10Be, double deltad_m, CRN_parameters& CRNp)
{
	// first find which depth interval the particle is in
	int depth_interval = int(dLoc/deltad_m);
	double d_top = double(depth_interval)*deltad_m*100;
	double d_bottom = double(depth_interval+1)*deltad_m*100;
	double deltad = deltad_m*100;
	// the factor of 100 is to convert to cm

	// convert density to g/cm^3
	double rho_s = rho_skg/1000;

	// get the cutoff depth
	double cutoff_depth = 5/(rho_s*k_f10Be);

	if (dLoc*100 > cutoff_depth)
	{
		Conc_f10Be += - Conc_f10Be*CRNp.lambda_10Be;
	}
	else
	{
		Conc_f10Be += dt*M_supply_surface*( exp(-k_f10Be*rho_s*d_top) -exp(-k_f10Be*rho_s*d_bottom) )/
		              (deltad*rho_s*one_min_exp_neg_5) - Conc_f10Be*CRNp.lambda_10Be;
//				cout << " added conc: " <<  dt*M_supply_surface*( exp(-k_f10Be*rho_s*d_top)
//								-exp(-k_f10Be*rho_s*d_bottom) )/
//				              (deltad*rho_s*one_min_exp_neg_5) << endl;
	}
}

// this updates fallout radionuclides
// NOTE!!!!
// the untis of M_supply are in atoms/cm^2
// and conc is in atoms per gram
// we convert rho into g/cm^3
// adn deltad into cm
// withing the CRN_tParticle object
// the units of k_f10Be here are cm^2/g
// chi_f10Be is the fraction of shallow fallout 10Be
void CRN_tParticle::update_fallout10Be_simple_density_2exp(double dt, double M_supply_surface,
					double rho_skg, double k1_f10Be, double k2_f10Be, double chi_f10Be, double deltad_m, CRN_parameters& CRNp)
{
	// first find which depth interval the particle is in
	int depth_interval = int(dLoc/deltad_m);
	double d_top = double(depth_interval)*deltad_m*100;
	double d_bottom = double(depth_interval+1)*deltad_m*100;
	double deltad = deltad_m*100;
	// the factor of 100 is to convert to cm

	// convert density to g/cm^3
	double rho_s = rho_skg/1000;

	// get the cutoff depth for k1
	double cutoff_depth1 = 5/(rho_s*k1_f10Be);
	double cutoff_depth2 = 5/(rho_s*k2_f10Be);

	if (dLoc*100 > cutoff_depth2)
	{
		Conc_f10Be += - Conc_f10Be*CRNp.lambda_10Be;
	}
	if (dLoc*100 < cutoff_depth2 && dLoc*100 > cutoff_depth1)
	{
		Conc_f10Be += dt*M_supply_surface*(1-chi_f10Be)*( exp(-k2_f10Be*rho_s*d_top) -exp(-k2_f10Be*rho_s*d_bottom) )/
		              (deltad*rho_s*one_min_exp_neg_5) - Conc_f10Be*CRNp.lambda_10Be;
	}
	else
	{
		Conc_f10Be += dt*M_supply_surface*(1-chi_f10Be)*( exp(-k2_f10Be*rho_s*d_top) -exp(-k2_f10Be*rho_s*d_bottom) )/
		               (deltad*rho_s*one_min_exp_neg_5) +
		              dt*M_supply_surface*chi_f10Be*( exp(-k1_f10Be*rho_s*d_top) -exp(-k1_f10Be*rho_s*d_bottom) )/
		               (deltad*rho_s*one_min_exp_neg_5) -
		              Conc_f10Be*CRNp.lambda_10Be;
		//		cout << " added conc: " <<  dt*M_supply_surface*( exp(-k_f10Be*rho_s*d_top)
		//						-exp(-k_f10Be*rho_s*d_bottom) )/
		//		              (deltad*rho_s*one_min_exp_neg_5) << endl;
	}
}

// this function resets the depth and effective depth
void CRN_tParticle::update_depths(double d, double ed)
{
	dLoc=d;
	effective_dLoc=ed;
}

// this function changes the effective depth (i.e., the shielding depth)
// for a constant rate of erosion
void CRN_tParticle::erode_mass_only(double dt, double mass_erosion_rate)
{
	effective_dLoc-= mass_erosion_rate*dt;
}

// this function changes the effective depth (i.e., the shielding depth)
// for erosion that is canging linearly in time
void CRN_tParticle::erode_mass_only_linear_increase(double dt, double mass_erosion_rate, double alpha)
{
	effective_dLoc-= 0.5*dt*(alpha*dt+2*mass_erosion_rate);
}

// this function simply resets zeta
void CRN_tParticle::update_zetaLoc(double new_zeta)
{
	zetaLoc = new_zeta;
}

// this function resets the zeta locations using an updated
// surface elevation that preserves the d and effective d locations
// it is only for use with a model that has no soil
void CRN_tParticle::update_zetaLoc_with_new_surface(double new_zeta)
{
	zetaLoc = new_zeta-dLoc;
}


// fuctions for dealing with volumes and surface areas
// it updates the volume based on the current mass, recalucaltes the
// surface area, and returns the volume
//
// volume is in m^3
// mass is in kg
// surface area is in m^2
// density is in kg/m^3
double CRN_tParticle::update_surface_area_and_get_volume(VolumeParticleInfo vpi)
{
	double new_surface_area;
	new_surface_area = vpi.return_surface_area(Type, GSDType, Mass);
	SurfaceArea = new_surface_area;

	double Volume = Mass/vpi.get_type_density(Type);
	return Volume;
}


// this weathers the particle using a listvec
// it updates surface area and mass
double CRN_tParticle::weather_particle(VolumeParticleInfo vpi,
										list< vector<double> >& loss_per_surface_area)
{
	list< vector<double> >::iterator loss_per_surface_area_iter;
	loss_per_surface_area_iter = loss_per_surface_area.begin();
	double mass_loss;
	for (int i = 0; i<Type; i++)
	{
		loss_per_surface_area_iter++;
	}
	mass_loss = (*loss_per_surface_area_iter)[CellIndex] * SurfaceArea;
	Mass = Mass - mass_loss;
	if(Mass < 0)
	{
    Mass = 0;
  }

	double new_surface_area;
	new_surface_area = vpi.return_surface_area(Type, GSDType, Mass);
	SurfaceArea = new_surface_area;

	double Volume = Mass/vpi.get_type_density(Type);
	return mass_loss;
	
	if (Mass == 0)
	{ 
    mass_loss = -99;
  }
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=





