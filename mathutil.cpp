/**************************************************************************\
**  mathutil.cpp:  Special math routines not in math libraries. Most or all
**                 from Numerical Recipes in C by Press et al.
**
**  $Id: mathutil.cpp,v 1.2 1999/01/12 21:03:01 gtucker Exp $
\**************************************************************************/
#ifndef MATHUTIL_CPP
#define MATHUTIL_CPP


#include "mathutil.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>
#include <math.h>
using namespace std;

/*********************************************************\
**  ran3
**
**  Random number generator from Numerical Recipes.
**  Returns a uniform random number between 0.0 and 1.0.
**  Set idum to any negative value to initialize or
**  reinitialize the sequence.
**
**  Parameters: idum - random seed
**
\*********************************************************/

#define MBIG  990000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum)
{
   //cout << &idum << endl;
   static int inext,inextp;
   static long ma[56];
   static int iff=0;
   long mj,mk;
   int i,ii,k;

   if (*idum < 0 || iff == 0) {
      iff=1;
      if(*idum>0)
       *idum = - *idum;
      mj=MSEED+ *idum;
      if (mj<0)
       mj = -mj;
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
         ii=(21*i) % 55;
         ma[ii]=mk;
         mk=mj-mk;
         if (mk < MZ) mk += MBIG;
         mj=ma[ii];
      }
      for (k=1;k<=4;k++)
          for (i=1;i<=55;i++) {
             ma[i] -= ma[1+(i+30) % 55];
             if (ma[i] < MZ) ma[i] += MBIG;
          }
      inext=0;
      inextp=31;
      *idum=1;
   }
   if (++inext == 56) inext=1;
   if (++inextp == 56) inextp=1;
   mj=ma[inext]-ma[inextp];
   if (mj < MZ) mj += MBIG;
   ma[inext]=mj;
   return fabs(mj*FAC);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran0(long *idum)
 {
  long k;
  double ans;

  *idum ^=MASK;
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  ans = AM*(*idum);
  *idum ^= MASK;
  return ans;
 }

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK


double expdev( long *idum)
 {
  double dum;

  do
   dum = ran3(idum);
  while (dum == 0.0);
  return -log(dum);
 }

int bin(double start,double end, double val, int Nbins)
{
 double binwidth = (end-start)/double(Nbins);
 //cout << "binwidth is: " << binwidth << endl;
 double bin_loc = val/binwidth;
 //cout << "val is: " << val << endl;
 int bl = int(bin_loc);
 if (bl == Nbins)
  bl = Nbins-1;

 //cout << "bin_loc is: " << bin_loc << " and the bin is: " << bl << endl;
 return bl;
}


// the vector that is returned by this function has N_bin nodes, each node is the location of the centerpoint
// of each bin. The width of the bin may be found by subtracing adjacent centerpoints.
vector<double> make_Hist_bins(double start_data, double end_data, int N_bins)
 {
  double range = end_data-start_data;
  double bin_width = range/double(N_bins);
  vector<double> bins(N_bins);
  for (int i = 0; i<N_bins; i++)
   bins[i] = start_data+i*bin_width+ bin_width*0.5;
  return bins;
 }


vector<int> make_Hist(vector<double> bins, vector<int> data)
 {
  // get the bin spacing
  double bin_width = bins[1]-bins[0];

  // find wher the bins start
  double data_start= bins[0]-bin_width*0.5;

  // get the size of teh data vector and the bin vector
  int data_sz = data.size();
  int N_bins = bins.size();

  // create the counting vector
  vector<int> binned_data(N_bins,0);

  double bin_num;
  int bn;
  // now bin the data
  for (int i = 0; i< data_sz; i++)
   {

    bin_num = (data[i]-data_start)/double(bin_width);
    bn = int(bin_num);					// truncate (the left border of the bin is inclusicve,
    							// the right border is not)

    if (bn<1)
     bn = 1;
    if (bn> N_bins)
     bn = N_bins;

    //cout << "age: " << data[i] << " bin: " << bn << endl;
    binned_data[bn-1]++;
   }

  return binned_data;
 }

string itoa(int num)
{
    stringstream converter;
    converter << num;
    return converter.str();
}

string dtoa(double num)
{
    stringstream converter;
    converter << num;
    return converter.str();
}

void split_string(const string& str, const string& delim, vector<string>& output)
{
    size_t start = 0, found = str.find_first_of(delim);

    while (found != string::npos)
    {
        if (found > start)
            output.push_back( str.substr(start, found - start) );

        start = ++found;
        found = str.find_first_of(delim, found);
    }
    if (start < str.size())
        output.push_back( str.substr(start) );
}

int pnpoly(vector<double>& vertx, vector<double>& verty, double testx, double testy)
{
  int i, j, c = 0;
  int nvert = vertx.size();
  for (i = 0, j = nvert-1; i < nvert; j = i++)
  {
	  if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	       (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
	  {
      	c = !c;
	  }
  }
  return c;
}


// calcualte the imaginary error function using trapezoid rule integration
double erfi(double tau)
{
	double sum = 0;
	double N_nodes = 100000;
	double dtau = tau/N_nodes;
	//double mini_dtau = dtau/4.0;
	//double term1,term2, term3, term4, term5;
	//double a;
	//double front_term = dtau/90.0;
	int ii;

	// trapezoidal rule
	double xloc_old;
	double xloc_new = 0.0;
	for (ii= 0; ii<N_nodes; ii++)
	{
		xloc_old = xloc_new;
		xloc_new = double(ii)*dtau;
		sum+= 0.5*dtau*(exp(xloc_old*xloc_old) + exp(xloc_new*xloc_new) );
	}

	// simpson's
	//	double xloc_old;
	//	double xloc_new = 0.0;
	//	double xloc_mid;
	//	front_term = dtau/6;
	//	for (ii= 0; ii<N_nodes; ii++)
	//	{
	//		xloc_old = xloc_new;
	//		xloc_new = double(ii)*dtau;
	//		xloc_mid = 0.5*(xloc_old+xloc_new);
	//		sum+= front_term*(exp(xloc_old*xloc_old) + 4*exp(xloc_mid*xloc_mid) + exp(xloc_new*xloc_new) );
	//	}


	// boole's
	//	for (ii = 0; ii<N_nodes; ii++);
	//	{
	//		a = dtau*double(ii);
	//		term1 =  7.0*exp(a*a);
	//		term2 = 32.0*exp((a+mini_dtau)*(a+mini_dtau));
	//		term3 = 12.0*exp((a+2.0*mini_dtau)*(a+2.0*mini_dtau));
	//		term4 = 32.0*exp((a+3.0*mini_dtau)*(a+3.0*mini_dtau));
	//		term5 =  7.0*exp((a+4.0*mini_dtau)*(a+4.0*mini_dtau));
	//		sum += front_term*(term1+term2+term3+term4+term5);
	//	}

	return 2*sum/sqrt(M_PI);
}

// get the least squared maximum likelihood estimator
double calculate_MLE(vector<double>& measured, vector<double>& modelled, vector<double>& sigma)
{
	// get the number of samples
	int n_samples = measured.size();
	double MLE_tot = 1;
	for (int i = 0; i<n_samples; i++)
	{
		//cout << "exp term: " << -0.5*(measured[i]-modelled[i])*(measured[i]-modelled[i])/
		//							 sigma[i]*sigma[i] << endl;
		MLE_tot = MLE_tot*exp(-0.5*(measured[i]-modelled[i])*(measured[i]-modelled[i])/
									 sigma[i]*sigma[i]);
	}
	return MLE_tot;
}

// get the least squared maximum likelihood estimator
double calculate_MLE(vector<double>& measured, vector<double>& modelled, double sigma)
{
	// get the number of samples
	int n_samples = measured.size();
	double MLE_tot = 1;
	for (int i = 0; i<n_samples; i++)
	{
		MLE_tot = MLE_tot*exp(-0.5* (measured[i]-modelled[i])*(measured[i]-modelled[i])/
									 sigma*sigma);
	}
	return MLE_tot;
}




// Comparison struct used by sort
// http://bytes.com/topic/c/answers/132045-sort-get-index
template<class T> struct index_cmp
{
  index_cmp(const T arr) : arr(arr) {}
  bool operator()(const size_t a, const size_t b) const
  {
    return arr[a] < arr[b];
  }
  const T arr;
};



void matlab_double_sort(
  vector<double> & unsorted,
  vector<double> & sorted,
  std::vector<size_t> & index_map)
{
  // Original unsorted index map
  index_map.resize(unsorted.size());
  for(size_t i=0;i<unsorted.size();i++)
  {
    index_map[i] = i;
  }
  // Sort the index map, using unsorted for comparison
  sort(
    index_map.begin(),
    index_map.end(),
    index_cmp<std::vector<double>& >(unsorted));

  sorted.resize(unsorted.size());
  matlab_double_reorder(unsorted,index_map,sorted);
}




// This implementation is O(n), but also uses O(n) extra memory
void matlab_double_reorder(
  std::vector<double> & unordered,
  std::vector<size_t> const & index_map,
  std::vector<double> & ordered)
{
  // copy for the reorder according to index_map, because unsorted may also be
  // sorted
  vector<double> copy = unordered;
  ordered.resize(index_map.size());
  for(int i = 0; i< int(index_map.size());i++)
  {
    ordered[i] = copy[index_map[i]];
  }
}

// calcualtes the area of a terahedron. The points have three data entries, the x, y and z locations
double calclulate_area_of_tetrahedron(vector<double> p1, vector<double> p2, vector<double> p3, vector<double> p4)
{
	double area;
	double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
	x1 = p1[0]; x2 = p2[0]; x3 = p3[0]; x4 = p4[0];
	y1 = p1[1]; y2 = p2[1]; y3 = p3[1]; y4 = p4[1];
	z1 = p1[2]; z2 = p2[2]; z3 = p3[2]; z4 = p4[2];
	double line1 =  x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3;
	double line2 = -x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3;
	double line3 =  x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2;
	double line4 = -x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2;

	area = (1.0/6.0)*(line1+line2+line3+line4);
	return area;
}


































#endif
