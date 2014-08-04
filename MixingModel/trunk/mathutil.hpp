//-*-c++-*-

/*********************************************************************\
**  mathutil.h:  Header file for special math utilities not in math.h.
**               All or most routines from Numerical Recipes in C by
**               Press et al.
**
**  $Id: mathutil.h,v 1.4 2002/07/08 17:21:49 arnaud Exp $
\*********************************************************************/

#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

#define PI 3.1415926

double ran3( long *idum );
double ran0( long *idum );
double expdev( long *idum);
int bin(double start, double end, double value, double Nbins);
vector<double> make_Hist_bins(double start_data, double end_data, int N_bins);
vector<int> make_Hist(vector<double> bins, vector<int> data);
string itoa(int num);
string dtoa(double num);
void split_string(const string& str, const string& delim, vector<string>& output);
int pnpoly(vector<double>& vertx, vector<double>& verty, double testx, double testy);
double erfi( double x);
double calculate_MLE(vector<double>& measured, vector<double>& modelled, vector<double>& sigma);
double calculate_MLE(vector<double>& measured, vector<double>& modelled, double sigma);
double calclulate_area_of_tetrahedron(vector<double> p1, vector<double> p2, vector<double> p3, vector<double> p4);

void matlab_double_sort(
  vector<double> & unsorted,
  vector<double> & sorted,
  vector<size_t> & index_map);

// Act like matlab's Y = X[I]
// where I contains a vector of indices so that after,
// Y[j] = X[I[j]] for index j
// this implies that Y.size() == I.size()
// X and Y are allowed to be the same reference

void matlab_double_reorder(
  vector<double> & unordered,
  vector<size_t> const & index_map,
  vector<double> & ordered);





#endif
