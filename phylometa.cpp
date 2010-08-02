//#include "stdafx.h"

#include <algorithm>
#include <functional>
#include <list>
#include <utility>
#include <stdexcept>
#include <fstream.h>
#include <ymath.h>
#include <stdlib.h>
#include <float.h>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <time.h>

//#include "statistics.h"

#define WANT_TIME
#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
                                     // newmatap.h will get include.h
#include "Matrix_Class/newmatap.h"                // need matrix applications
#include "Matrix_Class/newmatio.h"                                // need matrix output routines
#include "Random_Class/newran.h"


#include "Tree_Class/tree_msvc.h"

#include "TOOLS_Class/tools.h"
#include "statistics.h"
#include "phylogeny.h"
#include "comparativeAnalysis.h"
#include "metaAnalysis.h"

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30


struct MODGROUP
{
	ColumnVector y;
	ColumnVector X;
	Matrix V;
	std::vector < std::string > theSpecies;
};

struct THEQ
{
	double Qfixed;
	double Qrandom;
	int k;
};


double poolAll(ColumnVector y, ColumnVector var) 
{
	Matrix variance, V(y.Nrows(),y.Nrows()); V = 0.0;
	IdentityMatrix I(y.Nrows());
	ColumnVector pooled, X(y.Nrows()); X = 1.0;
	
	for(int i = 1; i < y.Nrows() + 1; i++) V(i,i) = var(i);

	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Matrix Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	// old Q cout << y.t() * V.i() * y - pooled.t() * (X.t() * V.i() * X) * pooled << endl;
	cout << "[fixed]  " << pooled(1) << " (" << variance(1,1) << ") | LVI = " << pooled(1) - 1.96 * sqrt(variance(1,1)) << ", UVI = "<< pooled(1) + 1.96 * sqrt(variance(1,1)) << " | Qw = " << Qw(1,1) << ", df = " << (y.Nrows() - 1) << ", p = " << CHIsq_Prob(Qw(1,1), y.Nrows() - 1) << endl;


	double c = (V.i()).Trace() - ((X.t() * (V * V).i() * X)*(X.t() * V.i() * X).i()).Trace();
	double between = (Qw(1,1) - (y.Nrows() - 1)) / c;
	if(Qw(1,1) < (y.Nrows() - X.Ncols())) between = 0.0;

	V = between * I + V;
	
	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	cout << y.t() * V.i() * y - pooled.t() * (X.t() * V.i() * X) * pooled << endl;
	cout << "[random] " << pooled(1) << " (" << variance(1,1) << ") | LVI = " << pooled(1) - 1.96 * sqrt(variance(1,1)) << ", UVI = "<< pooled(1) + 1.96 * sqrt(variance(1,1)) << " | between = " << between << endl << endl;

	Matrix temp; temp =  y.t() * V.i() * y - pooled.t() * (X.t() * V.i() * X) * pooled;
	double Qrand = temp(1,1);
	return Qrand;
}

void poolAll_metaRegression(ColumnVector y, ColumnVector var, ColumnVector predictor, bool includeIntercept) 
{
	int k = y.Nrows();

	Matrix variance, V(k,k), X; V = 0.0;
	IdentityMatrix I(k);
	ColumnVector pooled, P(k); 

	for(int i = 1; i < k + 1; i++) {
		V(i,i) = var(i);
		P(i) = predictor(i);
	}

	if (includeIntercept == true) { 
		ColumnVector INTERCEPT; INTERCEPT = 1.0;
		X = INTERCEPT | P;
	}
	else X = P;


	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Matrix Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	cout << "[fixed] intercept = " << pooled(1) << " (" << variance(1,1) << ") | LVI = " << pooled(1) - 1.96 * sqrt(variance(1,1)) << ", UVI = "<< pooled(1) + 1.96 * sqrt(variance(1,1)) << " | Z = " << pooled(1)/sqrt(variance(1,1)) << ", p = " << Zscore(pooled(1)/sqrt(variance(1,1))) << endl;
	cout << "            slope = " << pooled(2) << " (" << variance(2,2) << ") | LVI = " << pooled(2) - 1.96 * sqrt(variance(2,2)) << ", UVI = "<< pooled(2) + 1.96 * sqrt(variance(2,2)) << " | Z = " << pooled(2)/sqrt(variance(2,2)) << ", p = " << Zscore(pooled(1)/sqrt(variance(2,2))) << endl;

	/*double c = (V.i()).Trace() - ((X.t() * (V * V).i() * X)*(X.t() * V.i() * X).i()).Trace();
	double between = (Qw(1,1) - (k - 1)) / c;
	if(Qw(1,1) < (k - X.Ncols())) between = 0.0;

	V = between * I + V;
	
	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	cout << "[random] " << pooled(1) << " (" << variance(1,1) << ") | LVI = " << pooled(1) - 1.96 * sqrt(variance(1,1)) << ", UVI = "<< pooled(1) + 1.96 * sqrt(variance(1,1)) << " | between = " << between << endl << endl;
	*/
}

void poolAll_metaRegression_MATRIX(ColumnVector y, ColumnVector predictor, bool includeIntercept, Matrix V) 
{
	int k = y.Nrows();

	Matrix variance, X;
	IdentityMatrix I(k);
	ColumnVector pooled, P(k); 

	for(int i = 1; i < k + 1; i++) P(i) = predictor(i);

	if (includeIntercept == true) { 
		ColumnVector INTERCEPT(k); INTERCEPT = 1.0;
		X = INTERCEPT | P;
	}
	else X = P;

	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Matrix Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 

	cout << "[fixed] intercept = " << pooled(1) << " (" << variance(1,1) << ") | LCI = " << pooled(1) - 1.96 * sqrt(variance(1,1)) << ", UCI = " << pooled(1) + 1.96 * sqrt(variance(1,1)) << " | Z = " << absolute(pooled(1))/sqrt(variance(1,1)) << ", p = " << Zscore(absolute(pooled(1))/sqrt(variance(1,1))) << endl;
	cout << "            slope = " << pooled(2) << " (" << variance(2,2) << ") | LCI = " << pooled(2) - 1.96 * sqrt(variance(2,2)) << ", UCI = " << pooled(2) + 1.96 * sqrt(variance(2,2)) << " | Z = " << absolute(pooled(2))/sqrt(variance(2,2)) << ", p = " << Zscore(absolute(pooled(2))/sqrt(variance(2,2))) << endl;
}

int countNumberWithinGroups(std::vector< std::string > group, std::string theGroup)
{
	int N = 0;
	for(int i = 0; i < group.size(); i++) if(group[i] == theGroup) N++;
	return N;
}


std::vector< MODGROUP > seperateGroups(ColumnVector y, ColumnVector var, std::vector< std::string > group)
{
	std::vector < MODGROUP > allTheGroups;
	MODGROUP newGroup;
	std::vector< std::string > theGroups = cleanStringVector(group);
	
	for(int j = 0; j < theGroups.size(); j++) {
		int groupSize = countNumberWithinGroups(group, theGroups[j]), temp = 1;
		ColumnVector yG(groupSize), xG(groupSize); yG = 0.0; xG = 0.0;
		Matrix VG(groupSize, groupSize); VG = 0.0;
		for(int i = 0; i < group.size(); i++) {
			if(group[i] == theGroups[j]) {
				//cout << i << "->" << j << endl;
				yG(temp) = y(i+1);
				xG(temp) = 1;
				VG(temp, temp) = var(i+1);
				temp++;
			}
		}
		newGroup.y = yG;
		newGroup.X = xG;
		newGroup.V = VG;
		allTheGroups.push_back(newGroup);
	}
	return allTheGroups;
}


void poolAllByGroup(ColumnVector y, ColumnVector var, std::vector< std::string > group) 
{
	std::vector< std::string > theGroups = cleanStringVector(group);
	Matrix X(group.size(), theGroups.size()); X = 0.0;
	for(int i = 0; i < group.size(); i++) {
		for(int j = 0; j < theGroups.size(); j++) {
			if(group[i] == theGroups[j]) X(i+1,j+1) = 1;
			else X(i+1,j+1) = 0;
		}
	}
	
	Matrix variance, V(y.Nrows(),y.Nrows()); V = 0.0;
	IdentityMatrix I(y.Nrows());
	ColumnVector pooled, Ones(y.Nrows()); Ones = 1.0;

	for(i = 1; i < y.Nrows() + 1; i++) V(i,i) = var(i);

	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Matrix Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	Matrix Qt = y.t() * (V.i() - V.i() * Ones * (Ones.t() * V.i() * Ones).i() * Ones.t() * V.i()) * y; 

	cout << "[fixed] by group" << endl;
	for(i = 1; i <= pooled.Nrows(); i++) {
		std::cout << "  " << theGroups[i-1] << " ";
		cout << pooled(i) << " (" << variance(i,i) << ") | LVI = " << pooled(i) - 1.96 * sqrt(variance(i,i)) << ", UVI = "<< pooled(i) + 1.96 * sqrt(variance(i,i)) << endl;
	}
	cout << "  Qw = " << Qw(1,1) << ", df = " << (y.Nrows() - 1) << ", p = " << CHIsq_Prob(Qw(1,1), y.Nrows() - 1) << endl;
	cout << "  Qb = " << Qt(1,1) - Qw(1,1) << ", df = " << 1 << ", p = " << CHIsq_Prob(Qt(1,1) - Qw(1,1), 1) << endl;;

	double c = (V.i()).Trace() - ((X.t() * (V * V).i() * X)*(X.t() * V.i() * X).i()).Trace();
	double between = (Qw(1,1) - (y.Nrows() - X.Ncols())) / c;
	if(Qw(1,1) < (y.Nrows() - X.Ncols())) between = 0.0;

	V = between * I + V;
	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	Qt = y.t() * (V.i() - V.i() * Ones * (Ones.t() * V.i() * Ones).i() * Ones.t() * V.i()) * y; 

	cout << "[random] by group (between var = " << between << ")" << endl;
	for(i = 1; i <= pooled.Nrows(); i++) {
		std::cout << "  " << theGroups[i-1] << " ";
		cout << pooled(i) << " (" << variance(i,i) << ") | LVI = " << pooled(i) - 1.96 * sqrt(variance(i,i)) << ", UVI = "<< pooled(i) + 1.96 * sqrt(variance(i,i)) << endl;
	}
	cout << "  Qb = " << Qt(1,1) - Qw(1,1) << ", df = " << 1 << ", p = " << CHIsq_Prob(Qt(1,1) - Qw(1,1), 1) << endl;

	



}





double poolAll_MATRIX(ColumnVector y, Matrix V) 
{
	Matrix variance;
	IdentityMatrix I(y.Nrows());
	ColumnVector pooled, X(y.Nrows()); X = 1.0;

	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Matrix Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	cout << "[fixed]  " << pooled(1) << " (" << variance(1,1) << ") | LVI = " << pooled(1) - 1.96 * sqrt(variance(1,1)) << ", UVI = "<< pooled(1) + 1.96 * sqrt(variance(1,1)) << " | Qw = " << Qw(1,1) << ", df = " << (y.Nrows() - 1) << ", p = " << CHIsq_Prob(Qw(1,1), y.Nrows() - 1) << endl;

	double c = (V.i()).Trace() - ((X.t() * (V * V).i() * X)*(X.t() * V.i() * X).i()).Trace ();
	double between = (Qw(1,1) - (y.Nrows() - 1)) / c;
	if(Qw(1,1) < (y.Nrows() - X.Ncols())) between = 0.0;
	V = between * I + V;
	
	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	cout << y.t() * V.i() * y - pooled.t() * (X.t() * V.i() * X) * pooled << endl;
	cout << "[random] " << pooled(1) << " (" << variance(1,1) << ") | LVI = " << pooled(1) - 1.96 * sqrt(variance(1,1)) << ", UVI = "<< pooled(1) + 1.96 * sqrt(variance(1,1)) << " | between = " << between << endl << endl;
	
	Matrix temp; temp =  y.t() * V.i() * y - pooled.t() * (X.t() * V.i() * X) * pooled;
	double Qrand = temp(1,1);
	return Qrand;

}

void poolAllByGroup_MATRIX(ColumnVector y, Matrix V, std::vector< std::string > group) 
{
	std::vector< std::string > theGroups = cleanStringVector(group);
	Matrix X(group.size(), theGroups.size()); X = 0.0;
	for(int i = 0; i < group.size(); i++) {
		for(int j = 0; j < theGroups.size(); j++) {
			if(group[i] == theGroups[j]) X(i+1,j+1) = 1;
			else X(i+1,j+1) = 0;
		}
	}
	
	Matrix variance;
	IdentityMatrix I(y.Nrows());
	ColumnVector pooled, Ones(y.Nrows()); Ones = 1.0;

	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Matrix Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	Matrix Qt = y.t() * (V.i() - V.i() * Ones * (Ones.t() * V.i() * Ones).i() * Ones.t() * V.i()) * y; 

	cout << "[fixed] by group" << endl;
	for(i = 1; i <= pooled.Nrows(); i++) {
		std::cout << "  " << theGroups[i-1] << " ";
		cout << pooled(i) << " (" << variance(i,i) << ") | LVI = " << pooled(i) - 1.96 * sqrt(variance(i,i)) << ", UVI = "<< pooled(i) + 1.96 * sqrt(variance(i,i)) << endl;
	}
	cout << "  Qw = " << Qw(1,1) << ", df = " << (y.Nrows() - 1) << ", p = " << CHIsq_Prob(Qw(1,1), y.Nrows() - 1) << endl;
	cout << "  Qb = " << Qt(1,1) - Qw(1,1) << ", df = " << 1 << ", p = " << CHIsq_Prob(Qt(1,1) - Qw(1,1), 1) << endl;;


	double c = (V.i()).Trace() - ((X.t() * (V * V).i() * X)*(X.t() * V.i() * X).i()).Trace();
	double between = (Qw(1,1) - (y.Nrows() - X.Ncols())) / c;
	
	if(Qw(1,1) < (y.Nrows() - X.Ncols())) between = 0.0;

	V = between * I + V;
	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	Qt = y.t() * (V.i() - V.i() * Ones * (Ones.t() * V.i() * Ones).i() * Ones.t() * V.i()) * y; 

	cout << "[random] by group (between var = " << between << ")" << endl;
	for(i = 1; i <= pooled.Nrows(); i++) {
		std::cout << "  " << theGroups[i-1] << " ";
		cout << pooled(i) << " (" << variance(i,i) << ") | LVI = " << pooled(i) - 1.96 * sqrt(variance(i,i)) << ", UVI = "<< pooled(i) + 1.96 * sqrt(variance(i,i)) << endl;
	}
	cout << "  Qb = " << Qt(1,1) - Qw(1,1) << ", df = " << 1 << ", p = " << CHIsq_Prob(Qt(1,1) - Qw(1,1), 1) << endl;;
}

void singleLine(int space, int length) 
{
	for(int i = 0; i < space; i++) cout << " ";
	for(i = 0; i < length; i++) cout << (char)196;
	cout << endl;
}

void doubleLine(int length) 
{
	for(int i = 0; i < length; i++) cout << (char)205;
	cout << endl << endl;
}

THEQ homoSUMMARY(ColumnVector y, ColumnVector var, std::vector< std::string > group) 
{
	THEQ theQs;
	std::vector< std::string > theGroups = cleanStringVector(group);
	Matrix X(group.size(), theGroups.size()); X = 0.0;

	for(int i = 0; i < group.size(); i++) {
		for(int j = 0; j < theGroups.size(); j++) {
			if(group[i] == theGroups[j]) X(i+1,j+1) = 1;
			else X(i+1,j+1) = 0;
		}
	}

	theQs.k = y.Nrows();
	Matrix variance, V(y.Nrows(),y.Nrows()); V = 0.0;
	IdentityMatrix I(y.Nrows());
	ColumnVector pooled, Ones(y.Nrows()); Ones = 1.0;
	
	for(i = 1; i < y.Nrows() + 1; i++) V(i,i) = var(i);

	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Matrix Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	theQs.Qfixed = Qw(1,1);
	Matrix Qt = y.t() * (V.i() - V.i() * Ones * (Ones.t() * V.i() * Ones).i() * Ones.t() * V.i()) * y; 

	//SUMMARY HOMOGENEITY TABLE
	cout << "  " << "TABLE 1. Summary of fit statistics." << endl;
	singleLine(2, 42);
	cout << "  " << "  Source" << setw(17) << "Q" << setw(7) << "df" << setw(5) << "p" << endl;
	singleLine(2, 42);
	//std::cout << floatToString(CHIsq_Prob(Qt(1,1) - Qw(1,1), 1),4) << std::endl;
	//cout << ((CHIsq_Prob(Qt(1,1) - Qw(1,1), 1) < 0.0001) ? "<0.0001" : stringToChar(floatToString(CHIsq_Prob(Qt(1,1) - Qw(1,1), 1),4))) << endl;
	cout << "  " << setw(18) << setiosflags(ios::fixed) << "Between groups  " << setw(9) << setprecision(2) << Qt(1,1) - Qw(1,1) << setw(5) << setprecision(0) << 1 << setw(8) << setprecision(4) << CHIsq_Prob(Qt(1,1) - Qw(1,1), 1) << endl;
	cout << "  " << setw(18) << "Within groups   " << setw(9) << setprecision(2) << Qw(1,1)  << setw(5) << setprecision(0) << y.Nrows() - 1 << setw(8) << setprecision(4) << CHIsq_Prob(Qw(1,1), y.Nrows() - 1) << endl;
	
	std::vector< MODGROUP > allGroups = seperateGroups(y, var, group);
	for(i = 0; i < theGroups.size(); i++) {
		ColumnVector yg = allGroups[i].y;
		Matrix Vg = allGroups[i].V;
		ColumnVector Xg = allGroups[i].X;
		Matrix Qwg = yg.t() * (Vg.i() - Vg.i() * Xg * (Xg.t() * Vg.i() * Xg).i() * Xg.t() * Vg.i()) * yg; 
		std::cout << "  " << std::setw(17) << "Within group " << theGroups[i];
		cout << setw(9) << setprecision(2) << Qwg(1,1) << setw(5) << (allGroups[i].y.Nrows() - 1) << setw(8) << setprecision(4) << CHIsq_Prob(Qwg(1,1), allGroups[i].y.Nrows() - 1) << endl;
	}
	cout << "  " << setw(18) << "Total           " << setw(9) << setprecision(2) << Qt(1,1)  << endl;
	cout << "  " << endl;


	//RANDOM EFFECTS
	double c = (V.i()).Trace() - ((X.t() * (V * V).i() * X)*(X.t() * V.i() * X).i()).Trace();
	double between = (Qw(1,1) - (y.Nrows() - X.Ncols())) / c;
	if(Qw(1,1) < (y.Nrows() - X.Ncols())) between = 0.0;

	Matrix Vr = between * I + V;
	ColumnVector pooledr = (X.t() * Vr.i() * X).i() * X.t() * Vr.i() * y;
	Matrix variancer = (X.t() * Vr.i() * X).i();
	Matrix Qwr = y.t() * (Vr.i() - Vr.i() * X * (X.t() * Vr.i() * X).i() * X.t() * Vr.i()) * y; 
	theQs.Qrandom= Qwr(1,1);
	Matrix Qtr = y.t() * (Vr.i() - Vr.i() * Ones * (Ones.t() * Vr.i() * Ones).i() * Ones.t() * Vr.i()) * y; 

	//SUMMARY TABLE CONTINUED
	cout << "  " << setw(18) << "Between groups  "   << endl;
	cout << "  " << setw(18) << "random-effects" << setw(9) << setprecision(2) << Qtr(1,1) - Qwr(1,1) << setw(5) << 1 << setw(8) << setprecision(4) << CHIsq_Prob(Qtr(1,1) - Qwr(1,1), 1) << endl;
	singleLine(2, 42);
	cout << endl;

	//EFFECTS SUMMARY
	cout << endl << endl << "  " << "TABLE 2. Summary of pooled effect sizes ("<<(char)235<<") and variances (" <<(char)229 << (char)253 <<")." << endl;
	singleLine(2, 75);
	cout << "  " << "                                                           Is non-zero?" << endl;
	singleLine(57, 20);
	cout << "  " << "  Group          " << "k" << "      " << (char)235 << "       " << (char)229 << (char)253 << "         95%CI          Z     df     p" << endl;
	singleLine(2, 75);
	cout << "  " << setw(13) << "Fixed-effects" << endl << endl;

	ColumnVector totalPooled = (Ones.t() * V.i() * Ones).i() * Ones.t() * V.i() * y;
	Matrix totalVariance = (Ones.t() * V.i() * Ones).i();

	Matrix Qr = totalPooled.t() * Ones.t() * V.i() * Ones * totalPooled;
	cout << "  " << setw(13) << "  All studies"   << setprecision(3) << setw(5) << Ones.Nrows() << setw(9) << totalPooled(1) << setprecision(4) << setw(9) << totalVariance(1,1) << setw(5) << setprecision(3) << "(" << (totalPooled(1) - 1.96 * sqrt(totalVariance(1,1))) << "," << (totalPooled(1) + 1.96 * sqrt(totalVariance(1,1))) << ")" << setw(8) << setprecision(2) << Qr(1,1) << setw(5) << 1 << setw(9) << setprecision(4) << CHIsq_Prob(Qr(1,1), 1) << endl;
	for(i = 0; i < theGroups.size(); i++) {
		ColumnVector yg = allGroups[i].y;
		Matrix Vg = allGroups[i].V;
		ColumnVector Xg = allGroups[i].X;
		ColumnVector groupPooled = (Xg.t() * Vg.i() * Xg).i() * Xg.t() * Vg.i() * yg; 
		Matrix groupVariance = (Xg.t() * Vg.i() * Xg).i();
		Matrix groupQr = groupPooled.t() * Xg.t() * Vg.i() * Xg * groupPooled;
		std::cout << "  " << std::setw(10) << "  Group " << theGroups[i];
		cout << setprecision(3) << setw(7) << Xg.Nrows() << setw(9) << groupPooled(1) << setprecision(4) << setw(9) << groupVariance(1,1) << setw(17) << setw(5) << setprecision(3) << "(" << (groupPooled(1) - 1.96 * sqrt(groupVariance(1,1))) << "," << (groupPooled(1) + 1.96 * sqrt(groupVariance(1,1))) << ")" << setw(8) << setprecision(2) << groupQr(1,1) << setw(5) << 1 << setw(9) << setprecision(4) << CHIsq_Prob(groupQr(1,1), 1) << endl;
	}

	cout << endl << "  " << setw(13) << "Random-effects" << endl << endl;

	//RANDOM EFFECTS 2
	double allc = (V.i()).Trace() - ((X.t() * (V * V).i() * X)*(X.t() * V.i() * X).i()).Trace();
	Matrix Qwr2 = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	double between2 = (Qwr2(1,1) - (y.Nrows() - X.Ncols())) / allc;
	if(Qwr2(1,1) < (y.Nrows() - X.Ncols())) between2 = 0.0;
	Matrix allVr = between2 * I + V;
	ColumnVector allpooledr = (Ones.t() * allVr.i() * Ones).i() * Ones.t() * allVr.i() * y;
	Matrix allvariancer = (Ones.t() * allVr.i() * Ones).i();
	Matrix allQr = allpooledr.t() * Ones.t() * allVr.i() * Ones * allpooledr;
	cout << "  " << setw(13) << "  All studies"   << setprecision(3) << setw(5) << Ones.Nrows() << setw(9) << allpooledr(1) << setprecision(4) << setw(9) << allvariancer(1,1) << setw(5) << setprecision(3) << "(" << (allpooledr(1) - 1.96 * sqrt(allvariancer(1,1))) << "," << (allpooledr(1) + 1.96 * sqrt(allvariancer(1,1))) << ")" << setw(8) << setprecision(2) << allQr(1,1) << setw(5) << 1 << setw(9) << setprecision(4) << CHIsq_Prob(allQr(1,1), 1) << endl;
		for(i = 0; i < theGroups.size(); i++) {
		IdentityMatrix Ig(allGroups[i].y.Nrows()); 
		ColumnVector yg = allGroups[i].y;
		Matrix Vg = between * Ig + allGroups[i].V;
		ColumnVector Xg = allGroups[i].X;
		ColumnVector groupPooled = (Xg.t() * Vg.i() * Xg).i() * Xg.t() * Vg.i() * yg; 
		Matrix groupVariance = (Xg.t() * Vg.i() * Xg).i();
		Matrix groupQr = groupPooled.t() * Xg.t() * Vg.i() * Xg * groupPooled;
		std::cout << "  " << std::setw(10) << "  Group " << theGroups[i];
		cout << setprecision(3) << setw(7) << Xg.Nrows() << setw(9) << groupPooled(1) << setprecision(4) << setw(9) << groupVariance(1,1) << setw(17) << setw(5) << setprecision(3) << "(" << (groupPooled(1) - 1.96 * sqrt(groupVariance(1,1))) << "," << (groupPooled(1) + 1.96 * sqrt(groupVariance(1,1))) << ")" << setw(8) << setprecision(2) << groupQr(1,1) << setw(5) << 1 << setw(9) << setprecision(4) << CHIsq_Prob(groupQr(1,1), 1) << endl;
	}

	singleLine(2, 75);
	cout << endl;
	return theQs;
}

std::vector < Phylogeny > getPhyloModGroups(Phylogeny thePhylo)
{
	std::vector< std::string > theGroups = cleanStringVector(thePhylo.format_phylo_to_DATACOLUMN_moderators()), 
							   allMods = thePhylo.format_phylo_to_DATACOLUMN_moderators(),
							   theSpecies = thePhylo.getAllSpecies(),
							   tempCollection;
	std::vector <Phylogeny > thePhylos;
	
	for(int j = 0; j < theGroups.size(); j++) {
		for(int i = 0; i < theSpecies.size(); i++) if(theGroups[j] == allMods[i]) tempCollection.push_back(theSpecies[i]);
		thePhylos.push_back(thePhylo.treeSpeciesSubset(tempCollection));
		//printStringVector(tempCollection);
		tempCollection.clear();
	}
	
	return thePhylos;
}

THEQ homoPHYLOSUMMARY(ColumnVector y, ColumnVector var, std::vector< std::string > group, Phylogeny thePhylo) 
{
	THEQ theQs;
	std::vector< std::string > theGroups = cleanStringVector(group);
	Matrix X(group.size(), theGroups.size()); X = 0.0;
	std::vector< Phylogeny > thePhyloByGroups = getPhyloModGroups(thePhylo);

	for(int i = 0; i < group.size(); i++) {
		for(int j = 0; j < theGroups.size(); j++) {
			if(group[i] == theGroups[j]) X(i+1,j+1) = 1;
			else X(i+1,j+1) = 0;
		}
	}

	theQs.k = y.Nrows();
	Matrix variance, V(y.Nrows(),y.Nrows()), Vnormal(y.Nrows(),y.Nrows()); V = 0.0; Vnormal = 0.0;
	IdentityMatrix I(y.Nrows());
	ColumnVector pooled, Ones(y.Nrows()); Ones = 1.0;
	
	for(i = 1; i < y.Nrows() + 1; i++) Vnormal(i,i) = var(i);
	V = thePhylo.metaCovariance(true, false, 0.0);
	
	pooled = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	variance = (X.t() * V.i() * X).i();
	Matrix Qw = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	theQs.Qfixed = Qw(1,1);
	Matrix Qwnormal = y.t() * (Vnormal.i() - Vnormal.i() * X * (X.t() * Vnormal.i() * X).i() * X.t() * Vnormal.i()) * y; 
	Matrix Qt = y.t() * (V.i() - V.i() * Ones * (Ones.t() * V.i() * Ones).i() * Ones.t() * V.i()) * y; 

	//SUMMARY HOMOGENEITY TABLE
	cout << "  " << "TABLE 1. Summary of fit statistics." << endl;
	singleLine(2, 57);
	cout << "  " << "                                           Adjusted via" << endl;
	cout << "  " << "                                           # polytomies" << endl;
	singleLine(44, 15);
	cout << "  " << "  Source" << setw(17) << "Q" << setw(7) << "df" << setw(5) << "p" << setw(9) << "df" << setw(5) << "p" << endl;
	//	cout << "  " << "  Source" << setw(17) << "Q" << setw(7) << "df" << setw(5) << "p" << endl;
	singleLine(2, 57);
	cout << "  " << setw(18) << "Between groups  " << setw(9) << setprecision(2) << Qt(1,1) - Qw(1,1) << setw(5) << setprecision(0) << 1 << setw(8) << setprecision(4) << CHIsq_Prob(Qt(1,1) - Qw(1,1), 1) << endl;
	cout << "  " << setw(18) << "Within groups   " << setw(9) << setprecision(2) << Qw(1,1)  << setw(5) << setprecision(0) << y.Nrows() - 1 << setw(8) << setprecision(4) << CHIsq_Prob(Qw(1,1), y.Nrows() - 1) << setw(6) << setprecision(0) << (y.Nrows() - thePhylo.stats_numberOfPolytomies() - 1) << setw(8) << setprecision(4) << CHIsq_Prob(Qw(1,1), y.Nrows() - thePhylo.stats_numberOfPolytomies() - 1)<< endl;
	
	std::vector< MODGROUP > allGroups = seperateGroups(y, var, group);
	
	for(i = 0; i < theGroups.size(); i++) {
		ColumnVector yg = allGroups[i].y;
		Matrix Vg = thePhyloByGroups[i].metaCovariance(true, false, 0.0);
		ColumnVector Xg = allGroups[i].X;
		Matrix Qwg = yg.t() * (Vg.i() - Vg.i() * Xg * (Xg.t() * Vg.i() * Xg).i() * Xg.t() * Vg.i()) * yg; 
		std::cout << "  " << std::setw(17) << "Within group " <<  theGroups[i];
		cout << setw(9) << setprecision(2) << Qwg(1,1) << setw(5) << (allGroups[i].y.Nrows() - 1) << setw(8) << setprecision(4) << CHIsq_Prob(Qwg(1,1), allGroups[i].y.Nrows() - 1) << setw(6) << setprecision(0) << (allGroups[i].y.Nrows() - thePhyloByGroups[i].stats_numberOfPolytomies() - 1) << setw(8) << setprecision(4) << CHIsq_Prob(Qwg(1,1), allGroups[i].y.Nrows() - thePhyloByGroups[i].stats_numberOfPolytomies() - 1) << endl;
	}
	cout << "  " << setw(18) << "Total           " << setw(9) << setprecision(2) << Qt(1,1)  << endl;
	cout << "  " << endl;


	//RANDOM EFFECTS
	double c = (V.i()).Trace() - ((X.t() * (V * V).i() * X)*(X.t() * V.i() * X).i()).Trace();
	double between = (Qwnormal(1,1) - (y.Nrows() - X.Ncols())) / c;
	if(Qwnormal(1,1) < (y.Nrows() - X.Ncols())) between = 0.0;

	Matrix Vr = between * I + V;
	ColumnVector pooledr = (X.t() * Vr.i() * X).i() * X.t() * Vr.i() * y;
	Matrix variancer = (X.t() * Vr.i() * X).i();
	Matrix Qwr = y.t() * (Vr.i() - Vr.i() * X * (X.t() * Vr.i() * X).i() * X.t() * Vr.i()) * y; 
	theQs.Qrandom = Qwr(1,1);
	Matrix Qtr = y.t() * (Vr.i() - Vr.i() * Ones * (Ones.t() * Vr.i() * Ones).i() * Ones.t() * Vr.i()) * y; 

	//SUMMARY TABLE CONTINUED
	cout << "  " << setw(18) << "Between groups  "   << endl;
	cout << "  " << setw(18) << "random-effects" << setw(9) << setprecision(2) << Qtr(1,1) - Qwr(1,1) << setw(5) << 1 << setw(8) << setprecision(4) << CHIsq_Prob(Qtr(1,1) - Qwr(1,1), 1) << endl;
	singleLine(2, 57);
	cout << endl;

	//EFFECTS SUMMARY
	cout << endl << endl << "  " << "TABLE 2. Summary of pooled effect sizes ("<<(char)235<<") and variances (" <<(char)229 << (char)253 <<")." << endl;
	singleLine(2, 75);
	cout << "  " << "                                                           Is non-zero?" << endl;
	singleLine(57, 20);
	cout << "  " << "  Group          " << "k" << "      " << (char)235 << "       " << (char)229 << (char)253 << "         95%CI          Z     df     p" << endl;
	singleLine(2, 75);
	cout << "  " << setw(13) << "Fixed-effects" << endl << endl;

	ColumnVector totalPooled = (Ones.t() * V.i() * Ones).i() * Ones.t() * V.i() * y;
	Matrix totalVariance = (Ones.t() * V.i() * Ones).i();

	Matrix Qr = totalPooled.t() * Ones.t() * V.i() * Ones * totalPooled;
	cout << "  " << setw(13) << "  All studies"   << setprecision(3) << setw(5) << Ones.Nrows() << setw(9) << totalPooled(1) << setprecision(4) << setw(9) << totalVariance(1,1) << setw(5) << setprecision(3) << "(" << (totalPooled(1) - 1.96 * sqrt(totalVariance(1,1))) << "," << (totalPooled(1) + 1.96 * sqrt(totalVariance(1,1))) << ")" << setw(8) << setprecision(2) << Qr(1,1) << setw(5) << 1 << setw(9) << setprecision(4) << CHIsq_Prob(Qr(1,1), 1) << endl;
	for(i = 0; i < theGroups.size(); i++) {
		ColumnVector yg = allGroups[i].y;
		Matrix Vg = thePhyloByGroups[i].metaCovariance(true, false, 0.0);
		ColumnVector Xg = allGroups[i].X;
		ColumnVector groupPooled = (Xg.t() * Vg.i() * Xg).i() * Xg.t() * Vg.i() * yg; 
		Matrix groupVariance = (Xg.t() * Vg.i() * Xg).i();
		Matrix groupQr = groupPooled.t() * Xg.t() * Vg.i() * Xg * groupPooled;
		std::cout << "  " << std::setw(10) << "  Group " << theGroups[i];
		cout << setprecision(3) << setw(7) << Xg.Nrows() << setw(9) << groupPooled(1) << setprecision(4) << setw(9) << groupVariance(1,1) << setw(17) << setw(5) << setprecision(3) << "(" << (groupPooled(1) - 1.96 * sqrt(groupVariance(1,1))) << "," << (groupPooled(1) + 1.96 * sqrt(groupVariance(1,1))) << ")" << setw(8) << setprecision(2) << groupQr(1,1) << setw(5) << 1 << setw(9) << setprecision(4) << CHIsq_Prob(groupQr(1,1), 1) << endl;
	}

	cout << endl << "  " << setw(13) << "Random-effects" << endl << endl;

	//RANDOM EFFECTS 2
	double allc = (V.i()).Trace() - ((X.t() * (V * V).i() * X)*(X.t() * V.i() * X).i()).Trace();
	Matrix Qwr2 = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	double between2 = (Qwr2(1,1) - (y.Nrows() - X.Ncols())) / allc;
	if(Qwr2(1,1) < (y.Nrows() - X.Ncols())) between2 = 0.0;
	Matrix allVr = between2 * I + V;
	ColumnVector allpooledr = (Ones.t() * allVr.i() * Ones).i() * Ones.t() * allVr.i() * y;
	Matrix allvariancer = (Ones.t() * allVr.i() * Ones).i();
	Matrix allQr = allpooledr.t() * Ones.t() * allVr.i() * Ones * allpooledr;
	cout << "  " << setw(13) << "  All studies"   << setprecision(3) << setw(5) << Ones.Nrows() << setw(9) << allpooledr(1) << setprecision(4) << setw(9) << allvariancer(1,1) << setw(5) << setprecision(3) << "(" << (allpooledr(1) - 1.96 * sqrt(allvariancer(1,1))) << "," << (allpooledr(1) + 1.96 * sqrt(allvariancer(1,1))) << ")" << setw(8) << setprecision(2) << allQr(1,1) << setw(5) << 1 << setw(9) << setprecision(4) << CHIsq_Prob(allQr(1,1), 1) << endl;
	for(i = 0; i < theGroups.size(); i++) {
		IdentityMatrix Ig(allGroups[i].y.Nrows()); 
		ColumnVector yg = allGroups[i].y;
		Matrix Vg = between * Ig + thePhyloByGroups[i].metaCovariance(true, false, 0.0);
		ColumnVector Xg = allGroups[i].X;
		ColumnVector groupPooled = (Xg.t() * Vg.i() * Xg).i() * Xg.t() * Vg.i() * yg; 
		Matrix groupVariance = (Xg.t() * Vg.i() * Xg).i();
		Matrix groupQr = groupPooled.t() * Xg.t() * Vg.i() * Xg * groupPooled;
		std::cout << "  " << std::setw(10) << "  Group " << theGroups[i];
		cout << setprecision(3) << setw(7) << Xg.Nrows() << setw(9) << groupPooled(1) << setprecision(4) << setw(9) << groupVariance(1,1) << setw(17) << setw(5) << setprecision(3) << "(" << (groupPooled(1) - 1.96 * sqrt(groupVariance(1,1))) << "," << (groupPooled(1) + 1.96 * sqrt(groupVariance(1,1))) << ")" << setw(8) << setprecision(2) << groupQr(1,1) << setw(5) << 1 << setw(9) << setprecision(4) << CHIsq_Prob(groupQr(1,1), 1) << endl;
	}
	singleLine(2, 75);
	cout << endl;
	return theQs;
}

/*void homoPHYLOSUMMARY(void) 
{
	cout << "  " << "TABLE 1. Summary of fit statistics." << endl;
	singleLine(2, 57);
	cout << "  " << "                                          Adjusted via" << endl;
	cout << "  " << "                                          # polytomies" << endl;
	singleLine(43, 14);
	cout << "  " << "  Source" << setw(17) << "Q" << setw(7) << "df" << setw(5) << "p" << setw(9) << "df" << setw(5) << "p" << endl;
	singleLine(2, 57);
	cout << "  " << setw(18) << "Between groups  " << setw(9) << 20.30  << setw(5) << 2   << setw(8) << 0.2536 << setw(6) << 2   << setw(8) << 0.2536 <<endl;
	cout << "  " << setw(18) << "Within groups   " << setw(9) << 20.35  << setw(5) << 7   << setw(8) << 0.2536 << setw(6) << 2   << setw(8) << 0.2536 <<endl;
	cout << "  " << setw(18) << "Within group 1"   << setw(9) << 120.35 << setw(5) << 1   << setw(8) << 0.2536 << setw(6) << 2   << setw(8) << 0.2536 <<endl;
	cout << "  " << setw(18) << "Within group 2"   << setw(9) << 20.00  << setw(5) << 100 << setw(8) << 0.2536 << setw(6) << 100   << setw(8) << 0.2536 <<endl;
	cout << "  " << setw(18) << "Within group 3"   << setw(9) << 20.35  << setw(5) << 6   << setw(8) << 0.2536 << setw(6) << 2   << setw(8) << 0.2536 <<endl;
	cout << "  " << setw(18) << "Total           " << setw(9) << 20.35  << setw(5) << 10  << setw(8) << 0.2536 << setw(6) << 2   << setw(8) << 0.2536 <<endl;
	singleLine(2, 57);
	cout << endl;
}
*/

void effectsSUMMARY(void) 
{
	cout << "  " << "TABLE 2. Summary of pooled effect sizes ("<<(char)235<<") and variances (" <<(char)229 << (char)253 <<")." << endl;
	singleLine(2, 70);
	cout << "  " << "                                                     Is non-zero?" << endl;
	singleLine(51, 20);
	cout << "  " << "  Group          " << (char)235 << "       " << (char)229 << (char)253 << "         95%CI          Z     df     p" << endl;
	singleLine(2, 70);
	cout << "  " << setw(13) << "All studies"   << setw(7) << 0.335  << setw(9) << 0.02359 << setw(17) << "(0.336,-0.325)" << setw(8) << "35.3" << setw(5) << 10 << setw(9) << 0.2356 <<endl;
	cout << "  " << setw(13) << "  Group 1  "   << setw(7) << 0.335  << setw(9) << 0.02359 << setw(17) << "(0.336,-0.325)" << setw(8) << "5.3"  << setw(5) << 5 << setw(9) << 0.2356 << endl;
	cout << "  " << setw(13) << "  Group 2  "   << setw(7) << 0.335  << setw(9) << 0.02359 << setw(17) << "(0.336,-0.325)" << setw(8) << "35.3" << setw(5) << 3 << setw(9) << 0.2356 << endl;
	cout << "  " << setw(13) << "  Group 3  "   << setw(7) << 0.335  << setw(9) << 0.02359 << setw(17) << "(0.336,-0.325)" << setw(8) << "35.3" << setw(5) << 1 << setw(9) << 0.2356 << endl;
	singleLine(2, 70);
	cout << endl;
}

double getAIC(double Q, int k, bool isTraditional)
{
	double m = 1.0, AIC = 0.0; if(isTraditional == false) m = 2.0;
	//AIC = ((2.0 * m) + (2.0 * log((k /(2.0))*(log(2.0 * PI * (Q/k)) + 1.0))));
	AIC = (2.0 * m) + (k * (log(2.0 * PI * (Q / k)) + 1.0));
	return AIC;
}


void comparisonSUMMARY(THEQ theQtraditional, THEQ theQphylogenetic) 
{
	//cout << theQtraditional.Qfixed << " " <<  theQtraditional.Qrandom << " " << theQtraditional.k << endl;
	//cout << theQphylogenetic.Qfixed << " " <<  theQphylogenetic.Qrandom << " " << theQphylogenetic.k << " " << PI << endl;
		
	cout << "  " << "TABLE 1. Summary of model fit." << endl;
	singleLine(2, 51);
	cout << "  " << "  Meta-analysis" << setw(28) << "AIC" << endl;
	singleLine(35, 18);
	cout << "  " << "               " << setw(34) << "fixed    random" << endl;
	singleLine(2, 51);
	cout << "  " << setw(30) << "Traditional                 " << setprecision(2) << setw(9) << getAIC(theQtraditional.Qfixed, theQtraditional.k, true)  << setw(9) << getAIC(theQtraditional.Qrandom, theQtraditional.k, true) << endl;
	cout << "  " << setw(30) << "Phylogenetically-independent" << setprecision(2) << setw(9) << getAIC(theQphylogenetic.Qfixed, theQphylogenetic.k, false)  << setw(9) << getAIC(theQphylogenetic.Qrandom, theQphylogenetic.k, false) << endl;//	"*" << endl;
	singleLine(2, 51);
	cout << "  " << "Note. Lowest AIC is best fit." << endl;
	cout << endl;
}

double mean_meta_MATRIX(ColumnVector y, Matrix V) 
{
	ColumnVector pooled, X(y.Nrows()); X = 1.0;
	Matrix temp; temp = (X.t() * V.i() * X).i() * X.t() * V.i() * y;
	return temp(1,1);
}

double variance_meta_MATRIX(ColumnVector y, Matrix V) 
{
	ColumnVector pooled, X(y.Nrows()); X = 1.0;
	Matrix temp; temp = (X.t() * V.i() * X).i();
	return temp(1,1);
}

double SSE_meta_MATRIX(ColumnVector y, Matrix V) 
{
	ColumnVector pooled, X(y.Nrows()); X = 1.0;
	Matrix temp; temp = y.t() * (V.i() - V.i() * X * (X.t() * V.i() * X).i() * X.t() * V.i()) * y; 
	return temp(1,1);
}

double L_meta(double variance, double SSE, Matrix V)
{
	return exp(SSE * (-1.0 / (2.0 * variance))) * sqrt(pow(2.0 * PI * variance, V.Nrows()) * V.Determinant());
}

Matrix lambda_MATRIX(double lambda, Matrix P) 
{
	Matrix tempP; tempP = P; 
	for(int i = 1; i < P.Ncols() + 1; i++)
		for(int j = 1; j < P.Nrows() + 1; j++) if (i != j) tempP(i,j) = lambda * P(i,j);
	return tempP;
}

double lambda_meta(ColumnVector y, Matrix V, Matrix P)
{
	double MIN = 0.0001, MAX = 2, smallestL = 0.0, L = 0.0, var = 0.0, theSSE = 0.0, theLAmbda = 0.0;
	Matrix tempV, tempP, tempW;

	for(double i = MIN; i <= MAX; i+=MIN) {
		tempP = lambda_MATRIX(i, P);
		tempW = V * tempP * V;
		var = variance_meta_MATRIX(y, tempW);
		theSSE = SSE_meta_MATRIX(y, tempW);
		L = L_meta(var, theSSE, tempW);
		if(L > smallestL) {
			smallestL = L;
			theLAmbda = i;
		}
	}
	return theLAmbda;
}

int main (void)
{
	myBanner("phyloMeta", "1.0", "2009/06/20");
	
	//DATA SOURCES
	cout << "DATA FILES. Input two sources: (1) phylogeny and (2) effect size data." << endl;
	doubleLine(80);	
	Phylogeny thePhylo;
	char fileNamePhylogeny[100], fileNameMeta[100], pause[100];

	// // INPUT PHYLOGENY
	std::cout << "  (1) Please enter the filename containing the phylogeny." << std::endl << "      (e.g., phylo.txt): ";
	std::cin >> fileNamePhylogeny;
	if(!FileExist(fileNamePhylogeny)) {std::cout << std::endl << "ERROR: " << fileNamePhylogeny << " does not exist in directory..." << std::endl; std::cin >> pause; exit(0); }
	thePhylo.file_Open_NEWICKPhylogeny(fileNamePhylogeny);
	std::cout << std::endl;
	thePhylo.print_phylogeny();
	std::cout << std::endl;
	
	// // INPUT EFFECT SIZE DATA
	std::cout << std::endl << "  (2) Now the filename containing the effect size data." << std::endl << "      (e.g., meta.txt): ";
	std::cin >> fileNameMeta;
	if(!FileExist(fileNameMeta)) {std::cout << std::endl << "ERROR: " << fileNameMeta << " does not exist in directory..." << std::endl; std::cin >> pause; exit(0); }
	
//	cout << endl << endl;/
//	cout << "DATA SECTION B. Summary before analysis." << endl;
//	doubleLine(80);	
	
	std::vector< std::string > speciesList;
	char aSpecies[100];
	std::string tempSp;

	double anEffect, aVar;
	std::vector< double > theEffects, theVariances;
	int theMod;
	std::vector< int > modList;

	ifstream dataFile(fileNameMeta);
	while(dataFile >> aSpecies >> anEffect >> aVar >> theMod) {
		tempSp = aSpecies;
		speciesList.push_back(tempSp);
		theEffects.push_back(anEffect);
		theVariances.push_back(aVar);
		modList.push_back(theMod);
	}

	thePhylo.inputTaxaXdata(speciesList, theEffects, modList, theVariances);

	//thePhylo.print_taxaData();
	

	//poolAll(thePhylo.format_phylo_to_DATACOLUMN(), thePhylo.format_phylo_to_DATACOLUMN_VAR());
	
	//printStringVector(thePhylo.format_phylo_to_DATACOLUMN_moderators());
	//poolAllByGroup(thePhylo.format_phylo_to_DATACOLUMN(), thePhylo.format_phylo_to_DATACOLUMN_VAR(), thePhylo.format_phylo_to_DATACOLUMN_moderators());
	
	THEQ theQtraditional, theQphylogenetic;

	cout << endl << endl;
	// TRADITIONAL META-ANALYSIS
	cout << "RESULTS SECTION A. Traditional meta-analysis." << endl;
	doubleLine(80);
	theQtraditional = homoSUMMARY(thePhylo.format_phylo_to_DATACOLUMN(), thePhylo.format_phylo_to_DATACOLUMN_VAR(), thePhylo.format_phylo_to_DATACOLUMN_moderators());
	cout << endl;
	//effectsSUMMARY();

	// PHYLOGENETICALLY-INDEPENDENT META-ANALYSIS
	cout << endl << endl;
	cout << "RESULTS SECTION B. Phylogenetically-independent meta-analysis." << endl;
	doubleLine(80);
	theQphylogenetic = homoPHYLOSUMMARY(thePhylo.format_phylo_to_DATACOLUMN(), thePhylo.format_phylo_to_DATACOLUMN_VAR(), thePhylo.format_phylo_to_DATACOLUMN_moderators(), thePhylo);
	cout << endl;
	//effectsSUMMARY();
	
	// MODEL FIT BETWEEN TRADITIONAL AND PHYLOGENETICALLY-INDEPENDENT META-ANALYSIS
	cout << endl << endl;
	cout << "RESULTS SECTION C. Traditional vs. phylogenetically-independent meta-analysis." << endl;
	doubleLine(80);
	comparisonSUMMARY(theQtraditional, theQphylogenetic);
	
	
	cout << "press any key to exit... ";
	cin >>  pause;
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/*Phylogeny aPhylo;
	aPhylo.file_Open_NEWICKPhylogeny("Jennions_2008.phy");
	Matrix theP = aPhylo.matrix_standardize_UNQUAL(aPhylo.format_phylo_to_matrix());

	std::vector< std::string > theTAXA;
	std::vector< double > theDATA, theVar;

	double yVal, var, opt;
	int k = 0;
	char sp[100], filename[50] = "jennions_2009.dat";

	ifstream theFile(filename);
	while(theFile >> sp >> yVal >> opt >> var) {
		theTAXA.push_back(sp);
		theDATA.push_back(yVal);
		theVar.push_back(var);
		k++;
	}
	theFile.close();

	ColumnVector effects(k);
	Matrix theSD(k,k); theSD = 0.0;

	for(int i = 0; i < theDATA.size(); i++) {
		effects(i+1) = theDATA[i];
		theSD(i+1, i+1) = sqrt(theVar[i]);
	}
	
	cout << variance_meta_MATRIX(effects, theSD * theP * theSD) << endl;
	cout << SSE_meta_MATRIX(effects, theSD * theP * theSD) << endl;
	
	
	cout << lambda_meta(effects, theSD, theP) << endl;
	cout << aPhylo.stats_Pagel_lambda(effects, false, false) << endl;

*/
	return 0;
}
