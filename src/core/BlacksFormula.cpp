/*
 * BlackFormula.cpp
 *
 *  Created on: Nov 17, 2015
 *      Author: Yuntao Zhang
 */

#include <iostream>
#include <math.h>
#include <vector>
#include "BlacksFormula.h"

BlacksFormula::BlacksFormula(){}
BlacksFormula::~BlacksFormula() {}

/*********************************************************************************
 normalDistribution: Return N(x)
 [in]: double x
 [out]: double N(x)
 *********************************************************************************/
double BlacksFormula::normalDistribution(double x)
{
	double a1 =  0.254829592;
	   double a2 = -0.284496736;
	   double a3 =  1.421413741;
	   double a4 = -1.453152027;
	   double a5 =  1.061405429;
	   double p  =  0.3275911;

	   // Save the sign of x
	   int sign = 1;
	   if (x < 0)
	       sign = -1;
	   x = fabs(x)/sqrt(2.0);

	   // A&S formula 7.1.26
	   double t = 1.0/(1.0 + p*x);
	   double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

	   return 0.5*(1.0 + sign*y);
}

/*********************************************************************************
BlacksFormulaFunc : computes the price of cap using BlackÕs 1976 model
[in]: double f : forward rate
double P : price of pure discount bond
double L : principal amount of bond
double Rcap : cap rate
double vol : market volatility
double tau : length between payment times
double dtau : fixed tenor between reset times
[out] : double : price of cap
********************************************************************************/

double BlacksFormula::BlacksFormulaFunc(double f, double P, double L, double Rcap, double
vol, double tau, double dtau, bool isCallSwaption)
{
	//int w = isCallSwaption;
	//double d1 = (log(f / Rcap) + (vol*vol/2)*tau)/(vol*sqrt(tau));
	double part1 = log(f / Rcap);
	double part2 = (vol*vol/2)*tau;
	double part3 = vol*sqrt(tau);
	double d1 = (part1 + part2)/part3;
	double d2 = d1 - vol*sqrt(tau);
	if (isCallSwaption){
		double b=Rcap*normalDistribution(d2);
		double a=f*normalDistribution(d1);
		return P*dtau*L*(f*normalDistribution(d1) - Rcap*normalDistribution(d2));
	}else{
		return P*dtau*L*(Rcap*normalDistribution(-d2) - f*normalDistribution(-d1));
	}
}

/**********************************************************************************
priceBlackCap : computes the price of caplets using BlacksFormula
[in]: vector<double> capVol : vector of cap volatilities
vector<double> PDB : price of pure discount bonds
vector<double> maturity : vector of caplet maturities (payment times)
vector<double> Rcap : cap rate
double L : principal amount of loan
double tenor : length of time between payment times (reset dates)
bool isCallSwaption : 1 for call swaption, 0 for put swaption
[out]: vector<double> : caplets
*********************************************************************************/
std::vector<double> BlacksFormula::priceBlackCap(std::vector<double> capVol, std::vector<double> PDB,
		std::vector<double> maturity, double Rcap, double L, double tenor, bool isCallSwaption)
{
	unsigned int i;
	std::vector<double> f; // forward rates
	std::vector<double> R; // yield price
	std::vector<double> capV; // stores caplet volatilities
	std::vector<double> P; // stores pure bond prices
	std::vector<double> t; // payment dates
	std::vector<double>::iterator iter; // vector iterator
	std::vector<double> caplet; // stores caplets
	double faceValue = 0.0; // bond face value
	double tmp = 0.0; // temp variable
	// compute face value
	faceValue = L*(1 + Rcap*tenor);
	// store cap volatilities
	for (iter = capVol.begin(); iter != capVol.end(); iter++)
	{
		tmp = *iter;
		capV.push_back(tmp);
	}
	// compute pure discount bond prices
	for (iter = PDB.begin(); iter != PDB.end(); iter++)
	{
		tmp = *iter;
		P.push_back(tmp);
	}
	// store payment dates
	for (iter = maturity.begin(); iter != maturity.end(); iter++)
	{
		tmp = *iter;
		t.push_back(tmp);
	}
	// compute forward rates
	for (i = 0; i < capVol.size()-1; i++)
	{
		tmp = -(1/t[i])*(log(P[i]));
		R.push_back(tmp);
		tmp = -(1/tenor)*log(P[i+1]/P[i]);
		f.push_back(tmp);
	}
	// compute caplets with Blacks Formula
	for (i = 0; i < capVol.size()-1; i++)
	{
		tmp = BlacksFormulaFunc(f[i],P[i],faceValue,Rcap,capV[i],t[i],tenor,isCallSwaption);
		caplet.push_back(tmp);
	}
	return caplet;
}

