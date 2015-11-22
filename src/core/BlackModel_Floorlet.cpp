/*
 * BlackFormula.cpp
 *
 *  Created on: Nov 17, 2015
 *      Author: Yuntao Zhang
 */

#include <iostream>
#include <math.h>
#include <vector>
/*********************************************************************************
 normalDistribution: Return N(x)
 [in]: double x
 [out]: double N(x)
 *********************************************************************************/


double normalDistribution(double x)
{
  static const double RT2PI = sqrt(4.0*acos(0.0));
  static const double SPLIT = 10./sqrt(2);
  static const double a[] = {220.206867912376,221.213596169931,112.079291497871,33.912866078383,6.37396220353165,0.700383064443688,3.52624965998911e-02};
  static const double b[] = {440.413735824752,793.826512519948,637.333633378831,296.564248779674,86.7807322029461,16.064177579207,1.75566716318264,8.83883476483184e-02};

  const double z = fabs(x);
  double Nz = 0.0;

  // if z outside these limits then value effectively 0 or 1 for machine precision
  if(z<=37.0)
  {
    // NDash = N'(z) * sqrt{2\pi}
    const double NDash = exp(-z*z/2.0)/RT2PI;
    if(z<SPLIT)
    {
      const double Pz = (((((a[6]*z + a[5])*z + a[4])*z + a[3])*z + a[2])*z + a[1])*z + a[0];
      const double Qz = ((((((b[7]*z + b[6])*z + b[5])*z + b[4])*z + b[3])*z + b[2])*z + b[1])*z + b[0];
      Nz = RT2PI*NDash*Pz/Qz;
    }
    else
    {
      const double F4z = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
      Nz = NDash/F4z;
    }
  }
  return x>=0.0 ? 1-Nz : Nz;
}

/*********************************************************************************
BlacksFormula : computes the price of floorlet using BlackÕs 1976 model
[in]: double f : forward rate
double P : price of pure discount bond
double L : principal amount of bond
double Rcap : cap rate
double vol : market volatility
double tau : length between payment times
double dtau : fixed tenor between reset times
[out] : double : price of cap
********************************************************************************/

double BlacksFormula(double f, double P, double L, double Rcap, double
vol, double tau, double dtau)
{
	double d1 = (log(f) / Rcap) + ((vol*vol)/2)*tau/(vol*sqrt(tau));
	double d2 = d1 - vol*sqrt(tau);
	return P*dtau*L*(-f*normalDistribution(-d2) + Rcap*normalDistribution(-d1));
}

/**********************************************************************************
priceBlackFlo : computes the price of floorlets using BlacksFormula
[in]: vector<double> capVol : vector of cap volatilities
vector<double> PDB : price of pure discount bonds
vector<double> maturity : vector of caplet maturities (payment times)
vector<double> Rcap : cap rate
double L : principal amount of loan
double tenor : length of time between payment times (reset dates)
[out]: vector<double> : caplets
*********************************************************************************/
std::vector<double> priceBlackFlo(std::vector<double> capVol, std::vector<double> PDB,
		std::vector<double> maturity, double Rcap, double L, double tenor)
{
	int i;
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
	// compute yield price
	for (i = 0; i < capVol.size(); i++)
	{
		tmp = -(1/t[i])*(log(P[i]));
		R.push_back(tmp);
	}
	// compute forward rates
	for(i = 0; i < capVol.size()-1; i++)
		tmp = -(1/tenor)*log(P[i+1]/P[i]);
		f.push_back(tmp);
	}
	// compute caplets with Blacks Formula
	for (i = 0; i < capVol.size()-1; i++)
	{
		tmp = BlacksFormula(f[i],P[i],faceValue,Rcap,capV[i],t[i],tenor);
		caplet.push_back(tmp);
	}
	return caplet;
}
