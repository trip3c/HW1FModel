/*
 * BaseModule.h
 *
 *  Created on: Nov 11, 2015
 *      Author: Chandra
 */

#include <vector>
#include <map>

using namespace std;

#ifndef BASEMODULE_H_
#define BASEMODULE_H_

/**
 * This is the parent class for all modules. This class defines the core functionality
 * calculating values for intermediate equations
 */
class BaseModule {
public:
	/**
	 * This data structure holds the input parameters supplied by the user. It is
	 * assumed that all the vectors are of the same length, e.g., the 5th index of
	 * all vectors hold information about the corresponding data point for 5th
	 * time period index.
	 * @param time This variable holds the time period data point
	 * @param forward This variable holds the forward rates
	 * @param priceD This variable holds the discount factor for zero coupon bonds
	 * @param aMeanReversion This variable holds the mean reversion coefficient
	 * @param sigma This variable holds the volatility of the interest rate
	 * @param Et This variable holds the intermediate value of exponential integral
	 *           of the mean reversion
	 */
	struct DataStruct{

		vector<double> time;
		vector<double> forward;
		vector<double> priceD;
		vector<double> aMeanReversion;
		vector<double> sigma;
		vector<double> Et;
	};

	/**
	 * This data structure holds the input swaption volatility and its corresponding
	 * maturity and tenor
	 * @param vol This parameter holds the swaption volatility for one corresponding
	 * maturity and tenor. This vector is made up of matrix of maturity as rows and
	 * tenor as columns. The inner vector holds one row for one maturity, i.e. one
	 * maturity multiple tenors and the outer vector holds multiple rows for different
	 * maturity vectors.
	 * @param tenor This parameter holds the tenor of swaption
	 * @param maturity This parameter holds the maturity of swaption
	 */
	struct SwaptionVolStruct{
		vector<vector<double> > vol;
		vector<double> tenor;
		vector<double> maturity;
	};

	/**
	 * The constructor of the base module class
	 */
	BaseModule();

	/**
	 * The destructor of the base module class
	 */
	virtual ~BaseModule();

	/**
	 * The abstract method of the module, which its inherited class must implement
	 */
	virtual int moduleMainFunc() = 0;

	/**
	 * The reference for input data structure
	 */
	DataStruct data;

	/**
	 * The reference for the swaption volatility
	 */
	SwaptionVolStruct implVol;

	/**
	 * An object that maps time data point as key to index of the time as the value.
	 * This map is used for direct access of the data points
	 */
	map<double, int> timePos;

	void calculateEt();
	int locate(double t);


	double calculateB(double t, double T);
	double calculateA(double t, double T);
	double calculateVr(double t);
	double calculateVp(double Tf, double Tp);

	double erf(double x);
	double N(double x, double mu, double sigma);
	double ZBP(double Tf, double Tp, double X);
	double pSwaption(double K, double T0, double Tp);
	double fun_r(const vector<double>& c, const vector<double>& A, const vector<double>& B, double r);
	double fun_dev(const vector<double>& c, const vector<double>& A, const vector<double>& B, double r);
	double solveR(const vector<double>& c, const vector<double>& A, const vector<double>& B);

	vector<double> theta();
	double takeDev(const vector<double>& value, double point, int n);
};

#endif /* BASEMODULE_H_ */
