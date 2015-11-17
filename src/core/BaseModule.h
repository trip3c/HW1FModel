/*
 * BaseModule.h
 *
 *  Created on: Nov 11, 2015
 *      Author: Chandra
 */

#include <vector>
#include <map>
#include "Constants.h"

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
	 * @param weights This parameter holds the weights for the corresponding swaption
	 *                calibration
	 */
	struct SwaptionVolStruct{
		vector<vector<double> > vol;
		vector<double> tenor;
		vector<double> maturity;
		vector<vector<double> > weights;
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
	SwaptionVolStruct actualVol;

	/**
	 * The reference for constants file.
	 * TODO Move to properties file
	 */
	Constants constants;

	/**
	 * An object that maps time data point as key to index of the time as the value.
	 * This map is used for direct access of the data points
	 */
	map<double, int> timePos;

	/*
	 * This method calculates the exponential integration over the mean reversion for each period
	 */
	void calculateEt();

	/*
	 * This method generates the index of a particular time in the time vector
	 * @param timeP the time to search for
	 * @return the index of the time
	 */
	int locate(double timeP);

	/*
	 * This method calculates the value of B(t) function for Affine structure
	 * @param t The start of bond period
	 * @param T The end of bond period
	 * @return the value of B(t)
	 */
	double calculateB(double t, double T);

	/*
	 * This method calculates the value of A(t) function for Affine structure
	 * @param t The start of bond period
	 * @param T The end of bond period
	 * @return the value of A(t)
	 */
	double calculateA(double t, double T);

	/*
	 * This method calculates the variance of short rate from the initial time
	 * @param t the end of calculation period
	 * @return variance of short rate
	 */
	double calculateVr(double t);

	/*
	 * This method calculates the variance of bond ratio P(0,Tf)/P(0,Tp)
	 * @param Tf The maturity of denominator bond
	 * @param Tp The maturity of nominator bond
	 * @return the of variance bond ratio
	 */
	double calculateVp(double Tf, double Tp);

	/**
	 * This method calculates the value of Gauss error function, which is a helper function for calculating Normal cdf
	 * @param x the input of error function
	 * @return the value of error function
	 */
	double erf(double x);

	/**
	 * This method approximates the normal cumulative distribution function
	 * @param x the critical value
	 * @param mu the mean of normal distribution
	 * @param sigma the standard deviation of normal distribution
	 * @return the cdf probability
	 */
	double N(double x, double mu, double sigma);

	/**
	 * This method calculates the value of a zero-bond put option
	 * @param Tf the maturity of the option
	 * @param Tp the maturity of the bond
	 * @param X Strike of the option
	 * @return the value of put option
	 */
	double ZBP(double Tf, double Tp, double X);

	/**
	 * This method calculates the value of a payer swaption
	 * @param K the strike rate of the swaption
	 * @param T0 the maturity of the swaption
	 * @param Tp the swap tenor
	 * @return the value of swaption
	 */
	double pSwaption(double K, double T0, double Tp);

	/**
	 * This is a helper function for the solveR()
	 */
	double fun_r(const vector<double>& c, const vector<double>& A, const vector<double>& B, double r);

	/**
	 * This is another helper function for the solveR()
	 */
	double fun_dev(const vector<double>& c, const vector<double>& A, const vector<double>& B, double r);

	/**
	 * This method solves the r* in calculating the swaption price
	 * @param c the coupon payment of each period
	 * @param A the vector of A value for each period
	 * @param B the vector of B value for each period
	 * @return the solution of r*
	 */
	double solveR(const vector<double>& c, const vector<double>& A, const vector<double>& B);

	/**
	 * This method calculates the value of theta for each period
	 */
	vector<double> theta();

	/**
	 * This method calculates the first or second derivate of a function numerically with cubic spline smoothing method
	 * @param value the vector of value for the function for each period
	 * @param point where to take derivatives
	 * @param n specifies the 1st or 2nd derivative
	 * @return the value of the derivative
	 */
	double takeDev(const vector<double>& value, double point, int n);

	/**
	 * This method calculates the variance of the bond ratio, which is used as the approximated
	 * value for implied volatility in the estimation of mean reversion using Method 2
	 * @param T0 the maturity of the swaption
	 * @param Tn the tenor of the swaption
	 * @return the variance swap
	 */
	double calculateVswap(double T0, double Tn);

	/**
	 * This method assigns constant mean reversion across all the data points.
	 */
	void assignConstantMeanReversion();

	/**
	 * This method assigns constant volatility across all the data points.
	 */
	void assignConstantVol();

	/**
	 * This is the logistic function used to generate the mean reversion at time ti
	 * @param A0 the short term value
	 * @param A1 the mid/long term value
	 * @param A2 the slope at the transition
	 * @param A3 the time at the transition
	 * @param ti the piecewise time interval
	 * @return the mean reversion from logistic function
	 */
	double meanReversionLogisticFunction(double A0, double A1, double A2, double A3, double ti);

	/**
	 * This function calculates the variance ratios with same maturity but different
	 * tenors
	 * @param maturity the maturity timeperiod
	 * @param tenor1 the tenor1 timeperiod
	 * @param tenor2 the tenor2 timeperiod
	 * @return the implied variance ratios
	 */
	double impliedVarianceRatios(double maturity, double tenor1, double tenor2);

	/**
	 * This method initializes and assigns constant weights to volatility grid
	 */
	void initializeAndAssignConstantWeights();

	void meanReversionCalibrationFunctionF();

	double strikeRateForSwaptionATM(double maturity, double tenor);
};

#endif /* BASEMODULE_H_ */
