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

class BaseModule {
public:
	struct DataStruct{
		vector<double> t;
		vector<double> f;
		vector<double> p;
		vector<double> a;
		vector<double> Et;
		vector<double> A;
		vector<double> B;
		vector<double> sigma;
	};
	struct ImpliedVolStruct{
		//Vector of Rows(for one Maturity) of Columns(for one tenor)
		//Inner vector holds one row for one maturity, i.e. one maturity multiple tenors
		//Outer vector holds multiple rows for different maturity vectors
		vector<vector<double> > vol;
		vector<double> tenor;
		vector<double> maturity;
	};
	BaseModule();
	virtual ~BaseModule();

	virtual int moduleMainFunc() = 0;

	DataStruct data;
	ImpliedVolStruct implVol;
	map<double, int> timePos;

	void calculateEt();
	double calculateB(double t, double T);
	int locate(double t);
	double calculateVr(double t);
	double calculateVp(double Tf, double Tp);
	double calculateA(double t, double T);

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
