/*
 * BaseModule.cpp
 *
 *  Created on: Nov 11, 2015
 *      Author: Chandra
 */
#include <cassert>
#include <iostream>
#include <math.h>
#include <random>
#include <time.h>
#include <chrono>
//#include <fstream>
#include "BaseModule.h"
#include "spline.h"
#include "Constants.h"
#include "log.h"
#include <time.h>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
  using boost::math::normal; // typedef provides default type is double.
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
  using boost::variate_generator;
#include <boost/random.hpp>
#include "BlacksFormula.h"
#include <boost/math/distributions/inverse_gaussian.hpp>

BaseModule::BaseModule() {
	// TODO Auto-generated constructor stub
}

BaseModule::~BaseModule() {
	// TODO Auto-generated destructor stub
}


void BaseModule::calculateEt(){
	data.Et.assign(data.time.size(), 0);
	data.Et[0] = exp(data.aMeanReversion[0] * data.time[0]);
	for(unsigned int i = 1; i<data.time.size(); ++i){
		data.Et[i]  = exp(data.aMeanReversion[i] * (data.time[i] - data.time[i-1])) * data.Et[i-1];
	}
}

double BaseModule::calculateB(double t, double T){
    int position_t = locate(t);
    int position_T = locate(T);
    double summation = 0;
    for (int i = position_t+1; i <= position_T; i++){
        summation = summation + (data.time[i] - data.time[i-1]) / data.Et[i];
    }
    if(position_t){
        summation = summation * data.Et[position_t];
    }
    FILE_LOG(logDEBUG3) << "calculateB = " << summation <<"\tposition_t=\t"<<position_t << "\tposition_T=\t" << position_T ;
    return summation;
}


int BaseModule::locate(double t){
	if (t==0)
		return -1;
	std::map<double,int>::iterator it;
	it = timePos.find(t);
	assert(it!=timePos.end());
    return it->second;
}

double BaseModule::calculateVr(double t){
    //Vr always have s=0 in our Problem
    if (t == 0) {
        return 0;
    }
    int position_t = locate(t);
    double Vr = 0;
    for(int i = 0; i<=position_t; ++i){
        Vr += pow(data.Et[i], 2) * pow(data.sigma[i], 2) * (data.time[i] - data.time[i-1]);
    }
    Vr = Vr / pow(data.Et[position_t], 2);
    FILE_LOG(logDEBUG3) << "calculateVr = \t" << Vr ;
    return Vr;
}

double BaseModule::calculateVp(double Tf, double Tp){
    //Vp always have t=0 in our Problem
	double t = calculateVr(Tf) * pow(calculateB(Tf, Tp), 2);
	FILE_LOG(logDEBUG3) << "calculateVp = " << t ;
    return t;
}

double BaseModule::calculateA(double t, double T){
    int position_t = locate(t);
    int position_T = locate(T);
    double P0t=1;
    if(t){
        P0t=data.priceD[position_t];
    }
    double valB = calculateB(t, T);
    double x = log(data.priceD[position_T] / P0t) ;
    double y = valB * data.forward[position_t+1];
    double z = 0.5 * pow(valB, 2) * calculateVr(t);
    double t1 = x + y - z;
    FILE_LOG(logDEBUG3) << "calculateA = " << x << "+" << y << "-" << z << "=" << t1 ;
    return t1;
}

// Returns the erf() of a value (not super precise, but ok)
double BaseModule::erf(double x){
    double y = 1.0 / ( 1.0 + 0.3275911 * x);
    return 1 - (((((
                    + 1.061405429  * y
                    - 1.453152027) * y
                   + 1.421413741) * y
                  - 0.284496736) * y
                 + 0.254829592) * y)
    * exp (-x * x);
}

// Returns the probability of [-inf,x] of a gaussian distribution
double BaseModule::N(double x, double mu, double sigma){
    return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2.))));
}

double BaseModule::ZBP(double Tf, double Tp, double X){
    int position_Tf = locate(Tf);
    int position_Tp = locate(Tp);
    double l = log(data.priceD[position_Tf] * X / data.priceD[position_Tp]);
    double v = sqrt(calculateVp(Tf, Tp));
    double dP =  l / v + 0.5 * v;
    double dN = l / v - 0.5 * v;
    return X * data.priceD[position_Tf] * N(dP, 0, 1) - data.priceD[position_Tp] * N(dN, 0, 1);
}

//f(r)
double BaseModule::fun_r(const vector<double>& c, const vector<double>& A, const vector<double>& B, double r){
    double ret = 0;
    double t1=0, t2=0;
    for(unsigned int i = 0; i < c.size(); ++i){
    	t1=B[i] * r;
    	t2 = exp(A[i] - t1);
        ret += c[i] * t2;
        FILE_LOG(logDEBUG3) << "fun_r\t" << i << "\t" << "\t" << A[i]<< "\t" << B[i] <<"\t" << c[i] << "\t" << "\t"<<t1<<"\t"<<t2<<"\t"<<ret;
    }
    return ret - 1;
}

//f'(r)
double BaseModule::fun_dev(const vector<double>& c, const vector<double>& A, const vector<double>& B, double r){
    double d = 0;
    double t3 = 0, t1=0, t2=0 ;
    for(unsigned int i = 0; i < c.size(); ++i){
    	t1 = B[i] * r;
    	t2 = A[i] - t1;
    	t3 = exp(t2) * B[i];
        d = d - c[i] * t3;
        FILE_LOG(logDEBUG3) << "fun_dev\t" << i << "\t" << A[i]<< "\t" << B[i] <<"\t" << c[i] << "\t" << t1<<"\t"<<t2<<"\t"<<t3<<"\t"<<d;
    }
    return d;
}

//Newton's algorithm to solve r*
double BaseModule::solveR(const vector<double>& c, const vector<double>& A, const vector<double>& B){
    double r = 0.5;
    double fun_rValue = fun_r(c, A, B, r);
    double fun_devValue = 0;
    double temp=0;
    FILE_LOG(logDEBUG3) << "solveR-->\t" << fun_rValue << "\t" << fun_devValue <<"\t" << temp <<"\t"<<r;
    while(fabs(fun_rValue) > 1E-10){
    	fun_devValue = fun_dev(c, A, B, r);
    	temp = (fun_rValue/fun_devValue);
    	r = r-temp;
    	FILE_LOG(logDEBUG3) << "solveR\t" << fun_rValue << "\t" << fun_devValue <<"\t" <<temp <<"\t"<<r;
    	fun_rValue = fun_r(c, A, B, r);
    }
    return r;
}

double BaseModule::pSwaption(double strike, double T0, double Tp){
    int position_T0 = locate(T0);
    int position_Tp = locate(Tp);
    int cashFlowTimes = position_Tp-position_T0;
    vector<double> ci(cashFlowTimes);
    vector<double> Ti(cashFlowTimes);
    vector<double> Ai(cashFlowTimes);
    vector<double> Bi(cashFlowTimes);
    for(int i = 0; i < cashFlowTimes; ++i){
        Ti[i] = data.time[position_T0 + i + 1];
        ci[i] = strike * serviceLocator.getDefaultsPaymentFrequency();
        Ai[i]=calculateA(T0, Ti[i]);
        Bi[i]=calculateB(T0, Ti[i]);
        FILE_LOG(logDEBUG) << "Ai=" << Ai[i] << " Bi=" << Bi[i] << " ci=" << ci[i] << " Ti=" << Ti[i]
			 << " time=" << data.time[position_T0 + i] << " strike=" << strike;
    }
    ci.back() += 1;
    double r_star = solveR(ci, Ai, Bi);
    vector<double> Xi(cashFlowTimes);
    double price = 0;
    double zbprice = 0.0;
    for(int i = 0; i < cashFlowTimes; ++i){
    	zbprice = ZBP(T0, Ti[i], Xi[i]);
    	FILE_LOG(logDEBUG) << "Ai=" << Ai[i] << " Bi=" << Bi[i] << " ci=" << ci[i] << " r_star=" << r_star << " ZBP=" << zbprice ;
        Xi[i] = exp(Ai[i] - Bi[i] * r_star);
        price += ci[i] * ZBP(T0, Ti[i], Xi[i]);
        FILE_LOG(logDEBUG) << " Xi=" << Xi[i] << " price=" << price;
    }
    FILE_LOG(logDEBUG) << "FinalPrice=" <<"\t"<< price << "\t "<< T0 << " \t" << Tp ;
    return price;
}


double BaseModule::takeDev(const vector<double>& value, double point, int n){
    //take first or second derivative
    tk::spline s;
    s.set_points(data.time,value);
    if (n==1) {
        return (s(point+0.01)-s(point-0.01))/0.02;
    }
    else{
        return (s(point+0.01)-2*s(point)+s(point-0.01))/0.0001;
    }
}

vector<double> BaseModule::theta(){
    vector<double> V(data.time.size());
    vector<double> theta_vec(data.time.size());
    for(unsigned int i = 0; i < data.time.size(); ++i){
        vector<double> sig(i+1);
        for(unsigned int j=0; j<=i; j++){
            sig[j]=data.sigma[j] * calculateB(data.time[j], data.time[i]);
        }
        V[i]=pow(sig[0], 2) * data.time[0] ;
        for(unsigned int j = 1; j <= i; ++j){
            V[i]= V[i]+pow(sig[j], 2) * (data.time[j]-data.time[j-1]);
        }
    }
    for(unsigned int i = 0; i < data.time.size(); ++i){
        theta_vec[i]=takeDev(data.forward, data.time[i], 1)
        		+ data.aMeanReversion[i] * data.forward[i]
        		+ 0.5 * (takeDev(V, data.time[i], 2) + data.aMeanReversion[i] * takeDev(V, data.time[i], 1));
    }
    return theta_vec;
}

double BaseModule::calculateVswap(double T0, double Tn){
	double value=0.0;
	double Vp = calculateVp(T0, Tn);
	double priceT0 = data.priceD[locate(T0)];
	double priceTn = data.priceD[locate(Tn)];
	double val2 = pow((priceT0/(priceT0 - priceTn)),2.0);
	value = Vp * val2;
	FILE_LOG(logDEBUG) << "Vswap=\t" << value << "\tPriceT0=\t" << priceT0 << "\tPriceTn=\t" <<
			priceTn << "\tVp=\t" << Vp << "\tVswap=\t" << value;
	return value;
}


void BaseModule::assignConstantMeanReversion(double val){
	int x = data.time.size();
	data.aMeanReversion.assign(x, val);
}

void BaseModule::assignConstantVol(double val){
	int x = data.time.size();
	data.sigma.assign(x, val);
}

double BaseModule::meanReversionLogisticFunction(double A0, double A1, double A2, double A3, double ti){
	return (A0 + (A1-A0)/(1+exp(A2*(A3-ti))));
}

double BaseModule::impliedVarianceRatios(double maturity, double tenor1, double tenor2){
	int maturityPos = locate(maturity);
	int tenorPos1 = locate(tenor1);
	int tenorPos2 = locate(tenor2);
	double pM = data.priceD[maturityPos];
	double pT1 = data.priceD[tenorPos1];
	double pT2 = data.priceD[tenorPos2];
	double b1 = calculateB(maturity, tenor1);
	double b2 = calculateB(maturity, tenor2);
	double num = (pM-pT2)*b1;
	double den = (pM-pT1)*b2;
	double value = pow((num/den), 2.0);
	FILE_LOG(logDEBUG3) << "ImplVarianceRatio=\t" <<"Mat=\t"<<maturity<<"\ttenor1=\t"<<tenor1<<"\ttenor2=\t"
			<<tenor2<< "\tpM=\t" <<pM<<"\tpT1=\t" <<pT1 <<"\tpT2=\t"
			<<pT2<< "\tb1=\t" <<b1 <<"\tb2="<< b2<<"\tnum=\t" << num << "\tden=\t" << den << "\t" << value;
	return value;
}

void BaseModule::initializeAndAssignConstantWeights(){
	int rows = actualVol.maturity.size();
	int cols = actualVol.tenor.size();
	for (int i=0; i<rows; i++){
		vector<double> rowi (cols, serviceLocator.getDefaultsPaymentFrequency());
		actualVol.weights.push_back(rowi);
	}
}

double BaseModule::meanReversionCalibrationFunctionF(){
	int rows = actualVol.maturity.size();
	int cols = actualVol.tenor.size();
	vector<vector<double> > func;
	double funcTotal = 0.0;
	for (int i=0; i<rows; i++){
		vector<double> colsMinus1(cols-1);
		func.push_back(colsMinus1);
	}

	for (int i=0; i<rows; i++){
		for (int j=0; j<cols-1; j++){
			if (actualVol.weights[i][j] != 0.0){
				double tenor1 = actualVol.tenor[j]+actualVol.maturity[i];
				double tenor2 = actualVol.tenor[j+1]+actualVol.maturity[i];
				double weightRatio = actualVol.weights[i][j+1]/actualVol.weights[i][j];
				double volRatio = sqrt(impliedVarianceRatios(actualVol.maturity[i], tenor1, tenor2));
				double impliedVolRatio = actualVol.vol[i][j+1]/actualVol.vol[i][j];
				func[i][j] = weightRatio * pow((volRatio-impliedVolRatio),2.0);
				funcTotal += func[i][j];
				FILE_LOG(logDEBUG3) << "MeanReversionM2[" << i << "][" << j <<"]" <<"\t" << "Mat" << "\t" << actualVol.maturity[i] << "\t"
						<< "Ten1" << "\t" << actualVol.tenor[j] << "\t"  << "Ten2" << "\t" << actualVol.tenor[j+1] << "\t"
						<< "weightRatio=\t" << weightRatio << "\tVswapRatio=\t" << volRatio << "\tImplVolRatio=\t"
						<< impliedVolRatio << "\tFunc=\t" << func[i][j] << "\tFuncTotal=\t" << funcTotal;
			}
		}
	}
	return funcTotal;
}

double BaseModule::volatilityCalibrationFunctionG(){
	FILE_LOG(logDEBUG3) << "VolatilityCalibration-FunctionG starts";
	int rows = actualVol.maturity.size();
	int cols = actualVol.tenor.size();
	vector<vector<double> > func;
	double funcTotal = 0.0;
	for (int i=0; i<rows; i++){
		vector<double> colsMinus1(cols);
		func.push_back(colsMinus1);
	}

	BlacksFormula black;
	vector<vector<double> > blckPrices;

	for (int i=0; i<rows; i++){
		for (int j=0; j<cols; j++){
			if (actualVol.weights[i][j] != 0.0){
				double t = actualVol.maturity[i];
				double T = actualVol.tenor[j]+actualVol.maturity[i];
				double strike = actualVol.strikeRate[i][j];

				double priceModel = pSwaption(strike, t, T);
				double priceBlack = blackSwaptionPriceATM(t, actualVol.tenor[j], actualVol.vol[i][j], actualVol.strikeRate[i][j], true);
				double priceRatio = priceModel/priceBlack;
				func[i][j] = actualVol.weights[i][j]*pow((priceRatio - 1),2.0);
				funcTotal += func[i][j];
				FILE_LOG(logDEBUG3) << "VolatilityCalibration[" << i << "][" << j <<"]" <<"\t" << "Mat"
						<< "\t" << actualVol.maturity[i] << "\t" << "Ten" << "\t" << actualVol.tenor[j] << "\t"
						<< "PriceModel=\t" << priceModel << "\tPriceBlack=\t" << priceBlack
						<< "\tweight=\t" << actualVol.weights[i][j] << "\tPriceRatio=\t" << priceRatio
						<< "\tFunc=\t" << func[i][j] << "\tFuncTotal=\t" << funcTotal;
			}
		}
	}
	return funcTotal;
}

double BaseModule::strikeRateForSwaptionATM(double maturity, double tenor){
	int matPosition = locate(maturity);
	int tenPosition = locate(tenor);
	double cashFlowsSum = 0.0;
	for(int i=matPosition+1; i<=tenPosition; i++){
		cashFlowsSum += serviceLocator.getDefaultsPaymentFrequency()*data.priceD[i];
	}
	return (data.priceD[matPosition]-data.priceD[tenPosition])/cashFlowsSum;
}

void BaseModule::initializeStrikeRateForSwaptionATM(){
	for(vector<double>::iterator itMaturity = actualVol.maturity.begin(); itMaturity != actualVol.maturity.end(); ++itMaturity){
		vector<double> strikeRatesInner;
		for(vector<double>::iterator itTenor = actualVol.tenor.begin(); itTenor != actualVol.tenor.end(); ++itTenor){
			double strikeRate = strikeRateForSwaptionATM(*itMaturity, *itMaturity+*itTenor);
			strikeRatesInner.push_back(strikeRate);
		}
		actualVol.strikeRate.push_back(strikeRatesInner);
	}
}

vector<double> BaseModule::simulatedAnnealingFuncForMeanReversion(){
	double A0_min = serviceLocator.getMeanReversionA0(), A1_max = serviceLocator.getMeanReversionA1();
	double A0_0 = A0_min * 1.2, A1_0 = A1_max*0.8, A2_0 = serviceLocator.getMeanReversionA2(); //equivalent to x=A0, A1, A2
	double sd0_0= serviceLocator.getMeanReversionSDA0(),
			sd1_0=serviceLocator.getMeanReversionSDA1(),
			sd2_0=serviceLocator.getMeanReversionSDA2();

	assignVaryingMeanReversion(A0_0, A1_0, A2_0);
	double fLast = meanReversionCalibrationFunctionF();

	double gamma = serviceLocator.getMeanReversionGamma();
	int totalNoOfSimulation = serviceLocator.getMeanReversionSimulations();
	int iterationNo = 1;
	vector<double> fxVector, xrVector, sigmajVector;
	fxVector.push_back(fLast);
	double fx = fLast;
	double A0 = A0_0, A1=A1_0, A2=A2_0;
	int counter = 0;
	double minFx = fLast;
	double minA0 = A0_0;
	double minA1 = A1_0;
	double minA2 = A2_0;

	boost::mt19937 *rng2 = new boost::mt19937();
	rng2->seed(time(NULL));

	FILE_LOG(logDEBUG) << "SimulatedAnnMeanReversion\t" << iterationNo << "\t" << A0_0 << "\t" << sd0_0 << "\t" << A1_0 << "\t" << sd1_0 << "\t" << A2_0 << "\t" << sd2_0 << "\t" << 0 << "\t" << fLast;
	do {
		double sd0j = coolingMechanism(gamma, sd0_0, totalNoOfSimulation, iterationNo);
		normal_distribution<> distribution0(A0, sd0j);
		variate_generator<boost::mt19937&, normal_distribution<> > generate_next0(*rng2, distribution0);

		double sd1j = coolingMechanism(gamma, sd1_0, totalNoOfSimulation, iterationNo);
		normal_distribution<> distribution1(A1, sd1j);
		variate_generator<boost::mt19937&, normal_distribution<> > generate_next1(*rng2, distribution1);

		double sd2j = coolingMechanism(gamma, sd2_0, totalNoOfSimulation, iterationNo);
		normal_distribution<> distribution2(A2, sd2j);
		variate_generator<boost::mt19937&, normal_distribution<> > generate_next2(*rng2, distribution2);

		counter = 0;
		int simulationInnerCount = coolingMechanism(gamma, totalNoOfSimulation, totalNoOfSimulation, iterationNo);
		do{
			A0 = generate_next0();
			A1 = generate_next1();
			A2 = generate_next2();
			if((A0_min<A0) && (A0<A1_max) && (A0_min<A1) && (A1<A1_max) && (A2>0) && (A0<A1)){
				assignVaryingMeanReversion(A0, A1, A2);
				fx = meanReversionCalibrationFunctionF();
				FILE_LOG(logDEBUG) << "SimulatedAnnMeanReversion-Counter\t" << iterationNo << "\t" << A0 << "\t" << sd0j
						<< "\t" << A1 << "\t" << sd1j << "\t" << A2 << "\t" << sd2j << "\t" << counter << "\t" << fx;
				++counter;
			}
		}while(fx>=fLast && counter<simulationInnerCount);
		fxVector.push_back(fx);
		fLast = fx;
		++iterationNo;
		if (minFx > fx){
			minFx = fx;
			minA0 = A0;
			minA1 = A1;
			minA2 = A2;
		}

		FILE_LOG(logDEBUG) << "SimulatedAnnMeanReversion\t" << iterationNo << "\t" << A0 << "\t" << sd0j << "\t" << A1
				<< "\t" << sd1j << "\t" << A2 << "\t" << sd2j << "\t" << counter << "\t" << fx
				<< "\t" << minA0 << "\t" << minA1 << "\t" << minA2;
	}while(iterationNo < totalNoOfSimulation);
	FILE_LOG(logINFO) << "SimulatedAnnMeanReversion-Final\t" << minFx << "\t" << minA0 << "\t" << minA1 << "\t" << minA2;
	vector<double> retValue;
	retValue.push_back(minA0);
	retValue.push_back(minA1);
	retValue.push_back(minA2);
	retValue.push_back(minFx);
	return retValue;
}

void BaseModule::assignVaryingMeanReversion(double A0, double A1, double A2, double A3){
	FILE_LOG(logDEBUG3) << "VaryingMeanRev\tA3=\t"<<A3;
	for(vector<double>::iterator it = data.time.begin(), it1 = data.aMeanReversion.begin(); it!=data.time.end(); ++it, ++it1){
		*it1 = logisticFunc(A0, A1, A2, A3, *it);
	}
	for(vector<double>::iterator it = data.time.begin(), it1 = data.aMeanReversion.begin(); it!=data.time.end(); ++it, ++it1){
		FILE_LOG(logDEBUG3) << "VaryingMeanRev\t" << *it << "\t" << *it1 ;
	}
	calculateEt();
}

void BaseModule::assignVaryingMeanReversion(double A0, double A1, double A2){
	double A3 = (*(--(actualVol.maturity.end())) + *(--(actualVol.tenor.end())))/2;
	assignVaryingMeanReversion(A0, A1, A2, A3);
}

double BaseModule::cubicFunc(double a0, double a1, double a2, double a3, double ti){
	return a0+a1*ti+a2*pow(ti,2)+a3*pow(ti,3);
}

void BaseModule::assignVaryingVolatility(double a0, double a1, double a2, double a3){
	double x1 = (-a2+sqrt(pow(a2,2) - 3*a1*a3))/(3*a3);
	double x2 = -a2 / (3*a3);
	FILE_LOG(logDEBUG3) << "VaryingSigma-(x1,x2)\t" << x1 << "\t" << x2;
	for(vector<double>::iterator it = data.time.begin(), it1 = data.sigma.begin(); it!=data.time.end(); ++it, ++it1){
		if(*it < x2){
			*it1 = (a1 + 2 * a2 * x2 + 3 * a3 * pow(x2, 2)) * (*it - x2) + cubicFunc(a0, a1, a2, a3, x2);
		}
		else if(*it > x1){
			*it1 = cubicFunc(a0, a1, a2, a3, x1);
		}
		else{
			*it1 = cubicFunc(a0, a1, a2, a3, *it);
		}
	}
	for(vector<double>::iterator it = data.time.begin(), it1 = data.sigma.begin(); it!=data.time.end(); ++it, ++it1){
		FILE_LOG(logDEBUG3) << "VaryingSigma\t" << *it << "\t" << *it1 ;
	}
}

void BaseModule::assignVaryingVolatilityUpwardSloping(double a0, double a1, double a2, double a3){
	double x1 = (-a2-sqrt(pow(a2,2) - 3*a1*a3))/(3*a3);
	double x2 = -a2 / (3*a3);
	FILE_LOG(logDEBUG3) << "VaryingSigma-(x1,x2)\t" << x1 << "\t" << x2;
	for(vector<double>::iterator it = data.time.begin(), it1 = data.sigma.begin(); it!=data.time.end(); ++it, ++it1){
		if(*it < x2){
			*it1 = (a1 + 2 * a2 * x2 + 3 * a3 * pow(x2, 2)) * (*it - x2) + cubicFunc(a0, a1, a2, a3, x2);
		}
		else if(*it > x1){
			*it1 = cubicFunc(a0, a1, a2, a3, x1);
		}
		else{
			*it1 = cubicFunc(a0, a1, a2, a3, *it);
		}
	}
	for(vector<double>::iterator it = data.time.begin(), it1 = data.sigma.begin(); it!=data.time.end(); ++it, ++it1){
		FILE_LOG(logDEBUG3) << "VaryingSigma\t" << *it << "\t" << *it1 ;
	}
}

bool BaseModule::validCubicFunc(double a0, double a1, double a2, double a3){
	if(a3<=0){
		//assure the upward shape
		return false;
	}
	if(pow(a2, 2) - 3*a1*a3 <= 0 ){
		//assure two solutions of f'=0
		return false;
	}
	double x1 = (-a2+sqrt(pow(a2,2) - 3*a1*a3))/(3*a3);
	double x2 = -a2 / (3*a3);
	if(cubicFunc(a0,a1,a2,a3,x1)<0){
		//assure positive volatility
		return false;
	}
	if(a2>=0){
		//assure x2>0
		return false;
	}
//	if(x1>=*data.time.end()){
	if(x1>= *(--actualVol.maturity.end()) + *(--actualVol.tenor.end())){
		//assure x1>maximum maturity
		return false;
	}
	double fx2 = 2.0/100.0;
	if((a1 + 2 * a2 * x2 + 3 * a3 * pow(x2, 2)) * (0 - x2) + cubicFunc(a0, a1, a2, a3, x2)>=fx2){
		return false;
	}
	return true;
}

bool BaseModule::validCubicFuncUpwardSloping(double a0, double a1, double a2, double a3){
	if(a3>=0){
		//assure the upward shape
		return false;
	}
	if(pow(a2, 2) - 3*a1*a3 <= 0 ){
		//assure two solutions of f'=0
		return false;
	}
	double x1 = (-a2-sqrt(pow(a2,2) - 3*a1*a3))/(3*a3);
	double x2 = -a2 / (3*a3);

	if((a1 + 2 * a2 * x2 + 3 * a3 * pow(x2, 2)) * (0.25 - x2) + cubicFunc(a0, a1, a2, a3, x2)<0){
		//assure positive volatility
		return false;
	}
	if(a2<=0){
		//assure x2>0
		return false;
	}
//	if(x1>=*data.time.end()){
	if(x1>= *(--actualVol.maturity.end()) + *(--actualVol.tenor.end())){
		//assure x1>maximum maturity
		return false;
	}
	double fx2 = 2.0/100.0;
	if(cubicFunc(a0,a1,a2,a3,x1)>=fx2){
		return false;
	}
	return true;
}

double BaseModule::logisticFunc(double A0, double A1, double A2, double A3, double ti){
	return A0+(A1-A0)/(1+exp(A2*(A3-ti)));
}

double BaseModule::coolingMechanism(double gamma, double sigma0, int totalNoOfSimulation, int iterationNo){
	return sigma0 * exp(-gamma*iterationNo/(totalNoOfSimulation-1));
}

vector<vector<double> > BaseModule::transposeVector(vector<vector<double> > numbers){
	vector<vector<double> >::iterator vvi_iterator;
	vector<double>::iterator vi_iterator;

	int columns = numbers.front().size();
	vector<vector<double> > newNumbers;
	for(int i = 0; i<columns; i++){
		vector<double> newNumbersCols;
		int j;
		for(vvi_iterator = numbers.begin(); vvi_iterator!=numbers.end(); ++vvi_iterator) {
			for(vi_iterator = (*vvi_iterator).begin(), j=0; j<i; ++vi_iterator,++j){}
			newNumbersCols.push_back(*vi_iterator);
		}
		newNumbers.push_back(newNumbersCols);
	}
	return newNumbers;
}


void BaseModule::printVector(vector<double> vec){
	for(vector<double>::iterator iterator = vec.begin(); iterator!=vec.end(); ++iterator){
		cout << *iterator << "\t";
	}
	cout << endl;
}
void BaseModule::printVectorVector(vector<vector<double> > vec){
	for(vector<vector<double> >::iterator iterator = vec.begin(); iterator!=vec.end(); ++iterator){
		for(vector<double>::iterator iterator2 = (*iterator).begin(); iterator2!=(*iterator).end(); ++iterator2){
			cout << *iterator2 << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

void BaseModule::simulatedAnnealingFuncForVolatility(){
	double A0_0=volatilityParams.A0, A1_0=volatilityParams.A1, A2_0=volatilityParams.A2,A3_0=volatilityParams.A3;
	double sd0_0=fabs(A0_0)*serviceLocator.getVolatilitySDA0();
	double sd1_0=fabs(A1_0)*serviceLocator.getVolatilitySDA1();
	double sd2_0=fabs(A2_0)*serviceLocator.getVolatilitySDA2();
	double sd3_0=fabs(A3_0)*serviceLocator.getVolatilitySDA3();

	FILE_LOG(logDEBUG) << "SimulatedAnnVol-Local-1stCall";
	assignVaryingMeanReversion(serviceLocator.getMeanReversionA0(), serviceLocator.getMeanReversionA1(), serviceLocator.getMeanReversionA2());

	calculateEt();
	assignVaryingVolatility(A0_0, A1_0, A2_0, A3_0);
//	assignVaryingVolatilitySpline(A0_0, A1_0, A2_0, A3_0);
	double gLast = volatilityCalibrationFunctionG();

	double gamma = serviceLocator.getVolatilityGamma();
	int totalNoOfSimulation = serviceLocator.getVolatilitySimulations();
	int iterationNo = 1;
	vector<double> fxVector;//, xrVector, sigmajVector;
	fxVector.push_back(gLast);
	double gx = gLast;
	double A0 = A0_0, A1=A1_0, A2=A2_0, A3=A3_0;
	int counter = 0;
	double minGx = gLast;

	boost::mt19937 *rng2 = new boost::mt19937();
	rng2->seed(time(NULL));

	FILE_LOG(logINFO) << "SimulatedAnnVol-Local\t" << iterationNo << "\t" << A0 << "\t" << sd0_0 << "\t"
			<< A1 << "\t" << sd1_0 << "\t" << A2 << "\t" << sd2_0 << "\t" << A3 << "\t" << sd3_0 << "\t" << counter << "\t" << gx
			<< "\tGlobal\t" << volatilityParams.Gx << "\t" << volatilityParams.A0
			<< "\t" << volatilityParams.A1 << "\t" << volatilityParams.A2 << "\t" << volatilityParams.A3;
	do {
		double sd0j = coolingMechanism(gamma, sd0_0, totalNoOfSimulation, iterationNo);
		normal_distribution<> distribution0(A0, sd0j);
		variate_generator<boost::mt19937&, normal_distribution<> > generate_next0(*rng2, distribution0);

		double sd1j = coolingMechanism(gamma, sd1_0, totalNoOfSimulation, iterationNo);
		normal_distribution<> distribution1(A1, sd1j);
		variate_generator<boost::mt19937&, normal_distribution<> > generate_next1(*rng2, distribution1);

		double sd2j = coolingMechanism(gamma, sd2_0, totalNoOfSimulation, iterationNo);
		normal_distribution<> distribution2(A2, sd2j);
		variate_generator<boost::mt19937&, normal_distribution<> > generate_next2(*rng2, distribution2);

		double sd3j = coolingMechanism(gamma, sd3_0, totalNoOfSimulation, iterationNo);
				normal_distribution<> distribution3(A3, sd3j);
				variate_generator<boost::mt19937&, normal_distribution<> > generate_next3(*rng2, distribution3);

		counter=0;
		int simulationInnerCount = coolingMechanism(gamma, totalNoOfSimulation, totalNoOfSimulation, iterationNo);
		do{
			A0 = generate_next0();
			A1 = generate_next1();
			A2 = generate_next2();
			A3 = generate_next3();
			if(validCubicFunc(A0,A1,A2,A3)){
//			if(validCubicFuncSpline(A0,A1,A2,A3)){
				assignVaryingVolatility(A0, A1, A2, A3);
//				assignVaryingVolatilitySpline(A0, A1, A2, A3);
				gx = volatilityCalibrationFunctionG();
				FILE_LOG(logDEBUG) << "SimulatedAnnVol-Counter\t" << iterationNo << "\t" << counter << "\t"
						<< A0 << "\t" << sd0j << "\t" << A1
						<< "\t" << sd1j << "\t" << A2 << "\t" << sd2j << "\t" << A3
						<< "\t" << sd3j<< "\t" << gx;
				++counter;
			}
		}while(gx>=gLast && counter<=simulationInnerCount);
		fxVector.push_back(gx);
		gLast = gx;
		++iterationNo;
		if (minGx > gx){
			volatilityParams.A0 = A0;
			volatilityParams.A1 = A1;
			volatilityParams.A2 = A2;
			volatilityParams.A3 = A3;
			volatilityParams.Gx = gx;
		}

		FILE_LOG(logDEBUG) << "SimulatedAnnVol-Local\t" << iterationNo << "\t" << A0 << "\t" << sd0j << "\t"
				<< A1 << "\t" << sd1j << "\t" << A2 << "\t" << sd2j << "\t" << A3 << "\t" << sd3j << "\t" << counter << "\t" << gx
				<< "\tGlobal\t" << volatilityParams.Gx << "\t" << volatilityParams.A0
				<< "\t" << volatilityParams.A1 << "\t" << volatilityParams.A2 << "\t" << volatilityParams.A3;
	}while(iterationNo < totalNoOfSimulation);
	vector<double> retValue;
	FILE_LOG(logINFO) << "SimulatedAnnVol-GlobalMin-Gx\t" << volatilityParams.Gx << "\t" << volatilityParams.A0
			<< "\t" << volatilityParams.A1 << "\t" << volatilityParams.A2 << "\t" << volatilityParams.A3;
}

void BaseModule::checkMeanReversionConvergence(int loopCount){
	FILE_LOG(logDEBUG3) << "ConvergeMeanReversion\tCounter\tA0\tA1\tA2\tFx";
	for(int i=0; i<loopCount; ++i){
		vector<double> convergeVector = simulatedAnnealingFuncForMeanReversion();
		double A0 = convergeVector[0];
		double A1 = convergeVector[1];
		double A2 = convergeVector[2];
		double Fx = convergeVector[3];
		FILE_LOG(logDEBUG3) << "ConvergeMeanReversion\t" <<i+1<<"\t"<< A0 << "\t" << A1 <<"\t" << A2 << "\t" << Fx;
	}
}

double BaseModule::blackSwaptionPriceATM(double maturity, double tenor, double implyVol, double swapRate, bool isPayer){
    double w=-1;
    if(isPayer){
        w=1;
    }
    double d1=implyVol * sqrt(maturity)/2;
    double d2=-d1;
    double Bl = swapRate * w * (N(w*d1,0,1) - N(w*d2,0,1));
    double multiplier = 0;
    for(int i = locate(maturity)+1; i <= locate(maturity + tenor);i++)
        multiplier += serviceLocator.getDefaultsConstantWeights() * data.priceD[i];
    return Bl * multiplier;
}

double BaseModule::blackImpliedVolatilitySwaptionATM(double maturity, double tenor, double price, double swapRate){
	double multiplier = 0;
	for(int i = locate(maturity)+1; i <= locate(maturity + tenor);i++){
		multiplier += serviceLocator.getDefaultsConstantWeights() * data.priceD[i];
	}
	double val = ((price/multiplier)+swapRate)/2.0/swapRate;
	return RationalApproximation(val)*2;
}

double BaseModule::RationalApproximation(double p)
{
	boost::math::normal dist(0.0, 1.0);
	return boost::math::quantile(dist, p);
}

void BaseModule::assignVaryingVolatilitySpline(double a0, double a1, double a2, double a3){
	vector<double> value{a0,a1,a2,a3};
	vector<double> t{0.25,7,14,20};
	tk::spline s;
	s.set_points(t,value);
	double tme = actualVol.maturity.back() + actualVol.tenor.back();
	for(vector<double>::iterator it = data.time.begin(), it1 = data.sigma.begin(); *it<=tme; ++it, ++it1){
		*it1=s(*it);
	}
}

bool BaseModule::validCubicFuncSpline(double a0, double a1, double a2, double a3){
	double fx2=0.12;
	if(a0>fx2 || a1>fx2 || a2>fx2 || a3>fx2 || a0<=0 || a1<=0 || a2<=0 || a3<=0)
		return false;
	return true;
}
