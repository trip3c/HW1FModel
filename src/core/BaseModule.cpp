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


void main1(void);


BaseModule::BaseModule() {
	// TODO Auto-generated constructor stub
//	logger = Logger::getLogger();
//	logg = Logger::getLogger();
}

BaseModule::~BaseModule() {
	// TODO Auto-generated destructor stub
//	getLo.close();
//	logg.close();
}


void BaseModule::calculateEt(){
	data.Et.assign(data.time.size(), 0);
	data.Et[0] = exp(data.aMeanReversion[0] * data.time[0]);
	for(int i = 1; i<data.time.size(); ++i){
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
//    if(t){ //TODO Xiaohan to validate - Chandra
    if(position_t){
        summation = summation * data.Et[position_t];
    }
    FILE_LOG(logDEBUG) << "calculateB = " << summation <<"\tposition_t=\t"<<position_t << "\tposition_T=\t" << position_T ;
    return summation;
}


int BaseModule::locate(double t){
    //return index in vector for given value
//	std::cout << t <<" ";
	if (t==0)
		return -1;
	std::map<double,int>::iterator it;
	it = timePos.find(t);
	assert(it!=timePos.end());
//	std::cout << it->second << std::endl;
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
    FILE_LOG(logDEBUG) << "calculateVr = \t" << Vr ;
    return Vr;
}

double BaseModule::calculateVp(double Tf, double Tp){
    //Vp always have t=0 in our Problem
	double t = calculateVr(Tf) * pow(calculateB(Tf, Tp), 2);
	FILE_LOG(logDEBUG) << "calculateVp = " << t ;
    return t;
}

double BaseModule::calculateA(double t, double T){
    int position_t = locate(t);
    int position_T = locate(T);
    double P0t=1;
    if(t){
        P0t=data.priceD[position_t];
        //f0t=data.f[position_t];
    }
    double valB = calculateB(t, T);
    double x = log(data.priceD[position_T] / P0t) ;
    double y = valB * data.forward[position_t+1];
    double z = 0.5 * pow(valB, 2) * calculateVr(t);
    double t1 = x + y - z;
    FILE_LOG(logDEBUG) << "calculateA = " << x << "+" << y << "-" << z << "=" << t1 ;
    return t1;
}

//-----------------------


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
    for(int i = 0; i < c.size(); ++i){
        ret += c[i] * exp(A[i] - B[i] * r);
    }
    return ret - 1;
}

//f'(r)
double BaseModule::fun_dev(const vector<double>& c, const vector<double>& A, const vector<double>& B, double r){
    double d = 0;
    for(int i = 0; i < c.size(); ++i){
        d = d - c[i] * exp(A[i] - B[i] * r) * B[i];
    }
    return d;
}

//Newton's algorithm to solve r*
double BaseModule::solveR(const vector<double>& c, const vector<double>& A, const vector<double>& B){
    double r = 0.5;
    while (fabs(fun_r(c, A, B, r)) > 1E-10) {
        r = r - fun_r(c, A, B, r) / fun_dev(c, A, B, r);
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
        ci[i] = strike * (Ti[i] - data.time[position_T0 + i]);
        Ai[i]=calculateA(T0, Ti[i]);
        Bi[i]=calculateB(T0, Ti[i]);
        FILE_LOG(logDEBUG) << "Ai=" << Ai[i] << " Bi=" << Bi[i] << " ci=" << ci[i] << " Ti=" << Ti[i] << " time=" << data.time[position_T0 + i] << " strike=" << strike;
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
    for(int i = 0; i < data.time.size(); ++i){
        vector<double> sig(i+1);
        for(int j=0; j<=i; j++){
            sig[j]=data.sigma[j] * calculateB(data.time[j], data.time[i]);
        }
        V[i]=pow(sig[0], 2) * data.time[0] ;
        for(int j = 1; j <= i; ++j){
            V[i]= V[i]+pow(sig[j], 2) * (data.time[j]-data.time[j-1]);
        }
    }
    for(int i = 0; i < data.time.size(); ++i){
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
//	double a = constants.FIXED_MEAN_REVERSION;
	int x = data.time.size();
	data.aMeanReversion.assign(x, val);
}

void BaseModule::assignConstantVol(double val){
//	double s1 = constants.FIXED_VOLATILITY;
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
	FILE_LOG(logDEBUG) << "ImplVarianceRatio=\t" <<"Mat=\t"<<maturity<<"\ttenor1=\t"<<tenor1<<"\ttenor2=\t"
			<<tenor2<< "\tpM=\t" <<pM<<"\tpT1=\t" <<pT1 <<"\tpT2=\t"
			<<pT2<< "\tb1=\t" <<b1 <<"\tb2="<< b2<<"\tnum=\t" << num << "\tden=\t" << den << "\t" << value;
	return value;
}

void BaseModule::initializeAndAssignConstantWeights(){
	int rows = actualVol.maturity.size();
	int cols = actualVol.tenor.size();
	for (int i=0; i<rows; i++){
		vector<double> rowi (cols, constants.CONSTANT_WEIGHT);
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
				FILE_LOG(logDEBUG) << "MeanReversionM2[" << i << "][" << j <<"]" <<"\t" << "Mat" << "\t" << actualVol.maturity[i] << "\t"
						<< "Ten1" << "\t" << actualVol.tenor[j] << "\t"  << "Ten2" << "\t" << actualVol.tenor[j+1] << "\t"
						<< "weightRatio=\t" << weightRatio << "\tVswapRatio=\t" << volRatio << "\tImplVolRatio=\t"
						<< impliedVolRatio << "\tFunc=\t" << func[i][j] << "\tFuncTotal=\t" << funcTotal;
			}
		}
	}
	return funcTotal;
}

double BaseModule::strikeRateForSwaptionATM(double maturity, double tenor){
	int matPosition = locate(maturity);
	int tenPosition = locate(tenor);
	double priceAtMat = data.priceD[matPosition];
	double priceAtTen = data.priceD[tenPosition];
	int cashFlows = tenPosition - matPosition;
	double cashFlowsSum = 0.0;
	for(int i=matPosition+1; i<=tenPosition; i++){
		cashFlowsSum += constants.PAYMENT_FREQ*data.priceD[i];
	}
	return (data.priceD[matPosition]-data.priceD[tenPosition])/cashFlowsSum;
}

double BaseModule::blackFormula(double K, double F, double v, int w){

}

void BaseModule::simulatedAnnealingFunc(int calibMethod){
//	double x0 = 0.02;
//	double xmin = 0.00;
//	double xmax = 0.5;
	double A0_0 = 0.01, A1_0 = 0.04, A2_0 = 0.5; //equivalent to x=A0, A1, A2
	double sd0_0= 1, sd1_0=1, sd2_0=0.05;
	double A1minusA0 = A1_0 - A0_0;
	double A1minusA0MinTolerance = 0.75*A1minusA0;
	double A1minusA0MaxTolerance = 1.25*A1minusA0;

		assignVaryingMeanReversion(A0_0, A1_0, A2_0, 0.0);
		double fLast = meanReversionCalibrationFunctionF();

		double gamma = 7.0;
//		double sd0 = 1.5;
		int totalNoOfSimulation = 100;
		int iterationNo = 1;
		vector<double> fxVector, xrVector, sigmajVector;
		fxVector.push_back(fLast);
//		xrVector.push_back(x0);
//		sigmajVector.push_back(sd0);
//		double xr = x0;
//		double x = x0;
//		double xlast = x0;
		double fx = fLast;
		double A0r = A0_0, A1r=A1_0, A2r=A2_0;
		double A0 = A0_0, A1=A1_0, A2=A2_0;
		double A0last = A0_0, A1last=A1_0, A2last=A2_0;
		double A3=0;
		double aVal;
		int counter = 0;
		double minFx = fLast;

		boost::mt19937 *rng2 = new boost::mt19937();
		rng2->seed(time(NULL));
		// http://stackoverflow.com/questions/15747194/why-is-boosts-random-number-generation-on-a-normal-distribution-always-giving

//		aVal = logisticFunc(A0,A1,A2,A3,iterationNo);
		FILE_LOG(logDEBUG) << "SimulatedAnn-Final\t" << iterationNo << "\t" << A0_0 << "\t" << sd0_0 << "\t" << A1_0 << "\t" << sd1_0 << "\t" << A2_0 << "\t" << sd2_0 << "\t" << 0 << "\t" << fLast;
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
			do{
				A0 = generate_next0();
				A1 = generate_next1();
				A2 = generate_next2();
//			}while((A1-A0)>A1minusA0MinTolerance && (A1-A0)<A1minusA0MaxTolerance && (A1>A0));
				if(A0<A1){
					assignVaryingMeanReversion(A0, A1, A2, 0);
					fx = meanReversionCalibrationFunctionF();
					FILE_LOG(logDEBUG) << "SimulatedAnn-Counter\t" << iterationNo << "\t" << A0 << "\t" << sd0j << "\t" << A1 << "\t" << sd1j << "\t" << A2 << "\t" << sd2j << "\t" << counter << "\t" << fx;
					++counter;
				}
			}while(fx>=fLast && counter<(int)round(100/iterationNo));
			fxVector.push_back(fx);
			fLast = fx;
			++iterationNo;
			minFx = min(minFx, fx);
//			aVal = logisticFunc(A0,A1,A2,A3,iterationNo);

			FILE_LOG(logDEBUG) << "SimulatedAnn-Final\t" << iterationNo << "\t" << A0 << "\t" << sd0j << "\t" << A1 << "\t" << sd1j << "\t" << A2 << "\t" << sd2j << "\t" << counter << "\t" << fx;
		}while(iterationNo < totalNoOfSimulation);
		FILE_LOG(logDEBUG) << "SimulatedAnn-Min-Fx\t" << minFx;
}

void BaseModule::assignVaryingMeanReversion(double A0, double A1, double A2, double A3){
	A3 = *(data.time.end());
	for(vector<double>::iterator it = data.time.begin(), it1 = data.aMeanReversion.begin(); it!=data.time.end(); ++it, ++it1){
		double old = *it1;
		*it1 = logisticFunc(A0, A1, A2, A3, *it);
	}
	for(vector<double>::iterator it = data.time.begin(), it1 = data.aMeanReversion.begin(); it!=data.time.end(); ++it, ++it1){
		FILE_LOG(logDEBUG) << "VaryingMeanRev\t" << *it << "\t" << *it1 ;
	}
	calculateEt();
}

void BaseModule::simulatedAnnealingFunc2(int calibMethod){
//	double x0 = 0.02;
//	double xmin = 0.00;
//	double xmax = 0.5;
	double A0_0 = -0.01, A1_0 = 0.1, A2_0 = 0.5; //equivalent to x=A0, A1, A2
	double sd0_0= 1, sd1_0=1, sd2_0=0.05;
	double A1minusA0 = A1_0 - A0_0;
	double A1minusA0MinTolerance = 0.75*A1minusA0;
	double A1minusA0MaxTolerance = 1.25*A1minusA0;

//		assignConstantMeanReversion(x0);
//		double fr = meanReversionCalibrationFunctionF();

		double gamma = 5.0;
//		double sd0 = 1.5;
		int totalNoOfSimulation = 100;
		int iterationNo = 1;
		vector<double> fxVector, xrVector, sigmajVector;
//		fxVector.push_back(fr);
//		xrVector.push_back(x0);
//		sigmajVector.push_back(sd0);
//		double xr = x0;
//		double x = x0;
//		double xlast = x0;
//		double fx;
		double A0r = A0_0, A1r=A1_0, A2r=A2_0;
		double A0 = A0_0, A1=A1_0, A2=A2_0;
		double A0last = A0_0, A1last=A1_0, A2last=A2_0;
		double A3=100;
		double aVal;

		boost::mt19937 *rng2 = new boost::mt19937();
		rng2->seed(time(NULL));
		// http://stackoverflow.com/questions/15747194/why-is-boosts-random-number-generation-on-a-normal-distribution-always-giving

		aVal = logisticFunc(A0,A1,A2,A3,iterationNo);
		FILE_LOG(logINFO) << iterationNo << "\t" << A0_0 << "\t" << sd0_0 << "\t" << A1_0 << "\t" << sd1_0 << "\t" << A2_0 << "\t" << sd2_0 << "\t" << aVal;
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

			do{
				A0 = generate_next0();
				A1 = generate_next1();
				A2 = generate_next2();
			}while((A1-A0)>A1minusA0MinTolerance && (A1-A0)<A1minusA0MaxTolerance && (A1>A0));

			++iterationNo;
			aVal = logisticFunc(A0,A1,A2,A3,iterationNo);

			FILE_LOG(logINFO) << iterationNo << "\t" << A0 << "\t" << sd0j << "\t" << A1 << "\t" << sd1j << "\t" << A2 << "\t" << sd2j << "\t" << aVal;
		}while(iterationNo < totalNoOfSimulation);
}

double BaseModule::logisticFunc(double A0, double A1, double A2, double A3, double ti){
	return A0+(A1-A0)/(1+exp(A2*(A3-ti)));
}

void BaseModule::simulatedAnnealingFunc1(int calibMethod){
	double x0 = 0.02;
	double xmin = 0.00;
	double xmax = 0.5;

		assignConstantMeanReversion(x0);
		double fr = meanReversionCalibrationFunctionF();

		double gamma = 5.0;
		double sigma0 = 1.5;
		int totalNoOfSimulation = 100;
		int iterationNo = 1;
		vector<double> fxVector, xrVector, sigmajVector;
		fxVector.push_back(fr);
		xrVector.push_back(x0);
		sigmajVector.push_back(sigma0);
		double xr = x0;
		double x = x0;
		double xlast = x0;
		double fx;
		boost::mt19937 *rng2 = new boost::mt19937();
		rng2->seed(time(NULL));
		// http://stackoverflow.com/questions/15747194/why-is-boosts-random-number-generation-on-a-normal-distribution-always-giving

		FILE_LOG(logINFO) << x0 << "\t" << sigma0 ;
		do {
			double sigmaj = coolingMechanism(gamma, sigma0, totalNoOfSimulation, iterationNo);
			normal_distribution<> distribution(x, sigmaj);
			variate_generator<boost::mt19937&, normal_distribution<> > generate_next(*rng2, distribution);

			if (x == x0 || xlast == x){
				x = generate_next();
			}
			while (!(x > xmin && x < xmax)){
				x = generate_next();
			}
			xlast = x;

			FILE_LOG(logINFO) << x << "\t" << sigmaj << "\t";
			++iterationNo;
		}while(iterationNo < totalNoOfSimulation);
}

double BaseModule::coolingMechanism(double gamma, double sigma0, int totalNoOfSimulation, int iterationNo){
	return sigma0 * exp(-gamma*iterationNo/(totalNoOfSimulation-1));
}

std::mt19937 BaseModule::get_prng() {
//    std::uint_least32_t seed_data[std::mt19937::state_size];
//    std::random_device r;
//    std::generate_n(seed_data, std::mt19937::state_size, std::ref(r));
//    std::seed_seq q(std::begin(seed_data), std::end(seed_data));
//    return std::mt19937{q};
}

template<class T>
T gen_normal_4(T generator,
            std::vector<double> &res)
{
  for(size_t i=0; i<res.size(); ++i)
    res[i]=generator();
  // Note the generator is returned back
  return  generator;
}

//void main1(void)
//{
////  boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
////    generator(boost::mt19937(time(0)),
////              boost::normal_distribution<>());
////
////  std::vector<double> res(10);
////  // Assigning back to the generator ensures the state is advanced
////  generator=gen_normal_4(generator, res);
////
////  for(size_t i=0; i<10; ++i)
////    std::cout<<res[i]
////             <<std::endl;
////
////  generator=gen_normal_4(generator, res);
////
////  for(size_t i=0; i<10; ++i)
////    std::cout<<res[i]
////             <<std::endl;
//
//    // Construct a standard normal distribution s
//      normal s; // (default mean = zero, and standard deviation = unity)
//      cout << "Standard normal distribution, mean = "<< s.mean()
//        << ", standard deviation = " << s.standard_deviation() << endl;
//
//}

void main1(void)
{
	boost::mt19937 *rng2 = new boost::mt19937();
	rng2->seed(time(NULL));

	normal_distribution<> distribution(0, 1);
	variate_generator<boost::mt19937&, normal_distribution<> > resampler(*rng2, distribution);
//	boost::variate_generator< boost::mt19937, boost::normal_distribution<> > resampler(*rng2, distribution);
//
	for (int i = 0; i < 10; ++i)
	{
		cout << resampler() << endl ;
	}
}
