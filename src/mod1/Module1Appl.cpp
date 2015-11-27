/*
 * Module1Appl.cpp
 *
 *  Created on: Nov 10, 2015
 *      Author: Chandra
 */
#include <iostream>
#include <string>
#include <set>
#include <cassert>
#include <math.h>
#include <iomanip>
#include "core/spline.h"
#include "core/Constants.h"

#include "Module1Appl.h"

using namespace std;

Module1Appl::Module1Appl() {
}

Module1Appl::Module1Appl(BootstrapLoader sl) {
	serviceLocator = sl;
	helper = ModuleHelper(sl);
}

Module1Appl::~Module1Appl() {
}

int Module1Appl::moduleMainFunc(){
	map<char, vector<double> > mapVal = helper.initializeDataStructure();

	data.time = mapVal.find(helper.TIME)->second;
	data.priceD = mapVal.find(helper.PRICE)->second;
	data.forward = mapVal.find(helper.FORWARD)->second;

	actualVol.vol = helper.initializeVolatility();
	actualVol.maturity = helper.initializeSwaptionVolatilityMaturity();
	actualVol.tenor = helper.initializeSwaptionVolatilityTenor();

	for(int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	assignConstantMeanReversion(serviceLocator.getMeanReversionConstant());
	assignConstantVol(serviceLocator.getVolatilityConstant());
//	assignVaryingMeanReversion(0.0100008, 0.0788238, 0.52434, 0);/////
//	assignVaryingMeanReversion(0.0100008, 0.1088238, 0.5, 0);
//	assignVaryingMeanReversion(0.05500008, 0.1088238, 0.5, 0);//ok
//	assignVaryingMeanReversion(0.0100008, 0.1088238, 0.5, 0);
//	assignVaryingVolatility(0.0552728,-0.00404704,-0.000259549,0.0000485404);
//	assignVaryingVolatility(0.00977275,-0.000614657,-0.0000132694,0.00000547935);
//	assignVaryingVolatility(0.00976917,-0.000806147,-0.0000035288,0.0000138456);
//	assignVaryingVolatility(0.00974798,-0.000850625,-0.000000341365,0.0000153473);/////
//	assignVaryingVolatility(0.02074798,-0.00850625,-0.000000341365,0.0000753473);
//	assignVaryingVolatility(0.02074798,-0.00850625,-0.000000341365,0.0000953473);//Ok
//	assignVaryingVolatility(0.02074798,-0.00850625,-0.000000341365,0.0000953473);
//	assignVaryingVolatility(0.0120138,-0.001360769,-0.0000448605,0.0000049845);
//	assignVaryingVolatility(0.02074798,-0.000850625,-0.0000341365,0.0000953473);
//	assignVaryingVolatilitySpline(0.02074798,0.04074798,0.08074798,0.12074798);

//	assignVaryingMeanReversion(0.05500008, 0.1088238, 0.5, 0);//ok
//	assignVaryingVolatility(0.02074798,-0.00850625,-0.000000341365,0.0000953473);//Ok
	assignVaryingMeanReversion(serviceLocator.getMeanReversionA0(),
			serviceLocator.getMeanReversionA1(),
			serviceLocator.getMeanReversionA2());
	assignVaryingVolatility(serviceLocator.getVolatilityA0(),
			serviceLocator.getVolatilityA1(),
			serviceLocator.getVolatilityA2(),
			serviceLocator.getVolatilityA3());
	calculateEt();

	double t = 0.0;
	double T = 0.0;
	initializeStrikeRateForSwaptionATM();
	cout << "a\t" << *(data.aMeanReversion.begin())<<endl;
	cout << "sigma\t" << *(data.sigma.begin())<<endl;
	cout << "Maturity\\Tenor\t" ;
	for(vector<double>::iterator it_ten = actualVol.tenor.begin(); it_ten != actualVol.tenor.end(); ++it_ten){
		cout << *it_ten << "\t";
	}
	cout << endl;
	for(int i=0; i<actualVol.maturity.size(); ++i){
		cout << actualVol.maturity[i] << "\t";
		for(int j=0; j<actualVol.tenor.size(); ++j){
			t = actualVol.maturity[j];
			T = actualVol.tenor[j] + actualVol.maturity[i];
			double strike = actualVol.strikeRate[i][j];
			double price = pSwaption(strike, t, T);
//			double implVol = blackImpliedVolatilitySwaptionATM(t, T, price, strike);
//			cout << implVol << "\t";
			cout << price << "\t";
		}
		cout << endl;
	}

//	assignConstantMeanReversion(0.025); //dummy
//	assignVaryingMeanReversion(0.0100008, 0.0788238, 0.52434, 0);
//	//double A0_0=0.814719, A1_0=-0.0584122, A2_0=-0.00192567 ,A3_0=0.000213964;
//	double A0_0=0.14203958, A1_0=-0.0115965, A2_0=-0.0003823 ,A3_0=0.000042478;
//	assignConstantVol(0.025); //dummy
//	assignVaryingVolatility(A0_0, A1_0, A2_0, A3_0);
//	calculateEt();
//
//	initializeStrikeRateForSwaptionATM();
//	cout << pSwaption(actualVol.strikeRate[2][9], 0.75, 10.75);

	return 0;
}
