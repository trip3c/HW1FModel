/*
 * Module5Appl.cpp
 *
 *  Created on: Nov 22, 2015
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
#include "Module5Appl.h"
#include "core/BlacksFormula.h"

using namespace std;

Module5Appl::Module5Appl() {
}
Module5Appl::Module5Appl(BootstrapLoader sl) {
	serviceLocator = sl;
	helper = ModuleHelper(sl);
}
Module5Appl::~Module5Appl() {
}

int Module5Appl::moduleMainFunc(){
	map<char, vector<double> > mapVal = helper.initializeDataStructure();

	data.time = mapVal.find(helper.TIME)->second;
	data.priceD = mapVal.find(helper.PRICE)->second;
	data.forward = mapVal.find(helper.FORWARD)->second;

//	double A0_0=0.00977275, A1_0=-0.000614657, A2_0=-0.0000132694,A3_0=0.00000547935;
//	double A0_0=0.00974798, A1_0=-0.000850625, A2_0=-0.000000341365, A3_0=0.0000153473;
	double A0_0=0.02074798,A1_0=-0.00850625,A2_0=-0.000000341365,A3_0=0.0000953473;

//	double A0_0=0.008, A1_0=0.008, A2_0=0.008,A3_0=0.008;
	volatilityParams.A0=A0_0;
	volatilityParams.A1=A1_0;
	volatilityParams.A2=A2_0;
	volatilityParams.A3=A3_0;
	volatilityParams.Gx=0.0;

	assignConstantMeanReversion(0);
	assignConstantVol(0);
	actualVol.vol = helper.initializeVolatility();
	actualVol.maturity = helper.initializeSwaptionVolatilityMaturity();
	actualVol.tenor = helper.initializeSwaptionVolatilityTenor();
	for(int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	initializeStrikeRateForSwaptionATM();

	for(int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	initializeAndAssignConstantWeights();
	for(int i = 0; i<1; ++i){
		simulatedAnnealingFuncForVolatility();
	}
//	cout << RationalApproximation(0.95);
	return 0;
}
