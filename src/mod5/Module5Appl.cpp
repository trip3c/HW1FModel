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

	volatilityParams.A0=serviceLocator.getVolatilityA0();
	volatilityParams.A1=serviceLocator.getVolatilityA1();
	volatilityParams.A2=serviceLocator.getVolatilityA2();
	volatilityParams.A3=serviceLocator.getVolatilityA3();
	volatilityParams.Gx=0.0;

	assignConstantMeanReversion(0);
	assignConstantVol(0);
	actualVol.vol = helper.initializeVolatility();
	actualVol.maturity = helper.initializeSwaptionVolatilityMaturity();
	actualVol.tenor = helper.initializeSwaptionVolatilityTenor();
	for(unsigned int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	initializeStrikeRateForSwaptionATM();

	for(unsigned int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	initializeAndAssignConstantWeights();
	for(unsigned int i = 0; i<1; ++i){
		simulatedAnnealingFuncForVolatility();
	}
	return 0;
}
