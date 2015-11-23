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

Module5Appl::~Module5Appl() {
}

int Module5Appl::moduleMainFunc(){
	map<char, vector<double> > mapVal = helper.initializeDataStructure();

	data.time = mapVal.find(helper.TIME)->second;
	data.priceD = mapVal.find(helper.PRICE)->second;
	data.forward = mapVal.find(helper.FORWARD)->second;

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
	simulatedAnnealingFuncForVolatility();

	return 0;
}
