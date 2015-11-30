/*
 * Module6Appl.cpp
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
#include "Module6Appl.h"
#include "core/BlacksFormula.h"

using namespace std;

Module6Appl::Module6Appl() {
}
Module6Appl::Module6Appl(BootstrapLoader sl) {
	serviceLocator = sl;
	helper = ModuleHelper(sl);
}
Module6Appl::~Module6Appl() {
}

int Module6Appl::moduleMainFunc(){
	map<char, vector<double> > mapVal = helper.initializeDataStructure();

	data.time = mapVal.find(helper.TIME)->second;
	data.priceD = mapVal.find(helper.PRICE)->second;
	data.forward = mapVal.find(helper.FORWARD)->second;

	actualVol.vol = helper.initializeVolatility();
	actualVol.maturity = helper.initializeSwaptionVolatilityMaturity();
	actualVol.tenor = helper.initializeSwaptionVolatilityTenor();
	for(unsigned int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	assignConstantMeanReversion(serviceLocator.getMeanReversionConstant());
	assignConstantVol(serviceLocator.getVolatilityConstant());
	assignVaryingMeanReversion(serviceLocator.getMeanReversionA0(), serviceLocator.getMeanReversionA1(), serviceLocator.getMeanReversionA2());
	assignVaryingVolatility(serviceLocator.getVolatilityA0(),serviceLocator.getVolatilityA1(),
			serviceLocator.getVolatilityA2(),serviceLocator.getVolatilityA3());
	calculateEt();


	initializeStrikeRateForSwaptionATM();

	for(unsigned int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	initializeAndAssignConstantWeights();

	double t = 0.0;
	double T = 0.0;
	cout << "Maturity\\Tenor\t" ;
	for(vector<double>::iterator it_ten = actualVol.tenor.begin(); it_ten != actualVol.tenor.end(); ++it_ten){
		cout << *it_ten << "\t";
	}
	cout << endl;
	for(unsigned int i=0; i<actualVol.maturity.size(); ++i){
		cout << actualVol.maturity[i] << "\t";
		for(unsigned int j=0; j<actualVol.tenor.size(); ++j){
			t = actualVol.tenor[j];
			T = actualVol.tenor[j] + actualVol.maturity[i];
			double implVol = sqrt(calculateVswap(t, T));
			cout << implVol << "\t";
		}
		cout << endl;
	}

	return 0;
}
