/*
 * Module4Appl.cpp
 *
 *  Created on: Nov 21, 2015
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
#include "Module4Appl.h"
#include "core/BlacksFormula.h"

using namespace std;

Module4Appl::Module4Appl() {
}
Module4Appl::Module4Appl(BootstrapLoader sl) {
	serviceLocator = sl;
	helper = ModuleHelper(sl);
}
Module4Appl::~Module4Appl() {
}

int Module4Appl::moduleMainFunc(){
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
	calculateEt();

	initializeAndAssignConstantWeights();
	initializeStrikeRateForSwaptionATM();
	cout << "Maturity\\Tenor\t" ;
	for(vector<double>::iterator it_ten = actualVol.tenor.begin(); it_ten != actualVol.tenor.end(); ++it_ten){
		cout << *it_ten << "\t";
	}
	cout << endl;
	for(unsigned int i=0; i<actualVol.maturity.size(); ++i){
		cout << actualVol.maturity[i] << "\t";
		for(unsigned int j=0; j<actualVol.tenor.size(); ++j){
			double price = blackSwaptionPriceATM(actualVol.maturity[j], actualVol.tenor[j], actualVol.vol[i][j], actualVol.strikeRate[i][j], true);
			cout << price << "\t";
		}
		cout << endl;
	}

	return 0;
}
