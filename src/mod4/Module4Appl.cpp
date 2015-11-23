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

	for(int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	assignConstantMeanReversion(Constants::FIXED_MEAN_REVERSION);
	assignConstantVol(Constants::FIXED_VOLATILITY);
	calculateEt();

	double t = 0.0;
	double T = 0.0;
	cout << "a\t" << *(data.aMeanReversion.begin())<<endl;
	cout << "sigma\t" << *(data.sigma.begin())<<endl;

	initializeAndAssignConstantWeights();
	vector<vector<double> > marketVolatilityColwise = transposeVector(actualVol.vol);
	BlacksFormula black;
	vector<vector<double> > blckPrices;
	vector<double>::iterator itTenor = actualVol.tenor.begin();

	for (vector<vector<double> >::iterator itVol = marketVolatilityColwise.begin(); itVol!=marketVolatilityColwise.end(); ++itVol){
		vector<double> discountPrices;
		for (vector<double>::iterator itMaturity = actualVol.maturity.begin(); itMaturity!=actualVol.maturity.end(); ++itMaturity){
			int pos = locate(*itTenor + *itMaturity);
			discountPrices.push_back(data.priceD[pos]);
		}
		vector<double> pricesCol = black.priceBlackCap(*itVol, discountPrices, actualVol.maturity, constants.FIXED_STRIKE, 1.0, *itTenor, false);
		++itTenor;
		blckPrices.push_back(pricesCol);
	}
	blckPrices = transposeVector(blckPrices);
	cout << "Maturity\\Tenor\t" ;
	for(vector<double>::iterator it_ten = actualVol.tenor.begin(); it_ten != actualVol.tenor.end(); ++it_ten){
		cout << *it_ten << "\t";
	}
	cout << endl;
//	for(vector<double>::iterator it_mat = actualVol.maturity.begin(); it_mat!=actualVol.maturity.end(); ++it_mat){
//	}
	for(vector<vector<double> >::iterator blckit = blckPrices.begin(); blckit != blckPrices.end(); ++blckit){
		cout << "\t";
		for(vector<double>::iterator blckit2 = (*blckit).begin(); blckit2!=(*blckit).end(); ++blckit2){
			cout << *blckit2 << "\t";
		}
		cout << endl;
	}

	return 0;
}
