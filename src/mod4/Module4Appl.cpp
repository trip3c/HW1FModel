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

	cout << "a\t" << *(data.aMeanReversion.begin())<<endl;
	cout << "sigma\t" << *(data.sigma.begin())<<endl;

	initializeAndAssignConstantWeights();
	initializeStrikeRateForSwaptionATM();
//	vector<vector<double> > marketVolatilityColwise = transposeVector(actualVol.vol);
//	vector<vector<double> > strikeRatesTranspose = transposeVector(actualVol.strikeRate);
//	BlacksFormula black;
//	vector<vector<double> > blckPrices;
//	vector<double>::iterator itTenor = actualVol.tenor.begin();
//
//	int i=0;
//	for (vector<vector<double> >::iterator itVol = marketVolatilityColwise.begin(); itVol!=marketVolatilityColwise.end(); ++itVol){
//		vector<double> discountPrices;
//		int j=0;
//		for (vector<double>::iterator itMaturity = actualVol.maturity.begin(); itMaturity!=actualVol.maturity.end(); ++itMaturity){
//			int pos = locate(*itTenor + *itMaturity);
//			discountPrices.push_back(data.priceD[pos]);
//			++j;
//		}
//		vector<double> pricesCol = black.priceBlackCap(*itVol, discountPrices, actualVol.maturity,
//				strikeRatesTranspose[i][j], 1.0, *itTenor, false);
//		++itTenor;
//		blckPrices.push_back(pricesCol);
//		++i;
//	}
//	blckPrices = transposeVector(blckPrices);
//	cout << "Maturity\\Tenor\t" ;
//	for(vector<double>::iterator it_ten = actualVol.tenor.begin(); it_ten != actualVol.tenor.end(); ++it_ten){
//		cout << *it_ten << "\t";
//	}
//	cout << endl;
////	for(vector<double>::iterator it_mat = actualVol.maturity.begin(); it_mat!=actualVol.maturity.end(); ++it_mat){
////	}
//	for(vector<vector<double> >::iterator blckit = blckPrices.begin(); blckit != blckPrices.end(); ++blckit){
//		cout << "\t";
//		for(vector<double>::iterator blckit2 = (*blckit).begin(); blckit2!=(*blckit).end(); ++blckit2){
//			cout << *blckit2 << "\t";
//		}
//		cout << endl;
//	}

	vector<vector<double> > blckPrices;
	int i=0;
	for (vector<double>::iterator itMaturity = actualVol.maturity.begin(); itMaturity!=actualVol.maturity.end(); ++itMaturity){
		int j=0;
		vector<double> blckPricesInner;
		for (vector<double>::iterator itTenor = actualVol.tenor.begin(); itTenor!=actualVol.tenor.end(); ++itTenor){
			double prc = blackSwaptionPriceATM(*itMaturity, *itTenor, actualVol.vol[i][j], actualVol.strikeRate[i][j], true);
			blckPricesInner.push_back(prc);
			++j;
		}
		blckPrices.push_back(blckPricesInner);
		++i;
	}

	for(vector<vector<double> >::iterator iterator = blckPrices.begin(); iterator!=blckPrices.end(); ++iterator){
		for(vector<double>::iterator iterator2 = (*iterator).begin(); iterator2!=(*iterator).end(); ++iterator2){
			cout << *iterator2 << "\t";
		}
		cout << endl;
	}
	cout << endl;

	return 0;
}
