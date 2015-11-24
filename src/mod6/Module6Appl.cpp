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
	for(int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	assignConstantMeanReversion(0.05);
	assignConstantVol(0.008);
	assignVaryingMeanReversion(0.0100008, 0.0788238, 0.52434, 0);
//	assignVaryingVolatility(0.00977275,-0.000614657,-0.0000132694,0.00000547935);
//	assignVaryingVolatility(0.00976917,-0.000806147,-0.0000035288,0.0000138456);
//	assignVaryingVolatility(0.00951463,-0.000620164,-0.000000911601,0.00000742415);

	assignVaryingVolatilityUpwardSloping(0.0044705,0.000525428,0.0000173218,-0.00000192464);
	calculateEt();


	initializeStrikeRateForSwaptionATM();

	for(int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	initializeAndAssignConstantWeights();
//	simulatedAnnealingFuncForVolatility();

	double t = 0.0;
	double T = 0.0;
	cout << "Maturity\\Tenor\t" ;
	for(vector<double>::iterator it_ten = actualVol.tenor.begin(); it_ten != actualVol.tenor.end(); ++it_ten){
		cout << *it_ten << "\t";
	}
	cout << endl;
	for(int i=0; i<actualVol.maturity.size(); ++i){
		cout << actualVol.maturity[i] << "\t";
		for(int j=0; j<actualVol.tenor.size(); ++j){
			t = actualVol.tenor[j];
			T = actualVol.tenor[j] + actualVol.maturity[i];
			double implVol = sqrt(calculateVswap(t, T));
			cout << implVol << "\t";
		}
		cout << endl;
	}

	return 0;
}
