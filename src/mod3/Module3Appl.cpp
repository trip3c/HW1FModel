/*
 * Module3Appl.cpp
 *
 *  Created on: Nov 15, 2015
 *      Author: Chandra
 */
#include <iostream>
#include <string>
#include <set>
#include <cassert>
#include <math.h>
#include <iomanip>
#include "../core/spline.h"

#include "Module3Appl.h"

using namespace std;

Module3Appl::Module3Appl() {
}

Module3Appl::~Module3Appl() {
//	cout << "Module3Appl destructor" << endl;
}

int Module3Appl::moduleMainFunc(){
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

	assignConstantMeanReversion();
	assignConstantVol();
	calculateEt();

	double t = 0.0;
	double T = 0.0;
	cout << "a\t" << *(data.aMeanReversion.begin())<<endl;
	cout << "sigma\t" << *(data.sigma.begin())<<endl;
	cout << "Maturity\\Tenor\t" ;
	for(vector<double>::iterator it_ten = actualVol.tenor.begin(); it_ten != actualVol.tenor.end(); ++it_ten){
		cout << *it_ten << "\t";
	}
	cout << endl;
	for(vector<double>::iterator it_mat = actualVol.maturity.begin(); it_mat!=actualVol.maturity.end(); ++it_mat){
		cout << *it_mat << "\t";
		for(vector<double>::iterator it_ten = actualVol.tenor.begin(); it_ten != actualVol.tenor.end(); ++it_ten){

			t = *it_mat;
			T = *it_ten + t;
			double implVol = sqrt(calculateVswap(t, T));
			cout << implVol << "\t";
		}
		cout << endl;
	}
	return 0;
}


