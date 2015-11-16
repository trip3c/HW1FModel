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
#include "../core/spline.h"

#include "Module1Appl.h"

using namespace std;

Module1Appl::Module1Appl() {
}

Module1Appl::~Module1Appl() {
//	cout << "Module1Appl destructor" << endl;
}

int Module1Appl::moduleMainFunc(){
	map<char, vector<double> > mapVal = helper.initializeDataStructure();

	mapVal.find(helper.TIME);
	data.time = mapVal.find(helper.TIME)->second;
	data.priceD = mapVal.find(helper.PRICE)->second;
	data.forward = mapVal.find(helper.FORWARD)->second;

	implVol.vol = helper.initializeVolatility();
	implVol.maturity = helper.initializeSwaptionVolatilityMaturity();
	implVol.tenor = helper.initializeSwaptionVolatilityTenor();

	for(int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	assignMeanReversion();
	assignVol(0);
	double t = 0.0;
	double T = 0.0;
	double K = 0.04;
	cout << "a\t" << *(data.aMeanReversion.begin())<<endl;
	cout << "sigma\t" << *(data.sigma.begin())<<endl;
	cout << "Maturity\\Tenor\t" ;
	for(vector<double>::iterator it_ten = implVol.tenor.begin(); it_ten != implVol.tenor.end(); ++it_ten){
		cout << *it_ten << "\t";
	}
	cout << endl;
	for(vector<double>::iterator it_mat = implVol.maturity.begin(); it_mat!=implVol.maturity.end(); ++it_mat){
		cout << *it_mat << "\t";
		for(vector<double>::iterator it_ten = implVol.tenor.begin(); it_ten != implVol.tenor.end(); ++it_ten){
			calculateEt();

			t = *it_mat;
			T = *it_ten + t;
			double price = pSwaption(K, t, T);
			cout << price << "\t";
		}
		cout << endl;
	}
//
//	assignMeanReversion();
//	assignVol(0);
//	calculateEt();
////	calculateVr(1);
//	vector<double> th = theta();
//
//	cout << "a\t" << *(data.a.begin())<<"\t";
//	cout << "sigma\t" << *(data.sigma.begin())<<"\t";
//	for(vector<double>::iterator it = th.begin(); it != th.end(); ++it){
//		cout << *it << "\t" ;
//	}
//	cout << endl;

	return 0;
}

void Module1Appl::assignMeanReversion(){
	double a = 5.0/100;
	int x = data.time.size();
	data.aMeanReversion.assign(x, a);
}

void Module1Appl::assignVol(double s){
	double s1 = 70.0/100;
	int x = data.time.size();
	data.sigma.assign(x, s1);
}


