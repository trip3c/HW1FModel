/*
 * Module2Appl.cpp
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

#include "Module2Appl.h"

using namespace std;

Module2Appl::Module2Appl() {
}

Module2Appl::~Module2Appl() {
}

int Module2Appl::moduleMainFunc(){
	map<char, vector<double> > mapVal = helper.initializeDataStructure();

	mapVal.find(helper.TIME);
	data.time = mapVal.find(helper.TIME)->second;
	data.priceD = mapVal.find(helper.PRICE)->second;
	data.forward = mapVal.find(helper.FORWARD)->second;

	for(int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	assignMeanReversion();
	assignVol(0);
	calculateEt();
//	calculateVr(1);
	vector<double> th = theta();

	cout << "a\t" << *(data.aMeanReversion.begin())<<"\t";
	cout << "sigma\t" << *(data.sigma.begin())<<"\t";
	for(vector<double>::iterator it = th.begin(); it != th.end(); ++it){
		cout << *it << "\t" ;
	}
	cout << endl;

	return 0;
}

void Module2Appl::assignMeanReversion(){
	double a = 5.0/100;
	int x = data.time.size();
	data.aMeanReversion.assign(x, a);
}

void Module2Appl::assignVol(double s){
	double s1 = 70.0/100;
	int x = data.time.size();
	data.sigma.assign(x, s1);
}


