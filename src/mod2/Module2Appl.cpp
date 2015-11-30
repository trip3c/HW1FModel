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
#include "core/spline.h"
#include "core/Constants.h"

#include "Module2Appl.h"

using namespace std;

Module2Appl::Module2Appl() {
}
Module2Appl::Module2Appl(BootstrapLoader sl) {
	serviceLocator = sl;
	helper = ModuleHelper(sl);
}
Module2Appl::~Module2Appl() {
}

int Module2Appl::moduleMainFunc(){
	map<char, vector<double> > mapVal = helper.initializeDataStructure();

	mapVal.find(helper.TIME);
	data.time = mapVal.find(helper.TIME)->second;
	data.priceD = mapVal.find(helper.PRICE)->second;
	data.forward = mapVal.find(helper.FORWARD)->second;

	for(unsigned int i=0; i<data.time.size(); i++){
		timePos[data.time[i]] = i;
	}

	assignConstantMeanReversion(serviceLocator.getMeanReversionConstant());
	assignConstantVol(serviceLocator.getVolatilityConstant());
	assignVaryingMeanReversion(serviceLocator.getMeanReversionA0(), serviceLocator.getMeanReversionA1(),
			serviceLocator.getMeanReversionA2(), *(--(data.time.end())));
	assignVaryingVolatility(serviceLocator.getVolatilityA0(), serviceLocator.getVolatilityA1(),
			serviceLocator.getVolatilityA2(), serviceLocator.getVolatilityA3());
	calculateEt();
//	calculateVr(1);
	vector<double> th = theta();

//	cout << "a\t" << *(data.aMeanReversion.begin())<<"\t";
//	cout << "sigma\t" << *(data.sigma.begin())<<"\t";
	for(vector<double>::iterator it = th.begin(); it != th.end(); ++it){
		cout << *it << "\t" ;
	}
	cout << endl;

	return 0;
}
