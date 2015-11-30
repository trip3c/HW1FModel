/*
 * ModuleHelper.cpp
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
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include "core/spline.h"
#include "ModuleHelper.h"
#include "BootstrapLoader.h"

#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif
using namespace std;


ModuleHelper::ModuleHelper() {
}

ModuleHelper::ModuleHelper(BootstrapLoader pLoader) {
	serviceLocator = pLoader;
}
ModuleHelper::~ModuleHelper() {
}


map<char, vector<double> > ModuleHelper::interpolatePricesAndForwardRates(vector<vector<double> > newNumbers){
	vector<vector<double> >::iterator vvi_iterator;
	vector<double>::iterator vi_iterator;
	vvi_iterator = newNumbers.begin();
	vi_iterator = (*vvi_iterator).begin();
	vector<double> X((*vvi_iterator).size());
	X = *vvi_iterator;

	// Initializing forwards data - 2nd column
	++vvi_iterator;
	vi_iterator = (*vvi_iterator).begin();
	vector<double> Y((*vvi_iterator).size());
	Y = *vvi_iterator;

	// Inserting year into set
	vi_iterator = X.begin();
	double left = *vi_iterator;
	double right = 0.0;
	++vi_iterator;
	std::set<double> yearSet;
	for(; vi_iterator!=X.end(); ++vi_iterator){
		right = *vi_iterator;
		for(double i=left; i<=right; i=i+0.25){
			yearSet.insert(i);
		}
		left = *vi_iterator;
	}

	// Initializing year and forwards data
	tk::spline s;
	s.set_points(X,Y);    // currently it is required that X is already sorted

	vector<double> time;
	vector<double> forward;
	vector<double> priceD;
	map<char, vector<double> > returnMap;

	// Interpolating forwards data
	for (std::set<double>::iterator it=yearSet.begin(); it!=yearSet.end(); ++it){
		time.push_back(*it);
		forward.push_back(s(*it));
	}

	// Initializing discounts data - 3rd column
	++vvi_iterator;
	vi_iterator = (*vvi_iterator).begin();
	vector<double> Z((*vvi_iterator).size());
	Z = *vvi_iterator;

	// Interpolating discounts data
	s.set_points(X,Z);    // currently it is required that X is already sorted
	for (std::set<double>::iterator it=yearSet.begin(); it!=yearSet.end(); ++it){
		double t=s(*it);
		priceD.push_back(t);
	}
	returnMap[TIME] = time;
	returnMap[FORWARD] = forward;
	returnMap[PRICE] = priceD;
	return returnMap;
}

map<char, vector<double> > ModuleHelper::initializeDataStructure(){
//	string directoryPath = getDirectoryPath();
//	directoryPath += "\\src\\defs\\UCLA_Discounts.txt";
	string directoryPath = serviceLocator.getDiscounts();
	vector<vector<double> > importData = importDataFromFile(directoryPath);
	vector<vector<double> > forwardRateAndDiscountFactors =
			rowsToColumnTransposeVector(importData);
	map<char, vector<double> > returnMap =
			interpolatePricesAndForwardRates(forwardRateAndDiscountFactors);
	return returnMap;
}

vector<vector<double> > ModuleHelper::initializeVolatility(){
//	string directoryPath = getDirectoryPath();
//	directoryPath += "\\src\\defs\\SwaptionVolatilityOnly.txt";
	string directoryPath = serviceLocator.getSwaptionVolatility();
	vector<vector<double> > volData = importDataFromFile(directoryPath);
	return volData;
}

vector<double> ModuleHelper::initializeSwaptionVolatilityMaturity(){
//	string directoryPath = getDirectoryPath();
//	directoryPath += "\\src\\defs\\SwaptionVolatilityMaturity.txt";
	string directoryPath = serviceLocator.getSwaptionMaturity();
	vector<vector<double> > volMaturityData = importDataFromFile(directoryPath);
	return *(volMaturityData.begin());
}

vector<double> ModuleHelper::initializeSwaptionVolatilityTenor(){
//	string directoryPath = getDirectoryPath();
//	directoryPath += "\\src\\defs\\SwaptionVolatilityTenor.txt";
	string directoryPath = serviceLocator.getSwaptionTenor();
	vector<vector<double> > volTenorData = importDataFromFile(directoryPath);
	return *(volTenorData.begin());
}


vector<vector<double> > ModuleHelper::importDataFromFile(string filePath){
	vector<vector<double> > numbers;
	std::ifstream ifs(filePath.c_str());
	string temp;
	while (std::getline(ifs, temp)) {
		std::istringstream buffer(temp);
		vector<double> line((std::istream_iterator<double>(buffer)), std::istream_iterator<double>());
		numbers.push_back(line);
	}
	ifs.close();
//	printVectorOfVectors(numbers);
//	cout << "---------------------" << endl;
	return numbers;
}

vector<vector<double> > ModuleHelper::rowsToColumnTransposeVector(vector<vector<double> > numbers){
	vector<vector<double> >::iterator vvi_iterator;
	vector<double>::iterator vi_iterator;
	assert(numbers.size()>0);

	int columns = numbers.front().size();
	vector<vector<double> > newNumbers;
	for(int i = 0; i<columns; i++){
		vector<double> newNumbersCols;
		int j;
		for(vvi_iterator = numbers.begin(); vvi_iterator!=numbers.end(); ++vvi_iterator) {
			for(vi_iterator = (*vvi_iterator).begin(), j=0; j<i; ++vi_iterator,++j){}
			newNumbersCols.push_back(*vi_iterator);
		}
		newNumbers.push_back(newNumbersCols);
	}
	return newNumbers;
}

string ModuleHelper::getDirectoryPath(){
	char cCurrentPath[FILENAME_MAX];
	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
	{
		return "";
	}
	cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */
	return cCurrentPath;
}

void ModuleHelper::printVectorOfVectors(vector<vector<double> > vec){
	vector<vector<double> >::iterator vvi_iterator;
	vector<double>::iterator vi_iterator;

	for(vvi_iterator = vec.begin(); vvi_iterator!=vec.end(); ++vvi_iterator) {
		vector<double> X((*vvi_iterator).size());
		X = *vvi_iterator;
		for(vi_iterator = (*vvi_iterator).begin(); vi_iterator != X.end(); ++vi_iterator){
			cout << *vi_iterator << " " << endl;
		}
		cout << endl;
	}
}
