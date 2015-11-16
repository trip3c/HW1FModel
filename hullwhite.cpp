//============================================================================
// Name        : hullwhite.cpp
// Author      : chandra, james, alan, gloria
// Version     :
// Copyright   : Hull White One Factor Calibration
// Description : Hull White One Factor Calibration
//============================================================================

#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "spline.h"
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>

using namespace std;

vector<vector<double> > rowsToColumnVector(vector<vector<double> >);
//pair<vector<double>, vector<double> > termStructurePrice(vector<vector<double> > newNumbers);
void termStructurePrice(vector<vector<double> > newNumbers);

int main() {
	vector<vector<double> > numbers;
	string temp;

	std::ifstream ifs("E:\\work\\cpp_ws\\hullwhite\\src\\UCLA_Discounts.txt");

	while (std::getline(ifs, temp)) {
		std::istringstream buffer(temp);
		vector<double> line((std::istream_iterator<double>(buffer)), std::istream_iterator<double>());
		numbers.push_back(line);
	}
	ifs.close();
	vector<vector<double> > newNumbers = rowsToColumnVector(numbers);
//	pair<vector<double>, vector<double> > prices = termStructurePrice(newNumbers);
	termStructurePrice(newNumbers);
	return 0;
}

vector<vector<double> > rowsToColumnVector(vector<vector<double> > numbers){
	vector<vector<double> >::iterator vvi_iterator;
	vector<double>::iterator vi_iterator;

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

//pair<vector<double>, vector<double> > termStructurePrice(vector<vector<double> > newNumbers){
void termStructurePrice(vector<vector<double> > newNumbers){
	vector<vector<double> >::iterator vvi_iterator;
	vector<double>::iterator vi_iterator;
	vvi_iterator = newNumbers.begin();
	vi_iterator = (*vvi_iterator).begin();
	vector<double> X((*vvi_iterator).size());
	X = *vvi_iterator;

	++vvi_iterator;++vvi_iterator;
//	++vvi_iterator;
	vi_iterator = (*vvi_iterator).begin();
	vector<double> Y((*vvi_iterator).size());
	Y = *vvi_iterator;

	tk::spline s;
	s.set_points(X,Y);    // currently it is required that X is already sorted

//	vector<double> X1, Y1;
//	vi_iterator = X.begin();
//	double left = *vi_iterator;
//	double right = 0.0;
//	++vi_iterator;
//	for(; vi_iterator!=X1.end(); ++vi_iterator){
//		right = *vi_iterator;
//		for(double i=left; i<=right; i=i+0.25){
//			X1.push_back(i);
//			double t=s(i);
//			Y1.push_back(t);
//		}
//		left = *vi_iterator;
//	}

	cout << s(9.25) << endl;
//	vector<double>::iterator vi_iterator2;
//	for(vi_iterator = X1.begin(), vi_iterator2=Y1.begin(); vi_iterator!=X1.end(); ++vi_iterator,++vi_iterator2) {
//		cout << *vi_iterator << " " << *vi_iterator2 << endl;
//	}

	//return std::make_pair(X1, Y1);
}
