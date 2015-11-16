/*
 * FileUtils.cpp
 *
 *  Created on: Nov 10, 2015
 *      Author: Chandra
 */
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cassert>
#include "FileUtils.h"

#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif

using namespace std;
FileUtils::FileUtils() {
	// TODO Auto-generated constructor stub

}

FileUtils::~FileUtils() {
	// TODO Auto-generated destructor stub
}

vector<vector<double> > FileUtils::importData(string filePath){
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

vector<vector<double> > FileUtils::rowsToColumnTransposeVector(vector<vector<double> > numbers){
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

string FileUtils::getDirectoryPath(){
	char cCurrentPath[FILENAME_MAX];
	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
	{
		return "";
	}
	cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */
	return cCurrentPath;
}

void FileUtils::printVectorOfVectors(vector<vector<double> > vec){
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
