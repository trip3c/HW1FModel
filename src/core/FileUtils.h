/*
 * FileUtils.h
 *
 *  Created on: Nov 10, 2015
 *      Author: Chandra
 */
#include <vector>
#include <string>
using namespace std;

#ifndef FILEUTILS_H_
#define FILEUTILS_H_

class FileUtils {
private:
	void printVectorOfVectors(vector<vector<double> > vec);
public:
	FileUtils();
	virtual ~FileUtils();
	vector<vector<double> > importData(string filePath);
	vector<vector<double> > rowsToColumnTransposeVector(vector<vector<double> > numbers);
	string getDirectoryPath();
//	void importSwaptionVolatilityOnly(string filePath);
};

#endif /* FILEUTILS_H_ */
