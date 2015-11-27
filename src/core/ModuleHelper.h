/*
 * ModuleHelper.h
 *
 *  Created on: Nov 15, 2015
 *      Author: Chandra
 */

#include <vector>
#include <map>
#include <string>
#include "BootstrapLoader.h"

#ifndef MODULEHELPER_H_
#define MODULEHELPER_H_

class ModuleHelper {
private:
	BootstrapLoader serviceLocator;
	void printVectorOfVectors(std::vector<std::vector<double> > vec);
public:
	ModuleHelper();
	ModuleHelper(BootstrapLoader pLoader);
	virtual ~ModuleHelper();
	std::map<char, std::vector<double> > initializeDataStructure();
	std::vector<std::vector<double> > initializeVolatility();
	std::vector<double> initializeSwaptionVolatilityTenor();
	std::vector<double> initializeSwaptionVolatilityMaturity();

	std::map<char, std::vector<double> > interpolatePricesAndForwardRates(std::vector<std::vector<double> > newNumbers);
	std::vector<std::vector<double> > importDataFromFile(std::string filePath);
	std::vector<std::vector<double> > rowsToColumnTransposeVector(std::vector<std::vector<double> > numbers);
	std::string getDirectoryPath();

	enum structNames {
		FORWARD = 'f',
		TIME = 't',
		PRICE = 'p',
		VOLATILITY = 'v',
		MATURITY = 'm',
		TENOR = 'n'
	};
};

#endif /* MODULEHELPER_H_ */
