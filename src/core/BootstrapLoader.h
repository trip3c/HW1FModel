/*
 * BootstrapLoader.h
 *
 *  Created on: Nov 26, 2015
 *      Author: Chandra
 */

#include <map>
#include <string>
using namespace std;
#ifndef BOOTSTRAPLOADER_H_
#define BOOTSTRAPLOADER_H_

class BootstrapLoader {
public:
	BootstrapLoader();
	BootstrapLoader(std::string pConfigFile);
	virtual ~BootstrapLoader();

	string getApplicationLoggingLevel();

	double getDefaultsPaymentFrequency();
	double getDefaultsConstantWeights();

	double getMeanReversionConstant();
	double getMeanReversionA0();
	double getMeanReversionA1();
	double getMeanReversionA2();
	double getMeanReversionSDA0();
	double getMeanReversionSDA1();
	double getMeanReversionSDA2();
	double getMeanReversionGamma();
	double getMeanReversionSimulations();

	double getVolatilityConstant();
	double getVolatilityA0();
	double getVolatilityA1();
	double getVolatilityA2();
	double getVolatilityA3();
	double getVolatilitySDA0();
	double getVolatilitySDA1();
	double getVolatilitySDA2();
	double getVolatilitySDA3();
	double getVolatilityGamma();
	double getVolatilitySimulations();

	string getSwaptionVolatility();
	string getSwaptionMaturity();
	string getSwaptionTenor();
	string getSwaptionWeights();
	string getDiscounts();
private:
	std::map<const char*, std::string> filesMapPath;
	std::map<const char*, std::string> applicationMap;
	std::map<const char*, double> volatilityMap;
	std::map<const char*, double> meanReversionMap;
	std::map<const char*, double> defaultsMap;
	void loadConfigs(std::string pConfigFile);
	string getStringValue(map<const char*, std::string> mp, const char* key);
	double getDoubleValue(map<const char*, double> mp, const char* key);
};

#endif /* BOOTSTRAPLOADER_H_ */
