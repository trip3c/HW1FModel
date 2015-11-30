/*
 * BootstrapLoader.cpp
 *
 *  Created on: Nov 26, 2015
 *      Author: Chandra
 */

#include <core/BootstrapLoader.h>
#include <string>
using namespace std;
#include "core/Constants.h"
#include "core/ConfigFile.h"

BootstrapLoader::BootstrapLoader() {}

BootstrapLoader::BootstrapLoader(std::string pConfigFile) {
	loadConfigs(pConfigFile);
}

BootstrapLoader::~BootstrapLoader() {}

void BootstrapLoader::loadConfigs(std::string configFile){
	ConfigFile cf(configFile);

	applicationMap[Constants::CONFIG_APPLICATION_LOGGING_LEVEL] =
		(std::string) cf.Value(Constants::CONFIG_SECTIONS_APPLICATION, Constants::CONFIG_APPLICATION_LOGGING_LEVEL);

	defaultsMap[Constants::CONFIG_DEFAULTS_CONSTANT_WEIGHT] =
			cf.Value(Constants::CONFIG_SECTIONS_DEFAULTS, Constants::CONFIG_DEFAULTS_CONSTANT_WEIGHT);
	defaultsMap[Constants::CONFIG_DEFAULTS_PAYMENT_FREQ] =
			cf.Value(Constants::CONFIG_SECTIONS_DEFAULTS, Constants::CONFIG_DEFAULTS_PAYMENT_FREQ);

	meanReversionMap[Constants::CONFIG_MEANREVERSION_A0] =
			cf.Value(Constants::CONFIG_SECTIONS_MEANREVERSION, Constants::CONFIG_MEANREVERSION_A0);
	meanReversionMap[Constants::CONFIG_MEANREVERSION_CONSTANT] =
			cf.Value(Constants::CONFIG_SECTIONS_MEANREVERSION, Constants::CONFIG_MEANREVERSION_CONSTANT);
	meanReversionMap[Constants::CONFIG_MEANREVERSION_A1] =
			cf.Value(Constants::CONFIG_SECTIONS_MEANREVERSION, Constants::CONFIG_MEANREVERSION_A1);
	meanReversionMap[Constants::CONFIG_MEANREVERSION_A2] =
			cf.Value(Constants::CONFIG_SECTIONS_MEANREVERSION, Constants::CONFIG_MEANREVERSION_A2);
	meanReversionMap[Constants::CONFIG_MEANREVERSION_SD_A0] =
			cf.Value(Constants::CONFIG_SECTIONS_MEANREVERSION, Constants::CONFIG_MEANREVERSION_SD_A0);
	meanReversionMap[Constants::CONFIG_MEANREVERSION_SD_A1] =
			cf.Value(Constants::CONFIG_SECTIONS_MEANREVERSION, Constants::CONFIG_MEANREVERSION_SD_A1);
	meanReversionMap[Constants::CONFIG_MEANREVERSION_SD_A2] =
			cf.Value(Constants::CONFIG_SECTIONS_MEANREVERSION, Constants::CONFIG_MEANREVERSION_SD_A2);
	meanReversionMap[Constants::CONFIG_MEANREVERSION_GAMMA] =
			cf.Value(Constants::CONFIG_SECTIONS_MEANREVERSION, Constants::CONFIG_MEANREVERSION_GAMMA);
	meanReversionMap[Constants::CONFIG_MEANREVERSION_SIMULATIONS] =
			cf.Value(Constants::CONFIG_SECTIONS_MEANREVERSION, Constants::CONFIG_MEANREVERSION_SIMULATIONS);

	volatilityMap[Constants::CONFIG_VOLATILITY_CONSTANT] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_CONSTANT);
	volatilityMap[Constants::CONFIG_VOLATILITY_A0] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_A0);
	volatilityMap[Constants::CONFIG_VOLATILITY_A1] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_A1);
	volatilityMap[Constants::CONFIG_VOLATILITY_A2] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_A2);
	volatilityMap[Constants::CONFIG_VOLATILITY_A3] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_A3);
	volatilityMap[Constants::CONFIG_VOLATILITY_SD_A0] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_SD_A0);
	volatilityMap[Constants::CONFIG_VOLATILITY_SD_A1] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_SD_A1);
	volatilityMap[Constants::CONFIG_VOLATILITY_SD_A2] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_SD_A2);
	volatilityMap[Constants::CONFIG_VOLATILITY_SD_A3] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_SD_A3);
	volatilityMap[Constants::CONFIG_VOLATILITY_GAMMA] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_GAMMA);
	volatilityMap[Constants::CONFIG_VOLATILITY_SIMULATIONS] =
			cf.Value(Constants::CONFIG_SECTIONS_VOLATILITY, Constants::CONFIG_VOLATILITY_SIMULATIONS);

	filesMapPath[Constants::CONFIG_SWAPTION_MATURITY] =
		(std::string) cf.Value(Constants::CONFIG_SECTIONS_FILES, Constants::CONFIG_SWAPTION_MATURITY);
	filesMapPath[Constants::CONFIG_SWAPTION_TENOR] =
		(std::string) cf.Value(Constants::CONFIG_SECTIONS_FILES, Constants::CONFIG_SWAPTION_TENOR);
	filesMapPath[Constants::CONFIG_SWAPTION_VOLATILITY] =
		(std::string) cf.Value(Constants::CONFIG_SECTIONS_FILES, Constants::CONFIG_SWAPTION_VOLATILITY);
	filesMapPath[Constants::CONFIG_SWAPTION_WEIGHTS] =
		(std::string) cf.Value(Constants::CONFIG_SECTIONS_FILES, Constants::CONFIG_SWAPTION_WEIGHTS);
	filesMapPath[Constants::CONFIG_DISCOUNTS] =
		(std::string) cf.Value(Constants::CONFIG_SECTIONS_FILES, Constants::CONFIG_DISCOUNTS);
}

string BootstrapLoader::getStringValue(map<const char*, std::string> mp, const char* key){
	auto search = mp.find(key);
	if (search!=mp.end())
		return search->second;
	return "";
}

double BootstrapLoader::getDoubleValue(map<const char*, double> mp, const char* key){
	auto search = mp.find(key);
	if (search!=mp.end())
		return search->second;
	return 0.0;
}

string BootstrapLoader::getApplicationLoggingLevel(){
	return getStringValue(applicationMap, Constants::CONFIG_APPLICATION_LOGGING_LEVEL);
}

double BootstrapLoader::getDefaultsPaymentFrequency(){
	return getDoubleValue(defaultsMap, Constants::CONFIG_DEFAULTS_PAYMENT_FREQ);
}
double BootstrapLoader::getDefaultsConstantWeights(){
	return getDoubleValue(defaultsMap, Constants::CONFIG_DEFAULTS_CONSTANT_WEIGHT);
}

double BootstrapLoader::getMeanReversionConstant(){
	return getDoubleValue(meanReversionMap, Constants::CONFIG_MEANREVERSION_CONSTANT);
}
double BootstrapLoader::getMeanReversionA0(){
	return getDoubleValue(meanReversionMap, Constants::CONFIG_MEANREVERSION_A0);
}
double BootstrapLoader::getMeanReversionA1(){
	return getDoubleValue(meanReversionMap, Constants::CONFIG_MEANREVERSION_A1);
}
double BootstrapLoader::getMeanReversionA2(){
	return getDoubleValue(meanReversionMap, Constants::CONFIG_MEANREVERSION_A2);
}
double BootstrapLoader::getMeanReversionSDA0(){
	return getDoubleValue(meanReversionMap, Constants::CONFIG_MEANREVERSION_SD_A0);
}
double BootstrapLoader::getMeanReversionSDA1(){
	return getDoubleValue(meanReversionMap, Constants::CONFIG_MEANREVERSION_SD_A1);
}
double BootstrapLoader::getMeanReversionSDA2(){
	return getDoubleValue(meanReversionMap, Constants::CONFIG_MEANREVERSION_SD_A2);
}
double BootstrapLoader::getMeanReversionGamma(){
	return getDoubleValue(meanReversionMap, Constants::CONFIG_MEANREVERSION_GAMMA);
}
double BootstrapLoader::getMeanReversionSimulations(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_MEANREVERSION_SIMULATIONS);
}

double BootstrapLoader::getVolatilityConstant(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_CONSTANT);
}
double BootstrapLoader::getVolatilityA0(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_A0);
}
double BootstrapLoader::getVolatilityA1(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_A1);
}
double BootstrapLoader::getVolatilityA2(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_A2);
}
double BootstrapLoader::getVolatilityA3(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_A3);
}
double BootstrapLoader::getVolatilitySDA0(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_SD_A0);
}
double BootstrapLoader::getVolatilitySDA1(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_SD_A1);
}
double BootstrapLoader::getVolatilitySDA2(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_SD_A2);
}
double BootstrapLoader::getVolatilitySDA3(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_SD_A3);
}
double BootstrapLoader::getVolatilityGamma(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_GAMMA);
}
double BootstrapLoader::getVolatilitySimulations(){
	return getDoubleValue(volatilityMap, Constants::CONFIG_VOLATILITY_SIMULATIONS);
}

string BootstrapLoader::getSwaptionVolatility(){
	return getStringValue(filesMapPath, Constants::CONFIG_SWAPTION_VOLATILITY);
}
string BootstrapLoader::getSwaptionMaturity(){
	return getStringValue(filesMapPath, Constants::CONFIG_SWAPTION_MATURITY);
}
string BootstrapLoader::getSwaptionTenor(){
	return getStringValue(filesMapPath, Constants::CONFIG_SWAPTION_TENOR);
}
string BootstrapLoader::getSwaptionWeights(){
	return getStringValue(filesMapPath, Constants::CONFIG_SWAPTION_WEIGHTS);
}
string BootstrapLoader::getDiscounts(){
	return getStringValue(filesMapPath, Constants::CONFIG_DISCOUNTS);
}
