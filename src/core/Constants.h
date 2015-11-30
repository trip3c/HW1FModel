/*
 * Constants.h
 *
 *  Created on: Nov 15, 2015
 *      Author: Chandra
 */

#include <string>

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

namespace Constants {
//	const double PAYMENT_FREQ = 0.25;
//	static constexpr double FIXED_MEAN_REVERSION = 5.0/100;
//	static constexpr double FIXED_VOLATILITY = 1.5/100;
//	const double FIXED_STRIKE = 4.0/100;
//	const double CONSTANT_WEIGHT = 1;

	static constexpr const char* CONFIG_SECTIONS_APPLICATION = "Application";
	static constexpr const char* CONFIG_SECTIONS_DEFAULTS = "Defaults";
	static constexpr const char* CONFIG_SECTIONS_MEANREVERSION = "MeanReversion";
	static constexpr const char* CONFIG_SECTIONS_VOLATILITY = "Volatility";
	static constexpr const char* CONFIG_SECTIONS_FILES = "Files";

	static constexpr const char* CONFIG_APPLICATION_LOGGING_LEVEL = "LoggingLevel";

	static constexpr const char* CONFIG_DEFAULTS_PAYMENT_FREQ = "Payment_Freq";
	static constexpr const char* CONFIG_DEFAULTS_CONSTANT_WEIGHT = "constant_weight";

	static constexpr const char* CONFIG_MEANREVERSION_CONSTANT = "constant";
	static constexpr const char* CONFIG_MEANREVERSION_A0 = "a0";
	static constexpr const char* CONFIG_MEANREVERSION_A1 = "a1";
	static constexpr const char* CONFIG_MEANREVERSION_A2 = "a2";
	static constexpr const char* CONFIG_MEANREVERSION_SD_A0 = "sd_a0";
	static constexpr const char* CONFIG_MEANREVERSION_SD_A1 = "sd_a1";
	static constexpr const char* CONFIG_MEANREVERSION_SD_A2 = "sd_a2";
	static constexpr const char* CONFIG_MEANREVERSION_GAMMA = "gamma";
	static constexpr const char* CONFIG_MEANREVERSION_SIMULATIONS = "simulations";

	static constexpr const char* CONFIG_VOLATILITY_CONSTANT = "constant";
	static constexpr const char* CONFIG_VOLATILITY_A0 = "a0";
	static constexpr const char* CONFIG_VOLATILITY_A1 = "a1";
	static constexpr const char* CONFIG_VOLATILITY_A2 = "a2";
	static constexpr const char* CONFIG_VOLATILITY_A3 = "a3";
	static constexpr const char* CONFIG_VOLATILITY_SD_A0 = "sd_a0";
	static constexpr const char* CONFIG_VOLATILITY_SD_A1 = "sd_a1";
	static constexpr const char* CONFIG_VOLATILITY_SD_A2 = "sd_a2";
	static constexpr const char* CONFIG_VOLATILITY_SD_A3 = "sd_a3";
	static constexpr const char* CONFIG_VOLATILITY_GAMMA = "gamma";
	static constexpr const char* CONFIG_VOLATILITY_SIMULATIONS = "simulations";

	static constexpr const char* CONFIG_SWAPTION_VOLATILITY = "swaption_volatility";
	static constexpr const char* CONFIG_SWAPTION_MATURITY = "swaption_maturity";
	static constexpr const char* CONFIG_SWAPTION_TENOR = "swaption_tenor";
	static constexpr const char* CONFIG_SWAPTION_WEIGHTS = "swaption_weights";
	static constexpr const char* CONFIG_DISCOUNTS = "discounts";

};


#endif /* CONSTANTS_H_ */
