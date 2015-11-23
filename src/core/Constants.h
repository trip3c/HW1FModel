/*
 * Constants.h
 *
 *  Created on: Nov 15, 2015
 *      Author: Chandra
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

struct Constants {
	const double PAYMENT_FREQ = 0.25;
	static constexpr double FIXED_MEAN_REVERSION = 5.0/100;
	static constexpr double FIXED_VOLATILITY = 17.9/100;
//	const double FIXED_STRIKE = 4.0/100;
	const double CONSTANT_WEIGHT = 1;
};


#endif /* CONSTANTS_H_ */
