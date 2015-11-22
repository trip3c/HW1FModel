/*
 * BlacksFormula.h
 *
 *  Created on: Nov 21, 2015
 *      Author: Chandra
 */
#include <vector>

#ifndef BLACKSFORMULA_H_
#define BLACKSFORMULA_H_

class BlacksFormula{

public:
	BlacksFormula();
	virtual ~BlacksFormula();
	double normalDistribution(double x);
	double BlacksFormulaFunc(double f, double P, double L, double Rcap, double vol, double tau, double dtau, bool isCaplet);
	std::vector<double> priceBlackCap(std::vector<double> capVol, std::vector<double> PDB,
			std::vector<double> maturity, double Rcap, double L, double tenor, bool isCaplet);
};

#endif /* BLACKSFORMULA_H_ */
