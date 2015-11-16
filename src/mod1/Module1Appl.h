/*
 * Module1Appl.h
 *
 *  Created on: Nov 10, 2015
 *      Author: Chandra
 */
#include "..\core\FileUtils.h"
#include "..\core\BaseModule.h"

#ifndef MODULE1APPL_H_
#define MODULE1APPL_H_

class Module1Appl: public BaseModule{

private:
	FileUtils fileUtils;
	void initializeDataStructure();
	void interpolatePricesAndForwardRates(vector<vector<double> > newNumbers);
	void assignMeanReversion();
	void assignVol(double s);
	void printStructure();
public:
	Module1Appl();
	virtual ~Module1Appl();
	int moduleMainFunc();
};

#endif /* MODULE1APPL_H_ */
