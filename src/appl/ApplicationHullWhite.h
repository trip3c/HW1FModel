/*
 * ApplicationHullWhite.h
 *
 *  Created on: Nov 10, 2015
 *      Author: Chandra
 */

#include <string>

#ifndef APPLICATIONHULLWHITE_H_
#define APPLICATIONHULLWHITE_H_

class ApplicationHullWhite{
private:
	std::string defPath;
public:
	ApplicationHullWhite();
	ApplicationHullWhite(std::string sDefPath);

	int moduleChooserFunc();
};

#endif /* APPLICATIONHULLWHITE_H_ */
