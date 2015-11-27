/*
 * ApplicationHullWhite.h
 *
 *  Created on: Nov 10, 2015
 *      Author: Chandra
 */

#include <string>
#include "core/BootstrapLoader.h"

#ifndef APPLICATIONHULLWHITE_H_
#define APPLICATIONHULLWHITE_H_

class ApplicationHullWhite{
private:
	BootstrapLoader serviceLocator;
public:
	ApplicationHullWhite();
	ApplicationHullWhite(BootstrapLoader loader);

	int moduleChooserFunc();
};

#endif /* APPLICATIONHULLWHITE_H_ */
