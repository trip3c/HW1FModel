/*
 * Module2Appl.h
 *
 *  Created on: Nov 15, 2015
 *      Author: Chandra
 */

#include "core/BaseModule.h"
#include "core/ModuleHelper.h"

#ifndef MODULE2APPL_H_
#define MODULE2APPL_H_

class Module2Appl: public BaseModule{

private:
	ModuleHelper helper = ModuleHelper();
public:
	Module2Appl();
	virtual ~Module2Appl();
	int moduleMainFunc();
};

#endif /* MODULE2APPL_H_ */
