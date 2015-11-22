/*
 * Module4Appl.h
 *
 *  Created on: Nov 21, 2015
 *      Author: Chandra
 */

#include "core/BaseModule.h"
#include "core/ModuleHelper.h"

#ifndef MODULE4APPL_H_
#define MODULE4APPL_H_

class Module4Appl: public BaseModule{

private:
	ModuleHelper helper = ModuleHelper();
public:
	Module4Appl();
	virtual ~Module4Appl();
	int moduleMainFunc();
};

#endif /* MODULE4APPL_H_ */
