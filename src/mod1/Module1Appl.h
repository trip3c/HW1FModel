/*
 * Module1Appl.h
 *
 *  Created on: Nov 10, 2015
 *      Author: Chandra
 */

#include "core/BaseModule.h"
#include "core/ModuleHelper.h"

#ifndef MODULE1APPL_H_
#define MODULE1APPL_H_

class Module1Appl: public BaseModule{

private:
	ModuleHelper helper;
public:
	Module1Appl();
	Module1Appl(BootstrapLoader sl);
	virtual ~Module1Appl();
	int moduleMainFunc();
};

#endif /* MODULE1APPL_H_ */
