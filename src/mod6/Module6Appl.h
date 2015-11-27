/*
 * Module6Appl.h
 *
 *  Created on: Nov 24, 2015
 *      Author: Chandra
 */

#include "core/BaseModule.h"
#include "core/ModuleHelper.h"

#ifndef MODULE6APPL_H_
#define MODULE6APPL_H_

class Module6Appl: public BaseModule{

private:
	ModuleHelper helper = ModuleHelper();
public:
	Module6Appl();
	Module6Appl(BootstrapLoader sl);
	virtual ~Module6Appl();
	int moduleMainFunc();
};

#endif /* MODULE6APPL_H_ */
