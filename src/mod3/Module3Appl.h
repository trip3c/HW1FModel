/*
 * Module3Appl.h
 *
 *  Created on: Nov 15, 2015
 *      Author: Chandra
 */

#include "core/BaseModule.h"
#include "core/ModuleHelper.h"

#ifndef MODULE3APPL_H_
#define MODULE3APPL_H_

class Module3Appl: public BaseModule{

private:
	ModuleHelper helper = ModuleHelper();
public:
	Module3Appl();
	Module3Appl(BootstrapLoader sl);
	virtual ~Module3Appl();
	int moduleMainFunc();
};

#endif /* MODULE3APPL_H_ */
