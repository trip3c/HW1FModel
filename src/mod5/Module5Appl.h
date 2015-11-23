/*
 * Module5Appl.h
 *
 *  Created on: Nov 22, 2015
 *      Author: Chandra
 */

#include "core/BaseModule.h"
#include "core/ModuleHelper.h"

#ifndef MODULE5APPL_H_
#define MODULE5APPL_H_

class Module5Appl: public BaseModule{

private:
	ModuleHelper helper = ModuleHelper();
public:
	Module5Appl();
	virtual ~Module5Appl();
	int moduleMainFunc();
};

#endif /* MODULE5APPL_H_ */
