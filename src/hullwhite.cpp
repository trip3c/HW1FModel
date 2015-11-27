//============================================================================
// Name        : hullwhite.cpp
// Author      : chandra, james, alan, gloria
// Version     :
// Copyright   : Hull White One Factor Calibration
// Description : Hull White One Factor Calibration
//============================================================================

#include "appl/ApplicationHullWhite.h"
#include "core/log.h"
#include "core/ConfigFile.h"
#include <iostream>
#include <map>
using namespace std;
#include "core/Constants.h"
#include "core/BootstrapLoader.h"

int main(int argc, char * argv[]){
	std::string configFile = "config.txt";
	if (argc > 1){
		configFile = argv[1];
	}
	BootstrapLoader loader(configFile);
	ApplicationHullWhite app(loader);
	FILELog::ReportingLevel() = FILELog::FromString(loader.getApplicationLoggingLevel());
	FILE* log_fd = fopen( "output.txt", "w" );
	Output2FILE::Stream() = log_fd;
	app.moduleChooserFunc();
	return 0;
}
