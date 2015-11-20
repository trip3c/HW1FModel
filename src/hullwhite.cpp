//============================================================================
// Name        : hullwhite.cpp
// Author      : chandra, james, alan, gloria
// Version     :
// Copyright   : Hull White One Factor Calibration
// Description : Hull White One Factor Calibration
//============================================================================

#include "appl/ApplicationHullWhite.h"
#include "core/log.h"

int main(int argc, char * argv[]){
	ApplicationHullWhite app("E:\\work\\cpp_ws\\");
	FILELog::ReportingLevel() = FILELog::FromString(argv[1] ? argv[1] : "INFO");
	FILE* log_fd = fopen( "output.txt", "w" );
	Output2FILE::Stream() = log_fd;
	app.moduleChooserFunc();
	// TODO Close file
	return 0;
}
