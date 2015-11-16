//============================================================================
// Name        : hullwhite.cpp
// Author      : chandra, james, alan, gloria
// Version     :
// Copyright   : Hull White One Factor Calibration
// Description : Hull White One Factor Calibration
//============================================================================

#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "spline.h"
#include "appl\ApplicationHullWhite.h"

int main(int argc, char * argv[]){
	ApplicationHullWhite app("E:\\work\\cpp_ws\\");
	app.moduleChooserFunc();
	return 0;
}
