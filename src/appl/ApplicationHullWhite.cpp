#include <iostream>
#include <string>
#include "ApplicationHullWhite.h"
#include "mod1/Module1Appl.h"
#include "mod2/Module2Appl.h"
#include "mod3/Module3Appl.h"
#include "mod4/Module4Appl.h"
#include "mod5/Module5Appl.h"
#include "core/log.h"

using namespace std;

ApplicationHullWhite::ApplicationHullWhite(){
	ApplicationHullWhite("");
}

ApplicationHullWhite::ApplicationHullWhite(string sDefPath){
	defPath = sDefPath;
}

int ApplicationHullWhite::moduleChooserFunc(){
	cout << defPath << endl;
    FILE_LOG(logDEBUG) << "Testing logger " << 12 << " iterations";

	char choice;
	cout << "Choose Module to execute: " << endl;
	cout << "Press 1 for Module 1: Pricing European Swaption" << endl;
	cout << "Press 2 for Module 2: Calibrating theta to yield curve" << endl;
	cout << "Press 3 for Module 3: Calibrating mean reversion" << endl;
	cout << "Press 4 for Module 4: Calculating Black swaption price" << endl;
	cout << "Press 5 for Module 5: Calibrating volatility" << endl;
	cout << "Exit:     Press q" << endl;

	cin >> choice;
	if(choice=='1'){
		cout<< "Module 1 starting... " << endl;
		Module1Appl mod1 = Module1Appl();
		mod1.moduleMainFunc();
	}else if (choice=='2'){
		cout<< "Module 2 starting... " << endl;
		Module2Appl mod2 = Module2Appl();
		mod2.moduleMainFunc();
	}else if (choice=='3'){
		cout<< "Module 3 starting... " << endl;
		Module3Appl mod3 = Module3Appl();
		mod3.moduleMainFunc();
	}else if (choice=='4'){
		cout<< "Module 4 starting... " << endl;
		Module4Appl mod4 = Module4Appl();
		mod4.moduleMainFunc();
	}else if (choice=='5'){
		cout<< "Module 5 starting... " << endl;
		Module5Appl mod5 = Module5Appl();
		mod5.moduleMainFunc();
	}
	return 0;
}
