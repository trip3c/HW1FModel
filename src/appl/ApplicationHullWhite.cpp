#include <iostream>
#include <string>
#include "ApplicationHullWhite.h"
#include "..\mod1\Module1Appl.h"

using namespace std;

ApplicationHullWhite::ApplicationHullWhite(){
	ApplicationHullWhite("");
}

ApplicationHullWhite::ApplicationHullWhite(string sDefPath){
	defPath = sDefPath;
}

int ApplicationHullWhite::moduleChooserFunc(){
	cout << defPath << endl;

	char choice;
	do{
		cout << "Choose Module to execute: " << endl;
		cout << "Module 1: Press 1" << endl;
		cout << "Module 2: Press 2" << endl;
		cout << "Module 3: Press 3" << endl;
		cout << "Module 4: Press 4" << endl;
		cout << "Exit:     Press q" << endl;

		cin >> choice;
		if(choice=='1'){
			cout<< "Module 1 starting... " << endl;
			Module1Appl mod1 = Module1Appl();
			mod1.moduleMainFunc();
		}else if (choice=='2'){
			cout<< "choice 2 " << endl;
		}else if (choice=='3'){
			cout<< "choice 3 " << endl;
		}else if (choice=='4'){
			cout<< "choice 4 " << endl;
		}
	} while(choice == '1' || choice == '2' || choice == '3' || choice == '4');
	return 0;
}
