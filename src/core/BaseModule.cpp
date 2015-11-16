/*
 * BaseModule.cpp
 *
 *  Created on: Nov 11, 2015
 *      Author: Chandra
 */
#include <cassert>
#include <iostream>
#include <math.h>
#include "BaseModule.h"
#include "spline.h"

BaseModule::BaseModule() {
	// TODO Auto-generated constructor stub

}

BaseModule::~BaseModule() {
	// TODO Auto-generated destructor stub
}


void BaseModule::calculateEt(){
	data.Et.assign(data.t.size(), 0);
	data.Et[0] = exp(data.a[0] * data.t[0]);
	for(int i = 1; i<data.t.size(); ++i){
		data.Et[i]  = exp(data.a[i] * (data.t[i] - data.t[i-1])) * data.Et[i-1];
	}
}

double BaseModule::calculateB(double t, double T){
    int position_t = locate(t);
    int position_T = locate(T);
    double summation = 0;
    for (int i = position_t+1; i <= position_T; i++){
        summation = summation + (data.t[i] - data.t[i-1]) / data.Et[i];
    }
    if(t){
        summation = summation * data.Et[position_t];
    }
//    cout << "calculateB = " << summation << endl;
    return summation;
}


int BaseModule::locate(double t){
    //return index in vector for given value
//	std::cout << t <<" ";
	if (t==0)
		return -1;
	std::map<double,int>::iterator it;
	it = timePos.find(t);
	assert(it!=timePos.end());
//	std::cout << it->second << std::endl;
    return it->second;
}

double BaseModule::calculateVr(double t){
    //Vr always have s=0 in our Problem
    if (t == 0) {
        return 0;
    }
    int position_t = locate(t);
    double Vr = 0;
    for(int i = 0; i<=position_t; ++i){
        Vr += pow(data.Et[i], 2) * pow(data.sigma[i], 2) * (data.t[i] - data.t[i-1]);
    }
    Vr = Vr / pow(data.Et[position_t], 2);
//    cout << "calculateVr = " << Vr << endl;
    return Vr;
}

double BaseModule::calculateVp(double Tf, double Tp){
    //Vp always have t=0 in our Problem
	double t = calculateVr(Tf) * pow(calculateB(Tf, Tp), 2);
//	cout << "calculateVp = " << t << endl;
    return t;
}

double BaseModule::calculateA(double t, double T){
    int position_t = locate(t);
    int position_T = locate(T);
    double P0t=1;
    if(t){
        P0t=data.p[position_t];
        //f0t=data.f[position_t];
    }
    double valB = calculateB(t, T);
    double x = log(data.p[position_T] / P0t) ;
    double y = valB * data.f[position_t+1];
    double z = 0.5 * pow(valB, 2) * calculateVr(t);
    double t1 = x + y - z;
//    cout << "calculateA = " << x << "+" << y << "-" << z << "=" << t1 << endl;
    return t1;
}

//-----------------------


// Returns the erf() of a value (not super precice, but ok)
double BaseModule::erf(double x){
    double y = 1.0 / ( 1.0 + 0.3275911 * x);
    return 1 - (((((
                    + 1.061405429  * y
                    - 1.453152027) * y
                   + 1.421413741) * y
                  - 0.284496736) * y
                 + 0.254829592) * y)
    * exp (-x * x);
}

// Returns the probability of [-inf,x] of a gaussian distribution
double BaseModule::N(double x, double mu, double sigma){
    return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2.))));
}

double BaseModule::ZBP(double Tf, double Tp, double X){
    int position_Tf = locate(Tf);
    int position_Tp = locate(Tp);
    double l = log(data.p[position_Tf] * X / data.p[position_Tp]);
    double v = sqrt(calculateVp(Tf, Tp));
    double dP =  l / v + 0.5 * v;
    double dN = l / v - 0.5 * v;
    return X * data.p[position_Tf] * N(dP, 0, 1) - data.p[position_Tp] * N(dN, 0, 1);
}

//f(r)
double BaseModule::fun_r(const vector<double>& c, const vector<double>& A, const vector<double>& B, double r){
    double ret = 0;
    for(int i = 0; i < c.size(); ++i){
        ret += c[i] * exp(A[i] - B[i] * r);
    }
    return ret - 1;
}

//f'(r)
double BaseModule::fun_dev(const vector<double>& c, const vector<double>& A, const vector<double>& B, double r){
    double d = 0;
    for(int i = 0; i < c.size(); ++i){
        d = d - c[i] * exp(A[i] - B[i] * r) * B[i];
    }
    return d;
}

//Newton's algrithm to solve r*
double BaseModule::solveR(const vector<double>& c, const vector<double>& A, const vector<double>& B){
    double r = 0.5;
    while (fabs(fun_r(c, A, B, r)) > 1E-10) {
        r = r - fun_r(c, A, B, r) / fun_dev(c, A, B, r);
    }
    return r;
}

double BaseModule::pSwaption(double K, double T0, double Tp){
    int position_T0 = locate(T0);
    int position_Tp = locate(Tp);
    vector<double> ci(position_Tp-position_T0);
    vector<double> Ti(position_Tp-position_T0);
    vector<double> Ai(position_Tp-position_T0);
    vector<double> Bi(position_Tp-position_T0);
    for(int i = 0; i < position_Tp - position_T0; ++i){
        Ti[i] = data.t[position_T0 + i + 1];
        ci[i] = K * (Ti[i] - data.t[position_T0 + i]);
        Ai[i]=calculateA(T0, Ti[i]);
        Bi[i]=calculateB(T0, Ti[i]);
//        cout << "Ai=" << Ai[i] << " Bi=" << Bi[i] << " ci=" << ci[i] << " Ti=" << Ti[i] << endl;
    }
    ci.back() += 1;
    double r_star = solveR(ci, Ai, Bi);
    vector<double> Xi(position_Tp - position_T0);
    double price = 0;
    for(int i = 0; i < position_Tp - position_T0; ++i){
        Xi[i] = exp(Ai[i] - Bi[i] * r_star);
        price += ci[i] * ZBP(T0, Ti[i], Xi[i]);
//        cout << "Xi=" << Xi[i] << " price=" << price <<endl;
    }
    cout << "FinalPrice=" <<"\t"<< price << "\t "<< T0 << " \t" << Tp << endl;
    return price;
}


double BaseModule::takeDev(const vector<double>& value, double point, int n){
    //take first or second derivative
    tk::spline s;
    s.set_points(data.t,value);
    if (n==1) {
        return (s(point+0.01)-s(point-0.01))/0.02;
    }
    else{
        return (s(point+0.01)-2*s(point)+s(point-0.01))/0.0001;
    }
}

vector<double> BaseModule::theta(){
    vector<double> V(data.t.size());
    vector<double> theta_vec(data.t.size());
    for(int i = 0; i < data.t.size(); ++i){
        vector<double> sig(i+1);
        for(int j=0; j<=i; j++){
            sig[j]=data.sigma[j] * calculateB(data.t[j], data.t[i]);
        }
        V[i]=pow(sig[0], 2) * data.t[0] ;
        for(int j = 1; j <= i; ++j){
            V[i]= V[i]+pow(sig[j], 2) * (data.t[j]-data.t[j-1]);
        }
    }
    for(int i = 0; i < data.t.size(); ++i){
        theta_vec[i]=takeDev(data.f, data.t[i], 1)
        		+ data.a[i] * data.f[i]
        		+ 0.5 * (takeDev(V, data.t[i], 2) + data.a[i] * takeDev(V, data.t[i], 1));
    }
    return theta_vec;
}
