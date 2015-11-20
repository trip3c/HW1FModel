
# include <vector>
# include <time.h>
#include "MultiNormalDist.h"
#include <random>


using namespace std;

//Example of writing functor of F(a)
class func{
public:
    //calculate F(a)
    //parameter x should be vector of A0~A3
    double operator()(vector<double> x){
        //first calculate a(t) from A0~A3
        //calculate F(a(t))
        return 0;
    };
    
    //test whether A0~A3 make a_max<a(t)<a_min
    bool satisfyConstrain(vector<double> x){
        return true;
    }
    
};


//here F can be either F(a) or G(sigma)
//x0 is initial guess, N is the simulation time, r is the cool-down speed
//template<class F>
//vector<double> Calibrate(vector<double> x0, F f, double N, double r){
//    vector<double> s(x0.size() * x0.size());
//    //for now, we set the initial var=x0, and cov=0
//    for(int i=0;i<x0.size();++i){
//        for(int j=0;j<x0.size();++j){
//            if(j==i){
//                s[j+x0.size() * i]=fabs(x0[i]);
//            }
//            else{
//                s[j+x0.size() * i]=0;
//            }
//        }
//    }
//    for(int i=0;i<N;++i){
////    	srand (time(NULL));
////        int seed=rand();
//        time_t systime;                                // solutions of the three runs will be chosen as the final solution
//		time(&systime);
//		srand((unsigned int)systime);
//        int seed=rand();
//
//        vector<double> x1=multinormal_sample ( x0.size(), 1, s, x0, seed );
//        while(!f.satisfyConstrain(x1)){
//            x1=multinormal_sample ( x0.size(), 1, s, x0, seed );
//        }
//        if(f(x1)<f(x0)){
//            x0=x1;
//            //shrink the cov each step forward
//            for(int j=0;j<x0.size();++j){
//                s[j+x0.size() * j]=s[j+x0.size() * j]*exp(-r*1/(N-1));
//            }
//        }
//    }
//    return x0;
//};


