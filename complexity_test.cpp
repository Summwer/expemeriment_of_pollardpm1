#include "complexity_analysis.h"


int main(int argc, char **argv)
{
    //Cost estimator for IPP1 and our variant
    if(atoi(argv[1])==1){
        int bit_N = 1024, umax , ux, u; //suppose P-1 = 2*17*3^u
        long long pt = 1073741789, px; //max 30-bit prime
        vector<double> IPP1v2Ts = {}, DSPTs = {};
        for(u = 0; u <= min(100,ceil(512*log(2)/log(3))); u+=20){
            if(pt > pow(3,u)){
                px = pt;
                ux = 1;
                umax = max(1,u);
            }
            else{
                px = 3;
                ux = u;
                umax = u;
            }
            // cout<<"==============================="<<endl;
            // cout<<"u = "<<u<<endl;
            pair<double,double> Ts = complexity_test(bit_N, pt , umax, px, ux);
            IPP1v2Ts.insert(IPP1v2Ts.end(),Ts.first);
            DSPTs.insert(DSPTs.end(),Ts.second);
        }
        cout<<"u range:[0, "<<u<<" ]."<<endl;
        cout<<"IPP1v2Ts from estimator: ";
        print_vector(IPP1v2Ts,0, IPP1v2Ts.size());
        cout<<"DSPTs from estimator: ";
        print_vector(DSPTs,0, DSPTs.size());
    }

    //Cost estimator for fast multiplication
    if(atoi(argv[1])==2){
        for(int l = 21; l <= 26; l++){
            FM_complexity(1, l, prob_prime,true);
        }
    }

    // delete PrimeList;   //删除保存素数的数组，回收内存

}