#include "complexity_analysis.h"


int main()
{
    int bit_N = 1024, umax , ux, u; //suppose P-1 = 2*17*3^u
    long long pt = 1073741789, px; //max 30-bit prime
    vector<double> som22Ts = {}, DSPTs = {};
    for(u = 0; u <= ceil(512*log(2)/log(3)); u+=20){
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
        som22Ts.insert(som22Ts.end(),Ts.first);
        DSPTs.insert(DSPTs.end(),Ts.second);
    }
    cout<<"u range:[0, "<<u<<" ]."<<endl;
    cout<<"som22Ts from estimator: ";
    print_vector(som22Ts,0, som22Ts.size());
    cout<<"DSPTs from estimator: ";
    print_vector(DSPTs,0, DSPTs.size());


    // //ls = [21, 22, 23, 24, 25, 26]
    // for(int l = 21; l <= 26; l++)
    //     FM_complexiy(1, l, prob_prime, true);
    // // delete PrimeList;   //删除保存素数的数组，回收内存

}