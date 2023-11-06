//#include "pollard_pm1.h"
#include "utils.h"
// #include <bitset>


double prob_prime(int i);
double prob1(int i);
double normal_pow_complexity(int bit_a,int bit_pow);
pair<double,long long> FM_complexity(int n1, int n2, double (*pr)(int), bool verbose = false, double C1 = 0.5);
double slide_window_MP_cost( int bit_a, long long bit_pow, double C1 =0.5);
// double method1_complexity(int bit_N, int bit_a, long long &CntB,long long BList[],bool fix_bit);
// double method2_complexity(int bit_N, int bit_a, long long &CntB,long long BList[],bool fix_bit);
// void complexity_test(int bit_N, int bit_a, long long &CntB,long long BList[],bool fix_bit,int block);
double gen_prime_list_complexity(int B);
double original_pollard_p1_complexity(int bit_N, int bit_a, long long &CntB,long long BList[],bool fix_bit);
pair<double,double> complexity_test(int bit_N, long long pt , int umax, long long px, int ux);
double GCD_cost(long long n, double C2 = 0.5);