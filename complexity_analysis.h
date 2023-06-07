#include "pollard_pm1.h"
#include "utils.h"
// #include <bitset>

double normal_pow_complexity(int bit_a,int bit_pow);
pair<double,long long> FMP_complexiy(int n1, int n2, bool verbose = false);
double slide_window_pow_mod_complexity(int bit_a,int bit_pow);
// double method1_complexity(int bit_N, int bit_a, ll &CntB,ll BList[],bool fix_bit);
// double method2_complexity(int bit_N, int bit_a, ll &CntB,ll BList[],bool fix_bit);
// void complexity_test(int bit_N, int bit_a, ll &CntB,ll BList[],bool fix_bit,int block);
double gen_prime_list_complexity(int B);
double original_pollard_p1_complexity(int bit_N, int bit_a, ll &CntB,ll BList[],bool fix_bit);
void complexity_test(int bit_N, int pt , int umax, int px, int ux);