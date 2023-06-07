#include <gmp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <fplll.h>

using namespace std;
using namespace fplll;

//long long型对应的最大整数为9223372036854775807
#define ll long long

void CalPrime(ll MaxInt,ll &CntPrime,ll PrimeList[]);
void CalPrime_in_parallel(ll MaxInt,ll &CntPrime,ll PrimeList[], int threads);