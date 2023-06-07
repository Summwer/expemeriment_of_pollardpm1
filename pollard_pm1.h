#include <gmp.h>
#include <string>
#include "generate_primes.h"
#include "utils.h"


// Z_NR<mpz_t> fast_mul1(vector<Z_NR<mpz_t>> List);
Z_NR<mpz_t> fast_mul(vector<Z_NR<mpz_t>> List);
Z_NR<mpz_t> normal_mul(vector<Z_NR<mpz_t>> List);


void pollard_pm1_Pol74(mpz_t N, int bits, ll L, ll M, bool verbose);
void pollard_pm1_ref26(mpz_t N,int bits, ll &CntPrime,ll PrimeList[]);
void pollard_pm1_Bis03(mpz_t N, int bits,bool verbose);
void pollard_pm1_IPP1_V1(mpz_t N, int bits,ll &CntPrime,ll PrimeList[]);
void pollard_pm1_IPP1_V2(mpz_t N, int bits,ll &CntPrime,ll PrimeList[]);
void pollard_pm1_improved_IPP1_V2(mpz_t N, int bits, ll MaxInt, bool fm, bool verbose);



void generate_partial_primes(vector<Z_NR<mpz_t>> Primes, int k, vector<vector<Z_NR<mpz_t>>> &Partial_Primes);
void dynamic_scaling_pollard_pm1(mpz_t N, int bit_N, ll MaxInt,bool fm, bool verbose);
// void dynamic_scaling_pollard_pm12(mpz_t N, int bit_N, ll MaxInt,bool fm, bool verbose);



void dynamic_scaling_pollard_pm1_with_block_partition(mpz_t N, int bits,ll &CntPrime,ll PrimeList[], unsigned long D, int k, int blocksize);
