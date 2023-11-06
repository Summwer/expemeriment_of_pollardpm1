#include "pollard_pm1.h"
#include <fplll/threadpool.h>


void gen_random_int(mpz_t &p, int bit_num){
    clock_t time = clock();
    gmp_randstate_t grt;
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt, time);
    mpz_urandomb(p, grt, bit_num - 1);
    mpz_t tmp;
    mpz_init_set_ui(tmp,2);
    mpz_pow_ui(tmp,tmp,bit_num-1);
    mpz_add(p,p,tmp);
    mpz_clear(tmp);
}


void gen_prime(mpz_t &p,int bit_num){
    gen_random_int(p,bit_num - 1);
    mpz_nextprime(p, p);
}

//input: 
// :param p: the output
// :param p_bit_num: bit size of p
// :param factor_bit_num: the bit size of small prime factor. factor_bit_num|p_bit_num
//output: prime P such that P-1 = prod(pi), bit(pi) = bit_num.
void gen_prime_with_small_prime_factor(mpz_t &p,  int p_bit_num, int factor_bit_num){
    int amount =(int)floor(p_bit_num / factor_bit_num);
    int remain_bit_size = p_bit_num - factor_bit_num*(amount);
    mpz_t tmp;
    mpz_init(tmp);
    do{
        mpz_set_ui(p,1);
        // cout<<"p_bit_num = "<<p_bit_num<<", amount = "<<amount<<", remain_bit_size = "<<remain_bit_size<<endl;
        // cout<<factor_bit_num<<endl;
        for(int i = 0; i < amount; i++){
            // gen_prime(tmp,factor_bit_num);
            gen_random_int(tmp, factor_bit_num);
            mpz_mul(p,p,tmp);
        }
        if(remain_bit_size > 1){
            gen_random_int(tmp, remain_bit_size);
            mpz_mul(p,p,tmp);
        }
        mpz_add_ui(p,p,1);

        mpz_ui_pow_ui(tmp,2,p_bit_num - 1);
        mpz_cdiv_q(tmp, p, tmp); //tmp = ceil(p/c)
        if(mpz_cmp_ui(tmp,2) == -1)
            remain_bit_size+=1; //generate a larger p with bit size p_bit_num
        else if(mpz_cmp_ui(tmp,2) == 1)
            remain_bit_size-=1;
    }while( (not mpz_probab_prime_p(p, 100) ) or mpz_cmp_ui(tmp,2));
    mpz_clear(tmp);
}


//p  = 2 * pt * 3^u * p_1 * ··· * p_s
void gen_prime_for_high_index(mpz_t &p,  int p_bit_num, int factor_bit_num, int ind3, ll pt){
    int pre_bit;
    
    if(ind3 == 0)
        pre_bit = p_bit_num - floor(log2(pt)) - 2; //bit_p - bit(2 * pt)
    else
        pre_bit = p_bit_num - floor(log2(pt)+ind3 * log2(3)) - 2; //bit_p - bit(2 * pt * 3^u)
    int amount =(int)floor( pre_bit / factor_bit_num);
    int remain_bit_size = pre_bit - factor_bit_num*(amount);
    mpz_t tmp;
    mpz_init(tmp);
    do{
        mpz_ui_pow_ui(p,3,ind3);    
        mpz_mul_ui(p,p, 2 * pt);
        for(int i = 0; i < amount; i++){
            gen_prime(tmp,factor_bit_num);
            // gen_random_int(tmp, factor_bit_num);
            mpz_mul(p,p,tmp);
        }
        if(remain_bit_size > 1){
            gen_prime(tmp,remain_bit_size);
            // gen_random_int(tmp, remain_bit_size);
            mpz_mul(p,p,tmp);
        }
        else if(remain_bit_size==1)
            mpz_mul_ui(p,p,2);

        mpz_add_ui(p,p,1);
        mpz_ui_pow_ui(tmp,2,p_bit_num - 1);
        mpz_cdiv_q(tmp, p, tmp); //tmp = ceil(p/pow(2,bit_num))
        if(mpz_cmp_ui(tmp,2) == -1)
            remain_bit_size+=1; //generate a larger p with bit size p_bit_num
        else if(mpz_cmp_ui(tmp,2) == 1)
            remain_bit_size-=1;
        
        if(remain_bit_size == factor_bit_num){
            remain_bit_size = 0;
            amount++;
        }
        if(remain_bit_size < 0){
            remain_bit_size += factor_bit_num;
            amount--;
        }
        // gmp_printf("tmp = %Zd, p = %Zd\n", tmp, p);
    }while( (not mpz_probab_prime_p(p, 100) ) or mpz_cmp_ui(tmp,2));
    mpz_clear(tmp);
}

void test_cost_powm(int bit_N, int bit_a, int bit_n){
    mpz_t a, n, N;
    mpz_init(N);
    mpz_init(a);

    gen_random_int(N, bit_N);
    gen_random_int(a, bit_a);
    gen_random_int(n, bit_n);

    clock_t start = clock();
    mpz_powm(a,a,n,N);
    
    cout<<"bit_N = "<< bit_N <<", bit_a = "<<bit_a << ", bit_n = "<< bit_n <<", powm cost = "<<(double)(clock()-start)/1000/1000<<'s'<<endl;

    // start = clock();
    // mpz_powm_sec(a,a,n,N);
    // cout<<"bit_N = "<< bit_N <<", bit_a = "<<bit_a << ", bit_n = "<< bit_n <<", powm_sec cost = "<<(double)(clock()-start)/1000/1000<<'s'<<endl;
    
}


void test(mpz_t N, int bits, ll MaxInt, int umax, bool verbose){
    
    clock_t start, end;


    bool fm;
    
    start=clock();
    fm = true;
    dynamic_scaling_pollard_pm1(N, bits, MaxInt, fm, verbose);
    end=clock();
    printf("fm = true, time= %fs. \n\n",(double)(end-start)/1000/1000);
    
    
    start=clock();
    fm = true;
    pollard_pm1_improved_IPP1_V2(N,bits,MaxInt,fm,verbose);
    end=clock();
    printf("fast_mul = true, time= %fs. \n\n",(double)(end-start)/1000/1000);
   


    start=clock();
    ll L = MaxInt, M = L+1;
    pollard_pm1_Pol74( N, bits, L, M, verbose);
    end=clock();
    printf("time= %fs. \n",(double)(end-start)/1000/1000);

    
    start=clock();
    pollard_pm1_Bis03(N, bits,verbose);
    end=clock();
    printf("time= %fs. \n",(double)(end-start)/1000/1000);
    
}



void mul_test(int max_bit_size){
    ll MaxInt = (ll)pow(2,max_bit_size),CntPrime;
    clock_t start, end;

    ll PrimeListLen=MaxInt/2+1; 
    
    ll *PrimeList =new ll[PrimeListLen];

    CalPrime(MaxInt,CntPrime,PrimeList);
    vector<Z_NR<mpz_t>> Primes;
    Primes.resize(CntPrime);
    for(int i = 0; i < CntPrime; i++){
        Primes[i] = PrimeList[i];
    }
    printf("max_bit_size = %d.", (max_bit_size));
    start = clock();
    normal_mul(Primes);
    end=clock();
    printf("normal_mul: time= %fs. \n",(double)(end-start)/1000/1000);

    start = clock();
    fast_mul(Primes);
    end=clock();
    printf("fast_mul: time= %fs. \n",(double)(end-start)/1000/1000);
}


//Compare the cost for DSP, IPP1v2, pol74 and bis03.
void algorithm_comparison(){
    mpz_t p,q,N;
    bool verbose = false;
    int bits,  factor_bit_num = 30, umax = 32; //max_bit_size = 30,
    ll MaxInt = (ll)pow(2,30)-1;
    mpz_init(N);
    mpz_init(p);
    mpz_init(q);
    for(int j = 10; j <= 11; j++){
        cout<<"/////////////////////Pollard's P-1 Comparison//////////////////////"<<endl;
        cout<<"Start test: bit(N) = "<<pow(2,j)<<", bit(P) = bit(Q) = "<<pow(2,j-1)<<", MaxInt = "<<MaxInt<<endl;
        bits = pow(2,j);
        for(int i = 0; i < 10; i++){
            cout<<"-------------------------"<<endl;
            cout<<"--(index="<<i+1<<")."<<endl;
            // gen_prime(p,bits/2);
            gen_prime_with_small_prime_factor(p, bits/2, factor_bit_num);

            gen_prime(q,bits/2);

            //gen the randomized prime p such that (p-1) = prod of 32 bit primes.

            mpz_mul(N,p,q);

            gmp_printf ("%s = %Zd\n", "p", p);
            gmp_printf ("%s = %Zd\n", "q", q);
            gmp_printf ("%s = %Zd\n", "N", N);
            
            test(N, bits, MaxInt, umax, verbose);
            cout<<"--------------------------------------"<<endl;
        }

        cout<<"----------------end-------------------"<<endl;
    }
    mpz_clear(N);
    mpz_clear(p);
    mpz_clear(q);
}


void test_expo_growth(){
    mpz_t p,q,N;
    bool verbose = false;
    int bits = 1024,  factor_bit_num = 30, umax = 0; //max_bit_size = 30,
    ll MaxInt = (ll)pow(2,30)-1, pt = 1073741789;
    mpz_init(N);
    mpz_init(p);
    mpz_init(q);

    for(int u = 0; u <= min(ceil(bits/2*log(2)/log(3)) - 31, 100); u+=20){
        cout<<"///////////////////-index-growth-test-//////////////////////"<<endl;
        if(pt > pow(3,u)){
            umax = max(1,u);
        }
        else{
            umax = u;
        }
        cout<<"Start test: bit(N) = "<<bits<<", bit(P) = bit(Q) = "<<bits/2<<", MaxInt = "<<MaxInt<<", umax = "<<umax<<endl;

        gen_prime_for_high_index(p, bits/2, factor_bit_num, u, pt);
        gen_prime(q,bits/2);

        //gen the randomized prime p such that (p-1) = prod of 32 bit primes.
        mpz_mul(N,p,q);
        gmp_printf ("%s = %Zd\n", "p", p);
        gmp_printf ("%s = %Zd\n", "q", q);
        gmp_printf ("%s = %Zd\n", "N", N);
        test(N, bits, MaxInt, umax, verbose);
        cout<<"--------------------------------------"<<endl;

        cout<<"----------------end-------------------"<<endl;
    }
    mpz_clear(N);
    mpz_clear(p);
    mpz_clear(q);
}


int main(int argc, char **argv)
{   
    //Compare the cost for DSP, IPP1v2, pol74 and bis03.
    if(atoi(argv[1])==1)
        algorithm_comparison();
    
    //test DSP and IPP1v2 with the growth of exponent.
    if(atoi(argv[1])==2)
        test_expo_growth();
    
    //Test practical cost of fast multiplication and normal multiplicaiton.
    if(atoi(argv[1])==3){
        for(int max_bit_size = 21; max_bit_size <= 26; max_bit_size++)
            mul_test(max_bit_size);
    }
    

    return 1;

}