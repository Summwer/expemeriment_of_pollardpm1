#include "pollard_pm1.h"
#include "complexity_analysis.h"
#include <fplll/threadpool.h>


void gen_random_int(mpz_t &p, int bit_num){
    clock_t time = clock();
    gmp_randstate_t grt;
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt, time);
    
    mpz_init(p);
    mpz_urandomb(p, grt, bit_num);
}


void gen_prime(mpz_t &p,int bit_num){
    gen_random_int(p,bit_num);
    mpz_nextprime(p, p);
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
    
    // start=clock();
    // // MaxInt = pow(2,bit_num);
    // fm = true;
    // // dynamic_scaling_pollard_pm12(N, bits, MaxInt, fm, verbose);
    // end=clock();
    // printf("fm = true, time= %fs. \n",(double)(end-start)/1000/1000);
    // printf("==================end==================\n");
    
    start=clock();
    fm = true;
    dynamic_scaling_pollard_pm1(N, bits, MaxInt, fm, verbose);
    end=clock();
    printf("fm = true, time= %fs. \n",(double)(end-start)/1000/1000);
    

    // start=clock();
    // MaxInt = pow(2,bit_num);
    // fm = false;
    // dynamic_scaling_pollard_pm1(N, bits, umax, MaxInt, fm, verbose);
    // end=clock();
    // printf("fm = false, time= %fs. \n",(double)(end-start)/1000/1000);

    
    start=clock();
    fm = true;
    pollard_pm1_improved_IPP1_V2(N,bits,MaxInt,fm,verbose);
    end=clock();
    printf("fast_mul = true, time= %fs. \n",(double)(end-start)/1000/1000);


    // start=clock();
    // fm = false;
    // pollard_pm1_improved_IPP1_V2(N,bits,MaxInt,fm,verbose);
    // end=clock();
    // printf("fast_mul = false, time= %fs. \n",(double)(end-start)/1000/1000);



    // start=clock();
    // ll L = pow(2,bit_num -1 ), M = pow(2,bit_num);
    // pollard_pm1_Pol74( N, bits, L, M, verbose);
    // end=clock();
    // printf("time= %fs. \n",(double)(end-start)/1000/1000);

    
    // start=clock();
    // pollard_pm1_Bis03(N, bits,verbose);
    // end=clock();
    // printf("time= %fs. \n",(double)(end-start)/1000/1000);


    

    
    


    
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




int main()
{

    
    ll  MaxInt,CntPrime = 0;
    ll  PrimeListLen;
    mpz_t p,q,N;
    bool verbose = true, fm;
    int tmpval, bit_num, bits, max_bit_size = 32;
    clock_t start,end;

    
    mpz_init(N);
    

    //test powm
    // for(int i = 20; i < 40; i ++){
    //     int bit_n = pow(2,i);
    //     test_cost_powm(2048, 2048, bit_n);
    // }


    //test fast_mul
    // for(int i = 10; i < 40; i ++){
    //     // int max_bit_size = pow(2,i);
    //     mul_test(i);
    // }

    // return 1;
    // bit_num = 128;
    // bits = 2*bit_num;
    
    // gen_prime(p,bit_num); 
    // gen_prime(q,bit_num);
    // mpz_mul(N,p,q);


    // bits = 510;
    // max_bit_size = 26;
    // mpz_init_set_str(p,"37566992194796047979806151568744827439436683603956030969086239785000928018433",10); //292 bit
    // mpz_init_set_str(q,"38518541858296104433730326784561600872281175198196656219582697266475597561857",10); //292 bit
    // mpz_mul(N,p,q);

    // gmp_printf ("%s = %Zd\n", "p", p);
    // gmp_printf ("%s = %Zd\n", "q", q);
    // gmp_printf ("%s = %Zd\n", "N", N);


    // test(N, bits, max_bit_size, verbose);


    // bits = 599;
    // max_bit_size = 21;
    // mpz_init_set_str(p,"2388652345228407097432206444846081867499699894635672914157712512754453335437113605268237713409",10); //292 bit
    // mpz_init_set_str(q,"591890805088336671201434105631991650498306995542057146490103163398845227831930453491713",10); //292 bit
    // mpz_mul(N,p,q);

    // gmp_printf ("%s = %Zd\n", "p", p);
    // gmp_printf ("%s = %Zd\n", "q", q);
    // gmp_printf ("%s = %Zd\n", "N", N);


    // test(N, bits, max_bit_size, verbose);


    // bits = 53;
    // max_bit_size = 21;
    // mpz_init_set_str(p,"220152098017",10); //292 bit
    // mpz_init_set_str(q,"32615125633",10); //292 bit
    // mpz_mul(N,p,q);

    // gmp_printf ("%s = %Zd\n", "p", p);
    // gmp_printf ("%s = %Zd\n", "q", q);
    // gmp_printf ("%s = %Zd\n", "N", N);


    // test(N, bits, max_bit_size, verbose);


    // bits = 53;
    // ll pt = pow(2,30);
    // int umax = 20;
    // mpz_init_set_str(p,"25165993",10); //292 bit
    // mpz_init_set_str(q,"201327937",10); //292 bit
    // mpz_mul(N,p,q);

    // gmp_printf ("%s = %Zd\n", "p", p);
    // gmp_printf ("%s = %Zd\n", "q", q);
    // gmp_printf ("%s = %Zd\n", "N", N);


    // test(N, bits, max_bit_size, verbose);


    // bits = 510;
    // max_bit_size = 31;
    // ll pt = pow(2,30);
    // int umax = 20;
    // mpz_init_set_str(p,"73582318454575190699844337236574154960295547044435008995068140893793480595457",10); 
    // mpz_init_set_str(q,"38157961040901410371145312858719796183074013434013510604185954514707631046657",10); 
    // mpz_mul(N,p,q);

    // gmp_printf ("%s = %Zd\n", "p", p);
    // gmp_printf ("%s = %Zd\n", "q", q);
    // gmp_printf ("%s = %Zd\n", "N", N);


    // test(N, bits, max_bit_size, verbose);


    // bits = 2048;
    // max_bit_size = 30;
    // mpz_init_set_str(p,"173392280438546792683818136702600967448366743521806189327229178974568108822276784275748977768555394587467983783851310057943432666398430604814797445730219422150371757816773423848241811463624899203965492680109656320607902960451073955585153506896911302285757276495392556781174371601361441446529516356977237812613",10); 
    // mpz_init_set_str(q,"161042923079321432339269541488029212509662145416259730874164975663181566574089603989613633280826778615974303933101021932326645382943220593217219685574539852154132851900274310647897485529423805846103274841564609818123788384268009542517277530372498868209862982385157094124017311593031129384417153990751676475223",10); 
    // mpz_mul(N,p,q);

    // gmp_printf ("%s = %Zd\n", "p", p);
    // gmp_printf ("%s = %Zd\n", "q", q);
    // gmp_printf ("%s = %Zd\n", "N", N);


    // test(N, bits, max_bit_size, verbose);


    // bits = 2048;
    // // max_bit_size = 30;
    // ll pt  = pow(2,30);
    // int umax = 20;
    // mpz_init_set_str(p,"96019181396760255277219321255636470660379911633960919541408367924928819826729261587560034147788097999593233455913001087902487241254958922533745323178913281215089348502438371325934262305229015090485494552543043924445989448033506694713534997843124398291233619495690644844922731879409019117429936595417224316683",10); 
    // mpz_init_set_str(q,"161042923079321432339269541488029212509662145416259730874164975663181566574089603989613633280826778615974303933101021932326645382943220593217219685574539852154132851900274310647897485529423805846103274841564609818123788384268009542517277530372498868209862982385157094124017311593031129384417153990751676475223",10); 
    // mpz_mul(N,p,q);

    bits = 2048;
    ll pt  = 101;
    int umax = 20;
    mpz_init_set_str(p,"127536485251852200948792719358413977835689001990821945846331371069431856803663977077665891741031219640226701902675912278065655006334587759006317921108084180394985344070444038427275759797644767907385215525988295874014662027877979214582360355646885031050813144428751355509592929040502133722092036494049280000001",10); 
    mpz_init_set_str(q,"140990768895481085070549250007650337537512252816130254048645374656405043961839708059519401072749709654533451197982442973527184730413451804208662666818295713195344662303677930126334816927583615437870606894218326819246458188923151403479581694033586182897276133040799720605296512753672360533098139951348684665337",10); 
    mpz_mul(N,p,q);

    

    gmp_printf ("%s = %Zd\n", "p", p);
    gmp_printf ("%s = %Zd\n", "q", q);
    gmp_printf ("%s = %Zd\n", "N", N);

    // for(int i = 0; i < 10; i++)
    test(N, bits, pt, umax, verbose);
    // mpz_clear(N);
    return 1;

}