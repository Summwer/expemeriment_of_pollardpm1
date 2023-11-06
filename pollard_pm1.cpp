#include "pollard_pm1.h"
#include "fplll/nr/nr_Z.inl"
#include "gmp.h"




void mpz_inits(mpz_t a, mpz_t P, mpz_t res, mpz_t tmp, mpz_t q){
    mpz_init(a);
    mpz_init(P);
    mpz_init(tmp);
    mpz_init(res);
    mpz_init(q);
}

// void mpz_clears(mpz_t a, mpz_t P, mpz_t res, mpz_t tmp, mpz_t q){

//     mpz_clear(a);
//     mpz_clear(P);
//     mpz_clear(tmp);
//     mpz_clear(res);
//     mpz_clear(q);

// }


void print_res(mpz_t N, mpz_t res){
    mpz_t q;
    mpz_init(q);
    mpz_div(q,N,res);
    std::cerr<<std::endl;
    cerr<<"Factorize successfully!"<<endl;
    gmp_printf("N= %Zd \n",N);
    gmp_printf("p= %Zd \n",res);
    gmp_printf("q= %Zd \n",q);
}




Z_NR<mpz_t> normal_mul(vector<Z_NR<mpz_t>> List){
    int Cnt = List.size();

    Z_NR<mpz_t> result;
    result = 1;
    for(int i = 0; i< Cnt; i++){
        result.mul(result,List[i]);
    }
    
    return result;
}


Z_NR<mpz_t> fast_mul(vector<Z_NR<mpz_t>> List){
    int Cnt = List.size();
    
    if(Cnt == 0){
        Z_NR<mpz_t> result;
        result = 1;
        return result;
    }
    while(Cnt > 1){
        if(Cnt % 2 == 1){
            Cnt++;
            List.resize(Cnt);
            List[Cnt-1]=1;
        }
        int mul_size = Cnt/2;
        for(int i = 0; i< mul_size; i++){
            List[i].mul(List[i],List[Cnt-i-1]);
            // cerr<<i<<"/"<<mul_size<<":"<<List[i]<<"*"<< List[Cnt-i-1]<<"="<<result[i]<<endl;
        }
        List.resize(mul_size);
        Cnt = mul_size;
    }


    return List[0];

}


Z_NR<mpz_t> fast_mul2(vector<Z_NR<mpz_t>> List){

    vector<vector<Z_NR<mpz_t>>> classified_list;
    int bitsize = 2;
    classified_list.resize(1);
    classified_list[0].resize(0);
    for(int i = 0; i < int(List.size());i++){
        while(List[i] >= pow(2,bitsize)){
            bitsize++;
            classified_list.resize(bitsize-1);
            classified_list[bitsize-2].resize(0);
        }
        classified_list[bitsize-2].insert(classified_list[bitsize-2].end(),List[i]);
    }

    List.resize(bitsize);
    for(int i = 0; i < bitsize-1; i++){
        // cout<<i+2<<","<<classified_list[i].size()<<endl;
        List[i] = fast_mul(classified_list[i]);
    }
    return fast_mul(List);
}



void pollard_pm1_Pol74(mpz_t N, int bits, ll L, ll M, bool verbose){
    printf("=====================================\n");
    printf("Test pollard p-1 method in [Pol74]...\n");
    printf("L = %lld, M = %lld \n", L, M);

    
    //Genereate prime list
    ll LCntPrime, MCntPrime; 
    ll PrimeListLen=M/2+1; 
    ll *PrimeList =new ll[PrimeListLen];
    CalPrime(M,MCntPrime,PrimeList);

    for (int i = 0;i< MCntPrime;i++) {
        ll prime = PrimeList[i];
        if(prime > L){
            LCntPrime = i;
            break;
        }
    }

    // vector<Z_NR<mpz_t>> Primes;
    // Primes.resize(MCntPrime);
    // for(int i = 0; i < LCntPrime; i++){
    //     Primes[i] = PrimeList[i];
    // }

    mpz_t a,b,P,res;//
    
    mpz_init(a);
    mpz_init(b);
    mpz_init(P);
    mpz_init(res);

    mpz_set(res,N);
    
    gmp_randstate_t grt;
    gmp_randinit_default(grt);

    int prime; 
    int ind;

    do{
        gmp_randseed_ui(grt, clock());
        mpz_urandomb(a, grt, bits); //Generate a random int from 0 to 2^x-1
        mpz_set(b,a);
        for (int i = 0;i< LCntPrime;i++) {
            prime = PrimeList[i];
            ind= bits/2/log2(prime);//e = logB/logp
            mpz_ui_pow_ui(P,prime,ind); //P=p^e
            mpz_powm(b,b,P,N); //a=a^P(mod N)
            mpz_sub_ui(res,b,1);
            mpz_gcd(res,res,N);

            if(verbose){
                printf("\r ind = %3d, %10lld/%10lld", ind, PrimeList[i], PrimeList[LCntPrime-1]);
            }
            if(mpz_cmp_ui(res,1)!=0){
                print_res(N, res);
                cout<<endl;
                delete [] PrimeList;
                mpz_clear(a);
                mpz_clear(b);
                mpz_clear(P);
                mpz_clear(res);
                return;
            }
        }
        if( mpz_cmp(res,N)==0)
            LCntPrime--;
    }
    while( mpz_cmp(res,N)==0 and LCntPrime > 0);

    if(verbose)
        cerr<<endl;
    
    
    if(mpz_cmp_ui(res,1) == 0){
        mpz_set(a,b);
        if(MCntPrime > LCntPrime){
            ll p0 = PrimeList[LCntPrime-1], tmpDelta, Delta = 0;

            for(int i = LCntPrime ; i < MCntPrime; i++){
                ll prime = PrimeList[i];
                tmpDelta = (prime - p0)/2;
                if(tmpDelta > Delta){
                    Delta = tmpDelta;
                }
                p0 = prime;
            }

            cerr<<"Delta = "<<Delta<<endl;
   
        
            // Precompute a^2 ..., a^2Delta(mod N);
            mpz_t *aList =new mpz_t[Delta];
            for(int i = 0; i < Delta; i++){
                mpz_init(aList[i]);
                mpz_powm_ui(aList[i],a,2*(i+1),N);
            }

            p0 = PrimeList[LCntPrime-1];
            mpz_powm_ui(b,b,p0,N);

            for(ll i = LCntPrime; i < MCntPrime; i++){
                if(verbose){
                    cerr<<"\r"<< i+1 << "/" << MCntPrime << ": "<< PrimeList[i];
                }
                prime = PrimeList[i];
                tmpDelta = (prime - p0)/2;
                mpz_mul(b,b,aList[tmpDelta-1]);
                // gmp_printf("\n%Zd\n",aList[tmpDelta-1]);
                mpz_mod(b,b,N); //b=b*p^(p-p0)(mod N)
                mpz_sub_ui(res,b,1);
                mpz_gcd(res,res,N);
                
                // if(prime == 1048583)
                //     throw "";
                if(mpz_cmp_ui(res,1)!=0){
                    print_res(N, res);
                    delete [] aList;
                    delete [] PrimeList;
                    mpz_clear(a);
                    mpz_clear(b);
                    mpz_clear(P);
                    mpz_clear(res);
                    return;
                }
                p0 = prime;
            }
            if(verbose)
                cerr<<endl;
        }
    }
    else{
        print_res(N, res);
        delete [] PrimeList;
        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(P);
        mpz_clear(res);
        return;
    }
    delete [] PrimeList;
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(P);
    mpz_clear(res);
    printf("\nFail to factor N.\n");
}




// void pollard_pm1_ref26(mpz_t N, int bits, ll &CntPrime,ll PrimeList[]){
//     printf("=====================================\n");
//     printf("Test pollard p-1 method in [26]...\n");
//     mpz_t a,P,tmp,res,q;
    
//     mpz_inits(a, P, res, tmp, q);
    
    
//     gmp_randstate_t grt;
//     gmp_randinit_default(grt);

//     int prime; 
//     int ind;
//     do{
//         gmp_randseed_ui(grt, clock());
//         mpz_urandomb(a, grt, bits); //Generate a random int from 0 to 2^x-1

//         for (int i =0;i<CntPrime;i++) {
//             prime = PrimeList[i];
//             ind= std::ceil(log2(PrimeList[CntPrime-1])/log2(prime));//e = logB/logp
//             mpz_ui_pow_ui(P,prime,ind); //P=p^e
//             mpz_powm(a,a,P,N); //a=a^P(mod N)
//             mpz_sub_ui(tmp,a,1);
//             mpz_gcd(res,tmp,N);
//         }
//     }
//     while( (mpz_cmp_ui(res,1)==0 or mpz_cmp(res,N)==0 ));
    
//     if((mpz_cmp_ui(res,1)!=0 and mpz_cmp(res,N)!=0 )){
//         print_res(N,res);
//         return;
//     }
    
//     printf("\nFail to factor N through method [26].\n");
    
// }



void pollard_pm1_Bis03(mpz_t N, int bits,bool verbose){
    printf("=====================================\n");
    printf("Test pollard p-1 method in [Bis03]\n");
    mpz_t a,P,res;

    mpz_init(a);
    mpz_init(P);
    mpz_init(res);
    
    gmp_randstate_t grt;
    gmp_randinit_default(grt);
    

    ll i = 1;
    
    gmp_randseed_ui(grt, clock());
    mpz_urandomb(a, grt, bits); //Generate a random int from 0 to 2^x-1

    do{
        mpz_set_ui(P,i);
        mpz_powm(a,a,P,N); 
        mpz_sub_ui(res,a,1);
        mpz_gcd(res,res,N);
        i++;
        if(verbose)
            cerr<<"\r"<<i+1;
    }
    while( (mpz_cmp_ui(res,1)==0));
    if(verbose)
        cerr<<endl;
    
    if((mpz_cmp_ui(res,1)!=0 and mpz_cmp(res,N)!=0 )){
        print_res(N,res);
        mpz_clear(a);
        mpz_clear(P);
        mpz_clear(res);
        return;
    }
    if(mpz_cmp(res,N) == 0)
        printf("\ngcd(d,N) = N!!");
    printf("\nFail to factor N.\n");
    mpz_clear(a);
    mpz_clear(P);
    mpz_clear(res);
}

void pollard_pm1_IPP1_V1(mpz_t N, int bits,ll &CntPrime,ll PrimeList[]){
    printf("=====================================\n");
    printf("Test pollard p-1 method IPP1_V1...\n");
    mpz_t a,P,tmp,res,q;
    mpz_inits(a, P, res, tmp, q);
    int prime;
    int ind;

    clock_t start,start_for,end;
    start=clock();
    start_for=clock();

    mpz_set_ui(P,PrimeList[0]);

    for(int i =1;i<CntPrime;i++) {
        prime = PrimeList[i];
        mpz_mul_ui(P,P,prime); //P=P*p
        if(i % 100000==0 or i == CntPrime - 1){
            std::cerr<<"\r "<< i << "/" <<CntPrime<<", cost = "<<(double)(clock()-start_for)/1000/1000<<'s';
            start_for = clock();
        }
    }
    // gmp_printf("P= %Zd \n",P);

    gmp_randstate_t grt;
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt, clock());
    mpz_urandomb(a, grt, bits); //Generate a random int from 0 to 2^x-1


    do{
        mpz_powm(a,a,P,N); //a=a^P(mod N)
        mpz_sub_ui(tmp,a,1);
        mpz_gcd(res,tmp,N);
    }
    while(mpz_cmp_ui(res,1)==0);

    if((mpz_cmp_ui(res,1)!=0 and mpz_cmp(res,N)!=0 )){
        mpz_div(q,N,res);
        std::cerr<<std::endl;
        gmp_printf("N= %Zd \n",N);
        gmp_printf("p= %Zd \n",res);
        gmp_printf("q= %Zd \n",q);
    }
    else
        printf("Fail to factor N through method IPP1_V2.\n");

    end=clock();
    printf("time= %fs. \n",(double)(end-start)/1000/1000);

    mpz_clear(a);
    mpz_clear(P);
    mpz_clear(res);
}


void pollard_pm1_IPP1_V2(mpz_t N, int bits,ll &CntPrime,ll PrimeList[]){
    printf("=====================================\n");
    printf("Test pollard p-1 method IPP1_V2...\n");
    mpz_t a,P,tmp,res,q;
    mpz_inits(a, P, res, tmp, q);
    int prime;

    clock_t start,start_for,end;
    start=clock();
    start_for=clock();

    mpz_set_ui(P,PrimeList[0]);

    for(int i =1;i<CntPrime;i++) {
        prime = PrimeList[i];
        mpz_mul_ui(P,P,prime); //P=P*p
        if(i % 100000==0 or i == CntPrime - 1){
            std::cerr<<"\r "<< i << "/" <<CntPrime<<", cost = "<<(double)(clock()-start_for)/1000/1000<<'s';
            start_for = clock();
        }
    }
    // gmp_printf("P= %Zd \n",P);

    gmp_randstate_t grt;
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt, clock());
    mpz_urandomb(a, grt, bits); //Generate a random int from 0 to 2^x-1

    int ind_bound = bits/2;
    int ind = 0;
    
    do{
        mpz_powm(a,a,P,N); //a=a^P(mod N)
        mpz_sub_ui(tmp,a,1);
        mpz_gcd(res,tmp,N);
        ind++;
    }
    while(mpz_cmp_ui(res,1)==0 and ind < ind_bound - 1);

    

    if((mpz_cmp_ui(res,1)!=0 and mpz_cmp(res,N)!=0 )){
        mpz_div(q,N,res);
        std::cerr<<std::endl;
        gmp_printf("N= %Zd \n",N);
        gmp_printf("p= %Zd \n",res);
        gmp_printf("q= %Zd \n",q);
    }
    else if(mpz_cmp_ui(res,1)==0 ){
        ll px = PrimeList[CntPrime-1];
        ll max = 2*px;
        while(mpz_cmp_ui(res,1)==0){
            if(px < max)    
                px += 2;
            else
                px += 1;
            mpz_powm_ui(a,a,px,N); //a=a^px(mod N)
            mpz_sub_ui(tmp,a,1);
            mpz_gcd(res,tmp,N);
        }
    }
    else
        printf("Fail to factor N through method IPP1_V2.\n");
    
    end=clock();
    printf("time= %fs. \n",(double)(end-start)/1000/1000);




}



//IPP1v2
void pollard_pm1_improved_IPP1_V2(mpz_t N, int bits, ll MaxInt, bool fm, bool verbose){
    printf("=====================================\n");
    printf("[IPP1v2]Test pollard p-1 method improved IPP1_V2...\n");

    
    ll PrimeListLen=MaxInt/2+1, CntPrime; 
    ll *PrimeList =new ll[PrimeListLen];
    CalPrime(MaxInt,CntPrime,PrimeList);    

    vector<Z_NR<mpz_t>> Primes;
    Primes.resize(CntPrime);
    for(int i = 0; i < CntPrime; i++){
        Primes[i] = PrimeList[i];
    }


    mpz_t a,P,tmp,res,q;
    mpz_inits(a, P, res, tmp, q);
    mpz_set_ui(P,1);
    int prime;

    clock_t start,start_for,end;
    start=clock();

Narrow_prime_set:
    printf("Size of Prime List is %lld.\n", CntPrime);

    if(fm){
        Z_NR<mpz_t> tmpP = fast_mul(Primes);
        tmpP.get_mpz(P);
    }
    else{
        Z_NR<mpz_t> tmpP = normal_mul(Primes);
        tmpP.get_mpz(P);
    }
    Primes.clear();

    gmp_randstate_t grt;
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt, clock());
    mpz_urandomb(a, grt, bits); //Generate a random int from 0 to 2^x-1
    
    int ind_bound = floor((double)bits/2)+1;
    int ind = 0;
    
    do{
        mpz_powm(a,a,P,N); //a=a^P(mod N)
        mpz_sub_ui(res,a,1);
        mpz_gcd(res,res,N);
        ind++;
        if(verbose)
            cerr<<"\r"<<ind+1<<"/"<<ind_bound;
    }
    while(mpz_cmp_ui(res,1)==0 and ind < ind_bound - 1);
    if(verbose)
        cerr<<endl;

    if((mpz_cmp_ui(res,1)!=0 and mpz_cmp(res,N)!=0 )){
        print_res(N,res);
        delete [] PrimeList;
        return;
    }
    else if(mpz_cmp_ui(res,1)==0 ){
        ll px = PrimeList[CntPrime-1];
        ll max = 2*px;
        while(mpz_cmp_ui(res,1)==0){
            if(px < max)    
                px += 2;
            else
                px += 1;
            if(verbose)
                cerr<<"\r"<<px;
            mpz_powm_ui(a,a,px,N); //a=a^px(mod N)
            mpz_sub_ui(res,a,1);
            mpz_gcd(res,res,N);
            
        }
        if(verbose)
            cerr<<endl;
        if((mpz_cmp_ui(res,1)!=0 and mpz_cmp(res,N)!=0 )){
            print_res(N,res);
            delete [] PrimeList;
            return;
        }if(mpz_cmp(res,N)==0){
            cout<<"In last step of  [IPP1v2], gcd(d,N) = N unfortunately."<<endl;
        }
        
    }
    else if(mpz_cmp(res,N)==0){
        CntPrime -=1;
        goto Narrow_prime_set;
    }
    printf("\nFail to factor N through method improved IPP1_V2.\n");
    delete [] PrimeList;
}




void generate_partial_primes(vector<Z_NR<mpz_t>> &Primes, int k, vector<vector<Z_NR<mpz_t>>> &Partial_Primes){
    ll CntPrime = Primes.size(), i = 0;
    Partial_Primes.resize(k+1);
    Partial_Primes[0].resize(1);
    Partial_Primes[0][0] = 2;
    for(int j = 1; j <= k; j++){
        double D1 = pow(2,pow(2,j-1)), D2 = pow(2,pow(2,j));
        cerr<<j << "/" << k << ": D1 = " << D1 << ", D2 = " << D2 <<endl;
        Partial_Primes[j].resize(0);
        for(i = 0; i < Primes.size(); i++){
            if(Primes[i] > D2){
                break;
            }
        }
        Partial_Primes[j].insert(Partial_Primes[j].end(),Primes.begin(),Primes.begin()+i);
        Primes.erase(Primes.begin(),Primes.begin()+i);
        if(Primes.size()==0)
            break;
    }
    Primes.clear();
}


//ours 
void dynamic_scaling_pollard_pm1(mpz_t N, int bits, ll MaxInt, bool fm, bool verbose){
    printf("=====================================\n");
    printf("[ours]Test Dynamic Scaling Pollard's P-1 Algorithm...\n");

    ll CntPrime; 
    // MaxInt = min( (ll)pow(2,bits/2-1),MaxInt);

    ll PrimeListLen=MaxInt/2+1; 
    ll *PrimeList =new ll[PrimeListLen];
    CalPrime(MaxInt,CntPrime,PrimeList);

    vector<Z_NR<mpz_t>> Primes;
    Primes.resize(CntPrime);
    cout<<"CntPrime = "<<CntPrime<<endl;
    for(int i = 0; i < CntPrime; i++){
        Primes[i] = PrimeList[i];
    }
    delete [] PrimeList;

    // ll  D1, D2;
    int k = ceil(log2((double)bits/2));

    cout<<"k = "<<k<<endl;

    //Preclassify primes
    vector<vector<Z_NR<mpz_t>>> Partial_Primes;
    vector<Z_NR<mpz_t>> Pj, Pj1, Pj2;
    Pj.resize(k+1);
    Pj1.resize(k+1);
    Pj2.resize(k+1);
    generate_partial_primes(Primes, k, Partial_Primes);


    mpz_t a,b,res, n;
    bool initialized;
    Z_NR<mpz_t>  p0, Delta, tmpDelta;
 
    mpz_init_set_ui(a,2);
    mpz_init_set_ui(b,4); //b=a^2(mod N)
    mpz_init(res);
    mpz_init(n);
    
    // int prime,p0,Delta=0,ind,i0=0;
    // gmp_randstate_t grt;
    // gmp_randinit_default(grt);
    // gmp_randseed_ui(grt, clock());
    // mpz_urandomb(a, grt, bits); //Generate a random int from 0 to 2^x-1

    mpz_sub_ui(res,b,1);
    mpz_gcd(res,res,N); //d=gcd(b-1,N)
    if(mpz_cmp_ui(res,1)!=0){
        print_res(N, res);
        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(res);
        mpz_clear(n);
        Partial_Primes.clear();
        Pj.clear();
        Pj1.clear();
        Pj2.clear();
        return;
    }else{
        //set Pj, Pj1, Pj2
        Pj[0] = 2, Pj1[0] = 1, Pj2[0] = 1, Pj1[1] = 2, mpz_set(a,b);
        initialized = false;
JLOOP:  
        for(int j = 1; j <= k ; j++){
            if(not initialized){
                // D1 = pow(2,pow(2,j-1)), D2 = pow(2,pow(2,j));
                if(j < k){
                    if(fm)
                        Pj2[j] = fast_mul(Partial_Primes[j]);
                    else
                        Pj2[j] = normal_mul(Partial_Primes[j]);
                    Partial_Primes[j].clear();
                    Pj[j].mul(Pj1[j],Pj2[j]);
                    Pj1[j+1].mul(Pj1[j],Pj[j]);
                }
                else{
                    Pj[j] = Pj1[j];
                }
            }

            if(verbose)
                cerr<<"\r"<< j << "/" << k << ": log2(D1) = " << pow(2,j-1) << ", log2(D2) = " << pow(2,j);
            
            Pj[j].get_mpz(n);
            mpz_powm(b,b,n,N); 
            mpz_sub_ui(res,b,1);
            mpz_gcd(res,res,N);
        
            if(mpz_cmp_ui(res,1)!=0){
                if(verbose)
                    cout<<endl;
                while(mpz_cmp(res,N)==0){
                    if(j<1)
                        throw "Impossible situation appears, there's some error in code.";

                    //roll back to a and multi p separately
                    printf("\rgcd(d,N) = N, enter call back process, j = %d!!", j);
                
                    Pj1[j].get_mpz(n);
                    mpz_powm(b,a,n,N); 
                    mpz_sub_ui(res,b,1);
                    mpz_gcd(res,res,N);
                   
                    if(mpz_cmp_ui(res,1)==0){
                        for(ll i = 0; i < Partial_Primes[j].size(); i ++){
                            Partial_Primes[j][i].get_mpz(n);
                            mpz_powm(b,b,n,N); 
                            mpz_sub_ui(res,b,1);
                            mpz_gcd(res,res,N);
                            if(mpz_cmp_ui(res,1)!=0){
                                print_res(N, res);
                                mpz_clear(a);
                                mpz_clear(b);
                                mpz_clear(res);
                                mpz_clear(n);
                                Partial_Primes.clear();
                                Pj.clear();
                                Pj1.clear();
                                Pj2.clear();
                                return;
                            }
                        }
                    }
                    else if(mpz_cmp(res,N)!=0 and mpz_cmp_ui(res,1)!=0){
                        print_res(N, res);
                        mpz_clear(a);
                        mpz_clear(b);
                        mpz_clear(res);
                        mpz_clear(n);
                        Partial_Primes.clear();
                        Pj.clear();
                        Pj1.clear();
                        Pj2.clear();
                        return;
                    }
                    mpz_set(b,a);
                    j--;
                }
                if(verbose)
                    cout<<endl;
                if(mpz_cmp(res,N)!=0 and mpz_cmp_ui(res,1)!=0){
                    print_res(N, res);
                    mpz_clear(a);
                    mpz_clear(b);
                    mpz_clear(res);
                    mpz_clear(n);
                    Partial_Primes.clear();
                    Pj.clear();
                    Pj1.clear();
                    Pj2.clear();
                    return;
                }
            }
            else{
                mpz_set(a,b);
            }
        }
    }

    if(verbose)
        cerr<<endl;


    if(Partial_Primes[k].size() > 1){

        tmpDelta = 0;
        mpz_t tmp, tmp2; 
        mpz_init(tmp);
        mpz_init(tmp2);
        for(int i = 1 ; i < Partial_Primes[k].size(); i++){
            Z_NR<mpz_t> prime = Partial_Primes[k][i];
            tmpDelta.sub(prime, p0);
            tmpDelta.get_mpz(tmp);
            mpz_set_ui(tmp2,2);

            mpz_cdiv_q(tmp, tmp, tmp2);
            tmpDelta = tmp;
            
            if(tmpDelta > Delta){
                Delta = tmpDelta;
            }
            p0 = prime;
        }

        cerr<<"Delta = "<<Delta<<endl;

        //Precompute a^2 ..., a^2Delta(mod N);
        mpz_t *aList =new mpz_t[Delta.get_si()];
        for(int iter = 0; iter < Delta.get_si(); iter++){
            mpz_init(aList[iter]);
            mpz_powm_ui(aList[iter],a,2*(iter+1),N);
        }
    
        p0 = Partial_Primes[k][0];
        p0.get_mpz(n);
        mpz_powm(b,b,n,N);
        for(int i = 1 ; i < Partial_Primes[k].size(); i++){
            if(verbose){
                cerr<<"\r"<< i+1 << "/" << int(Partial_Primes[k].size()) << ": "<< Partial_Primes[k][i];
                // if( Partial_Primes[k][i] == 343267){
                //     cerr<<endl;
                //     // throw "";
                // }
            }
            Z_NR<mpz_t> prime = Partial_Primes[k][i];
            tmpDelta.sub(prime, p0);
            tmpDelta.get_mpz(tmp);
            mpz_set_ui(tmp2,2);
            mpz_cdiv_q(tmp, tmp, tmp2);
            tmpDelta = tmp;
        
            mpz_mul(b,b,aList[tmpDelta.get_si()-1]);
            // gmp_printf("\n%Zd\n",aList[tmpDelta-1]);
            mpz_mod(b,b,N); //b=b*p^(p-p0)(mod N)
            mpz_sub_ui(res,b,1);
            mpz_gcd(res,res,N);
            if(mpz_cmp_ui(res,1)!=0){
                if(mpz_cmp(res,N) == 0){
                    prime.get_mpz(n);
                    mpz_set_ui(a,4);
                    mpz_powm(b,a,n,N);
                    goto JLOOP;
                }
                print_res(N, res);
                mpz_clear(a);
                mpz_clear(b);
                mpz_clear(res);
                mpz_clear(n);
                Partial_Primes.clear();
                Pj.clear();
                Pj1.clear();
                Pj2.clear();
                return;
            }
            p0 = prime;
        }
        if(verbose)
            cerr<<endl;
    }
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(res);
    mpz_clear(n);
    Partial_Primes.clear();
    Pj.clear();
    Pj1.clear();
    Pj2.clear();
    cerr<< "Fail to find the solution"<<endl;
    return;
}






int find_index(ll &CntPrime,ll PrimeList[],unsigned long D){
    ll prime;
    for(int i = 0;i<CntPrime;i++) {
        prime = PrimeList[i];
        if(prime > D){
            return i-1;
        }
    }
    return CntPrime-1;
}



// void dynamic_scaling_pollard_pm1_with_block_partition(mpz_t N, int bits,ll &CntPrime,ll PrimeList[], unsigned long D, int k, int blocksize){
//     printf("=====================================\n");
//     printf("[ours]Test Dynamic Scaling Pollard's P-1 Algorithm with block partition...\n");
//     mpz_t a,P,tmp,res,q;
//     unsigned long new_D;
//     int PPsize = k;
//     mpz_t **PP = new mpz_t*[PPsize];
//     int *Psizes = new int[PPsize];
//     mpz_inits(a, P, res, tmp, q);
//     mpz_set_ui(P,1);
//     mpz_t b,Ptmp;
//     mpz_init(b);
//     mpz_init(Ptmp);
//     int prime,p0,Delta=0,ind,i=0;

//     clock_t start,start_for,end;
//     start=clock();

//     gmp_randstate_t grt;
//     gmp_randinit_default(grt);
//     gmp_randseed_ui(grt, clock());
//     mpz_urandomb(a, grt, bits); //Generate a random int from 0 to 2^x-1

//     mpz_powm_ui(b,a,2,N); //b=a^2(mod N)
//     mpz_sub_ui(tmp,b,1);
//     mpz_gcd(res,tmp,N); //d=gcd(b-1,N)
//     if(mpz_cmp_ui(res,1)!=0){
//         if(mpz_cmp(res,N)==0){
//             //return gcd(a+1,N)
//             mpz_add_ui(tmp,a,1);
//             mpz_gcd(res,tmp,N);
//             print_res(N, res);
//         }
//         else{
//             print_res(N, res);
//         }
//     }
//     else{
//         // mpz_sqrt(tmp,N);
//         // mpz_sub_ui(tmp,tmp,1);
//         // mpz_cdiv_q_ui(tmp, tmp,2);
//         // mpz_root(tmp,tmp,k); //(tmp)^(1/k)
//         // D = mpz_get_ui(tmp);
//         double lD = log2(D);
//         std::cerr<<"D = "<<D<<std::endl;

//         int D_index = find_index(CntPrime, PrimeList, D);
//         int pD_index = 0;
        
//         Psizes[0] = int(ceil((D_index-pD_index)/float(blocksize))); // |B \cap D|/blocksize
//         PP[0] = new mpz_t[Psizes[0]];
//         int s = -1;

//         //Multiply all primes in range of [0,D]
//         printf("Multiply all primes in range of [0,%lu) with blockwise partition.\n", D);
//         for(i=0;i<=D_index;i++) {
//             if(i%blocksize==0){
//                 s++;
//                 mpz_init_set_ui(PP[0][s],1);
//             }
//             prime = PrimeList[i];
//             ind= std::ceil(lD/log2(prime));//e = logD/logp
//             mpz_ui_pow_ui(tmp,prime,ind); //tmp=p^e
//             mpz_mul(PP[0][s],PP[0][s],tmp); //P=P*tmp
//             if(i % 100000==0 or i == D_index){
//                 std::cerr<<"\r "<< i+1 << "/" <<CntPrime<<", cost = "<<(double)(clock()-start_for)/1000/1000<<'s';
//                 start_for = clock();
//             }
//         }
//         for(int s = 0; s<Psizes[0]; s++){
//             mpz_set(a,b); //a = b
//             mpz_powm(b,b,PP[0][s],N); //b=b^P(mod N)
//             mpz_sub_ui(tmp,b,1);
//             mpz_gcd(res,tmp,N);
            
//             if(mpz_cmp(res,N)==0){//If primes are too many, then we should narrow the set.
//                 mpz_set(b,a); // b = a
//                 for(int iter = s*blocksize;iter<D_index;iter++) {
//                     prime = PrimeList[iter];
//                     ind= std::ceil(lD/log2(prime));
//                     for(int e=1; e <= ind; e++){
//                         mpz_powm_ui(b,b,prime,N);
//                         mpz_sub_ui(tmp,b,1);
//                         mpz_gcd(res,tmp,N);
//                         if(mpz_cmp_ui(res,1)!=0){
//                             print_res(N, res);
//                             end=clock();
//                             printf("time= %fs. \n",(double)(end-start)/1000/1000);
//                             return;
//                         }
//                     }
//                 }
//             }
//             else if(mpz_cmp_ui(res,1)!=0){
//                 print_res(N, res);
//                 end=clock();
//                 printf("time= %fs. \n",(double)(end-start)/1000/1000);
//                 return;
//             }
//             else{
//                 mpz_set(a,b);
//             }
//         }
//         printf("\n");
        
//         //Multiply all primes in range of [D^(2^(j-1)),D^(2^j)]
//         pD_index = D_index;
//         unsigned long previous_D;
//         for(int j = 1; j < k; j++){
//             previous_D = pow(D,pow(2,j-1));
//             new_D = pow(D,pow(2,j));

//             D_index = find_index(CntPrime, PrimeList, new_D);
//             Psizes[j] = int(ceil((D_index-pD_index)/float(blocksize))); // |B \cap D|/blocksize
//             PP[j] = new mpz_t[Psizes[j]];
//             s = -1;

//             printf("\nj = %d, D^(2^%d) = %lu, max(prime)= %lld \n",j,j,new_D,PrimeList[CntPrime-1]);

//             printf("Multiply all primes in range of [%lu,%lu) with blockwise partition.\n",previous_D, new_D);
//             //Ptmp = Pj*prod(..)
//             for(i = pD_index; i <= D_index; i++){
//                 if(i%blocksize==0){
//                     s++;
//                     mpz_init_set_ui(PP[j][s],1);
//                 }

//                 prime = PrimeList[i];
//                 ind= std::ceil(pow(2,j)*lD/log2(prime));//e = logD/logp
//                 mpz_ui_pow_ui(tmp,prime,ind); //tmp=p^e
//                 mpz_mul(PP[j][s],PP[j][s],tmp); //P=P*tmp
//                 if(i % 100000==0 or i == D_index){
//                     std::cerr<<"\r "<< i+1 << "/" <<CntPrime<<", cost = "<<(double)(clock()-start_for)/1000/1000<<'s';
//                     start_for = clock();
//                 }
//             }
       
//             for(int t = 0; t <= j; t++){
//                 for(int ii = 0; ii < max(1,int(pow(2,j-t-1)));ii++){
//                     for(int s = 0; s < Psizes[t]; s++){
//                         mpz_powm(b,b,PP[t][s],N); //b=b^Ptmp(mod N)
//                         mpz_sub_ui(tmp,b,1);
//                         mpz_gcd(res,tmp,N);

//                         if(mpz_cmp(res,N)==0){//If primes are too many, then we should narrow the set.
//                             mpz_set(b,a); // b = a
//                             for(int iter =pD_index; iter <D_index;iter++) {
//                                 prime = PrimeList[iter];
//                                 ind= std::ceil(pow(2,j)*lD/log2(prime));
//                                 for(int e=1; e <= ind; e++){
//                                     mpz_powm_ui(b,b,prime,N);
//                                     mpz_sub_ui(tmp,b,1);
//                                     mpz_gcd(res,tmp,N);
//                                     if(mpz_cmp_ui(res,1)!=0){
//                                         print_res(N, res);
//                                         return;
//                                     }
//                                 }
//                             }
//                         }
//                         else if(mpz_cmp_ui(res,1)!=0){
//                             print_res(N, res);
//                             end=clock();
//                             printf("time= %fs. \n",(double)(end-start)/1000/1000);
//                             return;
//                         }
//                         else{
//                             mpz_set(a,b);
//                         }
//                     }
//                 }
//             }
//             pD_index = D_index;
//         }

//         ////改
//         p0 = PrimeList[pD_index];

//         int tmpDelta = 0;
//         for(int iter=pD_index+1; iter<CntPrime;iter++){
//             prime = PrimeList[iter];
//             if(prime > pow(D,pow(2,k))){
//                 break;
//             }
//             tmpDelta = (prime - p0)/2;
//             if(tmpDelta > Delta){
//                 Delta = tmpDelta;
//             }
//             p0 = prime;
//         }

//         printf("\nDelta = %d\n", Delta);

//         //Precompute a^2 ..., a^2Delta(mod N);
//         mpz_t *aList =new mpz_t[Delta];
//         for(int iter = 0; iter < Delta; iter++){
//             mpz_init(aList[iter]);
//             mpz_powm_ui(aList[iter],a,2*(iter+1),N);
            
//         }
        
//         p0 = PrimeList[pD_index];
//         mpz_powm_ui(b,b,p0,N);
//         for(int iter=pD_index+1; iter<CntPrime;iter++){
//             prime = PrimeList[iter];
//             if(prime > pow(D,pow(2,k))){
//                 break;
//             }
//             tmpDelta = (prime - p0)/2;
//             mpz_mul(b,b,aList[tmpDelta-1]);
//             // gmp_printf("\n%Zd\n",aList[tmpDelta-1]);
//             mpz_mod(b,b,N); //b=b*p^(p-p0)(mod N)
//             mpz_sub_ui(tmp,b,1);
//             mpz_gcd(res,tmp,N);
//             if(mpz_cmp_ui(res,1)!=0){
//                 print_res(N, res);
//                 end=clock();
//                 printf("time= %fs. \n",(double)(end-start)/1000/1000);
//                 return;
//             }
//             p0 = prime;
//         }

//         if(mpz_cmp_ui(res,1)==0){
//             printf("\nFail to factor N through Dynamic Scaling Pollard's P-1 Algorithm with block partition...\n");
//         }
//     }

    
    
//     end=clock();
//     printf("time= %fs. \n",(double)(end-start)/1000/1000);

    
// }


