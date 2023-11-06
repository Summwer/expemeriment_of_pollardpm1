#include "complexity_analysis.h"
#include <algorithm>

//返回普通幂算法的时间复杂度, 
//输入：指数e, 底数A的比特数
double normal_pow_complexity(int bit_a,int bit_pow){
    return (double)2.*bit_pow*bit_a*log2(bit_a);
}

//fast multiplication complexity
// pair<double,double> FM_complexity(int bit){
//     return (double)bit*log2(bit);
// }


void gen_indices_and_bitsizes(vector<long long> &indices, vector<long long> &bitsizes, map<long long,long long> bitsizes_num){
    int index = 0;
    indices.resize(0);
    bitsizes.resize(0);
    indices.insert(indices.end(),0);

    for(auto it : bitsizes_num){
        bitsizes.insert(bitsizes.end(),it.first);
    }
    sort(bitsizes.begin(), bitsizes.end());

    for(auto it : bitsizes){
        index += bitsizes_num[it];
        indices.insert(indices.end(), index);
    }
}


double prob_prime(int i){
    return max(0.,min(1./((double)i*log(2.)),1.));
}

double prob1(int i){
    return 1.;
}

pair<double,long long> FM_complexity(int n1, int n2, double (*pr)(int), bool verbose, double C1){
    /*_summary_

    Args:
        n1 (int): lower bound bit size of list L
        n2 (int): upper bound bit size of list L
        prob (float): probability in number occurance, for all consecutive numbers, p = 1; for prime numbers, prob is a function prob = 1/(nln2)

    Returns:
        time cost, bit size of result
    */
    
    double T = 0.;
    long long amount = 0;
    map<long long, long long >  bitsizes_num;
    vector<long long> indices,bitsizes;
    indices.resize(0);
    indices.insert(indices.end(),0);
    bitsizes.resize(0);
    if(n1 < 1 or n1 != round(n1) or n2 != round(n2) or n1>n2)
        throw "n1 is not a valid value.";
    if(n1 == 1)
        n1 = 2;
    //bitsizes_num[0] is the amount of 1, since 1*a = a, so let the bitsize of 1 become 0 to avoid the increasement of bit size after product
    for(int i = n1; i < n2+1; i++){
        // bitsizes_num[i] = pow(2,(i-1));
        bitsizes_num[i] = round((*pr)(i) * pow(2,(i-1)));
        amount += bitsizes_num[i];
    }
    // print_map(bitsizes_num);
    while(amount > 1){
        if(amount % 2 != 0){
            if(bitsizes_num.find(0) == bitsizes_num.end())
                bitsizes_num[0] = 1;
            else
                bitsizes_num[0] ++;
            amount ++;
        }
        map<long long, long long> tmp_bitsizes_num;
        gen_indices_and_bitsizes(indices, bitsizes, bitsizes_num);
        // print_vector(bitsizes, 0, int(bitsizes.size()));
        for(int i = int(floor(amount/2)); i < amount; i++ ){
            //return the bitsize of the i-th number
            long long l1 = bitsizes[upper_bound(indices.begin(), indices.end(), amount-i-1) - indices.begin() - 1];
            long long l2 = bitsizes[upper_bound(indices.begin(), indices.end(), i) - indices.begin()-1]; //find the bit size of the position i, indices[j]<= i < indices[j+1]

            if(l2 > 0)
                T+=  l2 * log2(l2);
            if(tmp_bitsizes_num.find(l2+l1) == tmp_bitsizes_num.end())
                tmp_bitsizes_num[l2+l1] = 1;
            else
                tmp_bitsizes_num[l2+l1] += 1;
        }
        bitsizes_num = tmp_bitsizes_num;
        amount = int(floor(amount/2));
    }

    long long bitsize = 0; 
    for (map<long long, long long >::iterator it = bitsizes_num.begin();
		it != bitsizes_num.end(); it++) {
        bitsize = it->first;
    }
    if(verbose)
        printf("FMP time cost = %3.7f from bit size of %d to %d, get bit size with %3lld \n", T, n1, n2, bitsize);
    return make_pair(C1*T, bitsize);
}


//生成B以内全体素数集的时间复杂度估计
double gen_prime_list_complexity(int B){
    return 2*sqrt((double)B)/(log(B)-3);
}


//返回slide window模幂算法的时间复杂度
//选择最优的窗口
//输入：指数e的比特数, 底数A的比特数
double slide_window_MP_cost( int bit_a, long long bit_pow, double C1 ){

    long long L = bit_pow;

    double mincost = -1;
    double cost;
    long long mink;
    for(long long k =1; k<=L;k++){
        cost = pow(2,k-1)-1+(L-1)/(double)(k+1);
        if(mincost == -1 or mincost > cost){
            mincost = cost;
            mink = k;
        }
    }
    // cout<<"mincost = "<<mincost<<endl;
    mincost = mincost * C1 * bit_a * log2(bit_a);

    return mincost;
}


//cost of gcd
double GCD_cost(long long n, double C2){
    return (double) C2*pow(n,2);
}


//cost of a multiplication 
double M_cost(long long n){
    return (double) n*log2((double) n);
}

double IPP1v2_complexity(int bit_N, long long pt, int umax){
    int bit_pt = floor(log2(pt))+1;
    pair<double,long long> FMP_cost = FM_complexity(1, bit_pt,prob_prime);
    // cout<<umax * GCD_cost(bit_N) <<","<<umax * slide_window_MP_cost(bit_N,FMP_cost.second)<<endl;
    // cout<<endl;
    // cout<<"FMP_cost = "<<FMP_cost.first<<", GCD_cost = "<<GCD_cost(bit_N)<<", MP_cost = "<<slide_window_MP_cost(bit_N,FMP_cost.second)<<", pt = "<< FMP_cost.second<<endl;
    // cout<<"============="<<endl;

    return  FMP_cost.first + (double) umax * GCD_cost(bit_N) + (double) umax * slide_window_MP_cost(bit_N,FMP_cost.second);
}


double dynamic_scaling_pollard_pm1_complexity(int bit_N, long long pt,long long px, int ux){
    double T = 0;
    int k0 = ceil(log2(ux*log2(px))), bit_pt = floor(log2(pt))+1;
    long long bit_Pj1 = 1;
    pair<double,long long> Pj2_cost = make_pair(0,0);
    for(int j=1; j < k0+1; j++){
        if(pt >= pow(2,pow(2,j)))
            Pj2_cost = FM_complexity(pow(2,(j-1)), pow(2,j), prob_prime);
        else if(pt > pow(2,pow(2,(j-1))) and pt < pow(2,pow(2,j)))
            Pj2_cost = FM_complexity(pow(2,(j-1)), bit_pt, prob_prime);
        else
            Pj2_cost = make_pair(0,0);
        // cout<<T<<endl;
        // printf("Pj2_cost.first = %3.7f, Pj2_cost.second = %lld \n", Pj2_cost.first, Pj2_cost.second);
        T += Pj2_cost.first;
        long long bit_Pj = bit_Pj1 + Pj2_cost.second;
        bit_Pj1 += bit_Pj;
        T += M_cost(bit_Pj1) + M_cost(bit_Pj) + slide_window_MP_cost(bit_N, bit_Pj) + GCD_cost(bit_N);
        // printf("M_cost_Pj1 = %3.7f, M_cost_Pj = %3.7f, \n", M_cost(bit_Pj1),M_cost(bit_Pj));
    }
    return T;
}



pair<double,double> complexity_test(int bit_N, long long pt , int umax, long long px, int ux){
   double T1 = IPP1v2_complexity(bit_N, pt, umax);
   double T2 = dynamic_scaling_pollard_pm1_complexity(bit_N, pt, px, ux);
   printf("bit(N) = %d, pt = %lld, umax = %d, px = %lld, ux = %d\n", bit_N, pt, umax, px, ux);
   printf("Cost for IPP1v2: %e operations\n", T1);
   printf("Cost for dynamic scaling pollard P-1: %e operations\n", T2);
   printf("Ratio = %3.7f\n", T1/T2);
   return make_pair(T1,T2);
}


