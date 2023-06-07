#include "generate_primes.h"

//CalPrime函数计算列出小于等于MaxInt的所有素数，并把找到的素数保存到PrimeList数组
void CalPrime(ll MaxInt,ll &CntPrime,ll PrimeList[])
{
    //MaxInt 要查找素数的范围上限;
    //CntPrime 小于等于MaxInt的素数个数（不包含1）;
    //PrimeList 找到的素数存储到这个数组列表;
    ll i,j;
    clock_t start, end;

    bool *valid=new bool[MaxInt+1];
    CntPrime=0;

    //开始计算，获取当前系统时间，以便计算结束时统计耗时
    start = clock();
    for(i=2; i<=MaxInt; i++){  //将要查找素数范围内的valid全部赋值为true;
        valid[i]=true;
    }
    for(i=2; i<=MaxInt; i++)
    {
        if(valid[i])
        {
            if(MaxInt/i<i)
                break;
            for(j=i*i; j<=MaxInt; j+=i)    //将不是素数的数所对应的valid[i]赋值为false;
                valid[j]=false;
        }
    }
    
    //将[2至MaxInt]范围内的所有素数存入PrimeList数组中;
    for(i=2; i<=MaxInt; i++)
    {
        if(valid[i])
            PrimeList[CntPrime++]=i;
    }
    delete[] valid;   //释放bool数组，回收内存
    //显示计算所用的时间,我计算999999999(9亿)内的素数，发现有50847534个,计算用时[15391.218000]毫秒
    //这个算法和目前最快的算法比起来不算突出，我用同样的机器，用另一种更快但也更复杂的算法计算查找999999999(9亿)内的
    //素数耗时是80.968000毫秒
    end = clock() - start;
    printf("Generate prime list from 2 to %lld ( amount: %lld), cost = %lf s\n",MaxInt,CntPrime,(double)end/1000/1000);
}




//CalPrime函数计算列出小于等于MaxInt的所有素数，并把找到的素数保存到PrimeList数组
void CalPrime_in_parallel(ll MaxInt,ll &CntPrime,ll PrimeList[], int threads)
{
    thread_pool::thread_pool threadpool;

    threadpool.resize(threads);

    //MaxInt 要查找素数的范围上限;
    //CntPrime 小于等于MaxInt的素数个数（不包含1）;
    //PrimeList 找到的素数存储到这个数组列表;
    clock_t start, end;

    bool *valid=new bool[MaxInt+1];
    CntPrime=0;

    //开始计算，获取当前系统时间，以便计算结束时统计耗时
    start = clock();
    for(ll i=2; i<=MaxInt; i++){  //将要查找素数范围内的valid全部赋值为true;
        valid[i]=true;
    }

    ll len = MaxInt + 1;
    ll block = len/threads;
    vector<int>  t_id_begins(threads+1,0), departs(threads,0);
    t_id_begins[0] = 2;
    for(int t_id = 0; t_id < threads; t_id++){
        if(len%(threads)!=0 && ((t_id/(len%threads)) ? 0 : 1)){
            departs[t_id] = block + 1;
        }
        else{
            departs[t_id] = block;
            }
        if(t_id > 0)
            t_id_begins[t_id] = t_id_begins[t_id-1]+departs[t_id-1];
    } 
    t_id_begins[threads] = len;

    for(int t_id = 0; t_id <threads; t_id ++){
        threadpool.push([&valid, t_id, t_id_begins, MaxInt](){
            for(ll i = t_id_begins[t_id]; i < t_id_begins[t_id+1]; i++){
                if(valid[i]){
                    if(MaxInt/i<i)
                        break;
                    for(ll j=i*i; j<=MaxInt; j+=i)    //将不是素数的数所对应的valid[i]赋值为false;
                        valid[j]=false;
                }
            }
        });
    }
    threadpool.wait_work(); 
    // for(i=2; i<=MaxInt; i++)
    // {
    //     if(valid[i])
    //     {
    //         if(MaxInt/i<i)
    //             break;
            
    //         for(j=i*i; j<=MaxInt; j+=i)    //将不是素数的数所对应的valid[i]赋值为false;
    //             valid[j]=false;
    //     }
    // }
    
    //将[2至MaxInt]范围内的所有素数存入PrimeList数组中;
    for(ll i=2; i<=MaxInt; i++)
    {
        if(valid[i])
            PrimeList[CntPrime++]=i;
    }
    delete[] valid;   //释放bool数组，回收内存
    //显示计算所用的时间,我计算999999999(9亿)内的素数，发现有50847534个,计算用时[15391.218000]毫秒
    //这个算法和目前最快的算法比起来不算突出，我用同样的机器，用另一种更快但也更复杂的算法计算查找999999999(9亿)内的
    //素数耗时是80.968000毫秒
    end = clock() - start;
    printf("%lld(%lld亿)内的素数个数为%lld,计算用时[%lf]秒\n",MaxInt,MaxInt/100000000,CntPrime,(double)end/1000/1000);
}







