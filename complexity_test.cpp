#include "pollard_pm1.h"
#include "complexity_analysis.h"


int main()
{
    // ll  MaxInt,CntPrime;
    // ll  PrimeListLen;
    // int tmpval;
    // CntPrime=0;
    // //输入范围，即查找所有小于n的素数
    // // printf("\n请输入需计算素数的范围上限，输入0则结束运行：");
    // // cin>>MaxInt;
    // MaxInt = 30000;//1010489929;//1110489929;
    // if(MaxInt<1) return 0;

    // PrimeListLen=MaxInt/2+1; //其实可以用Pi(x)计算更准确的个数，以节约运行所需内存，这里懒得写了
    // //存放找到素数的动态数组
    // ll *PrimeList =new ll[PrimeListLen];
    // //调用getPrime函数搜索小于n的所有素数
    // CalPrime(MaxInt,CntPrime,PrimeList);

    // int bit_a = 2048; //a的比特数最大应该是和N的比特数一样的大的。
    // int bit_N = 2048;
    // int bit_pow=32;
    

    // complexity_test(bit_N,bit_a,CntPrime,PrimeList,true,1);//fix_bit, block =1 ,to find optimal cost1
    // printf("=============================\n");
    // complexity_test(bit_N,bit_a,CntPrime,PrimeList,false,1);//accurate_bit
    
    // return 0 ;

    const char *N_str = "27923599681213021407034926611637566469874523487438156297829021871854267185336835555834417699638124062081862065978435358275789755449922797880942372540343149490351082513378216121657620510652337065296119029212679608293116863721614934226652937927332211856317133563520857800844724599974454798616386551458666742419594273729122176671137165710858043854085715611469987477104346783638814358091242759009090423077469924079120414624179345633165061041336969051478949273930390780624527965481748666910510642670099841373434666890288335574645730875375871875051442416833323104303435247592133718723309813270278494345959197697574411387699";

    int bits = 2048;

    
    // pollard_pm1_ref26(N_str,bits,CntPrime,PrimeList);
    // pollard_pm1_ref28(N_str,bits);

    // pollard_pm1_IPP1_V1(N_str,bits,CntPrime,PrimeList);
    
    // pollard_pm1_IPP1_V2(N_str,bits,CntPrime,PrimeList);

    // pair<double,int> cost = FMP_complexiy(1,32);
    // printf("FMP time cost = %3.7f from bit size of %d to %d, get bit size with %3d \n", cost.first, 1, 32, cost.second);

    //example 1: table 1 in the article
    // int bit_N = 510, pt = pow(2, 31), umax = 110, px = 2, ux = 121;
    // complexity_test(bit_N, pt , umax, px, ux);

    int bit_N = 2048, pt = 101, umax = 20, px = 91, ux = 20;
    complexity_test(bit_N, pt , umax, px, ux);

   


    // delete PrimeList;   //删除保存素数的数组，回收内存

}