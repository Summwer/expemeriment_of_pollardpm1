
int main()
{
    ll  MaxInt,CntPrime;
    ll  PrimeListLen;
    int tmpval;
    CntPrime=0;
    //输入范围，即查找所有小于n的素数
    // printf("\n请输入需计算素数的范围上限，输入0则结束运行：");
    // cin>>MaxInt;
    MaxInt = 1010489929;//1110489929;
    if(MaxInt<1) return 0;

    PrimeListLen=MaxInt/2+1; //其实可以用Pi(x)计算更准确的个数，以节约运行所需内存，这里懒得写了
    //存放找到素数的动态数组
    ll *PrimeList =new ll[PrimeListLen];
    //调用getPrime函数搜索小于n的所有素数
    CalPrime(MaxInt,CntPrime,PrimeList);


    mpz_t N,C,a,n,tmp,res,q;
    clock_t start,start_for,end;
    mpz_init(N);
    mpz_init(C);
    mpz_init(a);
    mpz_init(n);
    mpz_init(tmp);
    mpz_init(res);
    mpz_init(q);
    

    mpz_init_set_str(N,"27923599681213021407034926611637566469874523487438156297829021871854267185336835555834417699638124062081862065978435358275789755449922797880942372540343149490351082513378216121657620510652337065296119029212679608293116863721614934226652937927332211856317133563520857800844724599974454798616386551458666742419594273729122176671137165710858043854085715611469987477104346783638814358091242759009090423077469924079120414624179345633165061041336969051478949273930390780624527965481748666910510642670099841373434666890288335574645730875375871875051442416833323104303435247592133718723309813270278494345959197697574411387699",10);
    mpz_init_set_str(C,"7258947159508409964909201618198961485182511031330607093079469549227678636135691682866478388493813646006864746292381173045446385694437294637795198571685259602028737851962419273515782445360555545371785964483745775531692007849878818692035286995242763797416589990718815073370648852296574083130631074948733755677112003977685049443159953005105534533139681526871158638633442181528614895492079737712019891897806731978662614356734978881539561916230344201162198387028160173765253478660937724271362155417573773453551001589223547458682592804691070465038135810589207890689897698325010489999724210465785326942473397327273738540840",10);
    


    // for(ll i=0;i<CntPrime;i++){
    //     int prime = PrimeList[i];
    //     // std::cerr<<(prime)<<' ';


    // }

    mpz_init_set_str(a,"4",10);
    // mpz_init_set_str(n,"2",10);

    start=clock();
    start_for=clock();
    int prime; //max_prime = 1010489929
    int ind;
    int ind_bound = 1023;
    int k = 1; //set a lower bound of power degree.
    // mpz_powm(a,a,n,N); 
    for (int i =0;i<CntPrime;i++)
    {
        prime = PrimeList[i];
        ind_bound = std::ceil(1023/log2(prime));
        ind = min(ind_bound, k);
        mpz_ui_pow_ui(n,prime,ind);
        mpz_powm(a,a,n,N); 
        if(i % 100000==0 or i == CntPrime - 1){
            std::cerr<<"\r "<< i << "/" <<CntPrime<<", cost = "<<(double)(clock()-start_for)/1000/1000<<'s';
            start_for = clock();
            mpz_sub_ui(tmp,a,1);
            mpz_gcd(res,tmp,N);
            if (mpz_cmp_ui(res,1)!=0)
            {
                if (mpz_cmp(res,N)!=0)
                {
                    mpz_div(q,N,res);
                    std::cerr<<std::endl;
                    gmp_printf("n= %Zd \n",n);
                    gmp_printf("p= %Zd \n",res);
                    gmp_printf("q= %Zd \n",q);

                    end=clock();
                    printf("time= %fs \n",(double)(end-start)/1000/1000);
                    return 0;
                }
            }
        }
    }
    




    //打印输出找到的质数个数，提问是否保存到本地文件
    //999999999内的素数存到CSV文件，大小是478MB
    // printf("\n需把小于等于[%lld]的全部[%lld]个素数保存到CSV文件吗？（不保存请输入0）",MaxInt,CntPrime);
    // char filename[128];
    // cin>>tmpval;
    // if(tmpval<1){
    //     delete PrimeList;
    //     return 0;
    // }
    // sprintf(filename,"./prime-in-%lld.csv",MaxInt);
    // printf("\n下面开始把所有合计[%lld]个小于等于[%lld]的质数写入文件【%s】.....",CntPrime,MaxInt,filename);
    // ofstream out(filename,ios::app);//app表示每次操作前均定位到文件末尾
    // if(out.fail()){
    //     cout<<"error\n";
    // }
    // for(ll i=0;i<CntPrime;i++)
    //     out<<PrimeList[i]<<"\n";
    // out.close();
    // printf("\n小于【%lld】的[%lld]个质数已全部写入到文件[%s]\n",MaxInt,CntPrime,filename);
    delete PrimeList;   //删除保存素数的数组，回收内存
    return 0;
}
//---------------------------------------------------------------------------------------