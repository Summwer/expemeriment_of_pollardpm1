# Guidance

We provide the experimental results, which display in form of figures or tables in the article,  in our supplement materials. Our implementations about this article have been uploaded to the website [https://github.com/Summwer/expemeriment_of_Pollardpm1](https://github.com/Summwer/expemeriment_of_Pollardpm1).



## Figure1&Figure6

Figure1 and Figure6 are same. In order to help readers to understand the efficiency improvement of our algorithm more intuitively in the introduction, we have placed Figure 6 in the introduction as Figure1. The theoretical cost of our variant and IPP1 can be computed by Theorem 4.2 and Theorem 4.1. (We've implemented the fuction in the file `draw_figure_for_growth_all.py` in the folder `Figure1&Figure6`.) 



We implement their cost estimator in the file  `complexity_analysis.cpp`, one can obtain the estimated cost of our variant and IPP1 by the command:

```bash
g++ -O3 -funroll-loops -std=c++14 complexity_test.cpp complexity_analysis.cpp  utils.cpp -lgmp -lfplll -pthread -o complexity_test
./complexity_test 1
```

We've stored the result in the file `pollardpm1_cost_in_estimator.log`. 



The practical cost of our variant and IPP1 in Figure1 is obtained by running the command:

```bash
g++ -O3 -funroll-loops -std=c++14 algorithm_test.cpp generate_primes.cpp pollard_pm1.cpp utils.cpp -lgmp -lfplll -pthread -o algorithm_test
./algorithm_test 1
```

We've stored the result in the file `expo_growth_test.log`. 





## Figure 2

The practical cost of fast multiplication and normal mulitiplication shown in Figure 2 can be obtained by running the command:

```bash
g++ -O3 -funroll-loops -std=c++14 algorithm_test.cpp generate_primes.cpp pollard_pm1.cpp utils.cpp -lgmp -lfplll -pthread -o algorithm_test
./algorithm_test 3
```

We've stored the result in the log file `practical_cost_of_fm_and_nm.log`. 

To obtain the estimated cost of Fast Multiplication through its estimator, one can run the following command:

```bash
g++ -O3 -funroll-loops -std=c++14 complexity_test.cpp complexity_analysis.cpp  utils.cpp -lgmp -lfplll -pthread -o complexity_test
./complexity_test 2
```

We've stored the result in the file `fm_cost_in_estimator.log`

The theoretical cost analysis shown in Corollary 3.2 and Corollary 3.4 are implemented in the file `fm_nm_cost_comarison.py`. One can also obatined the Figure 2 by running the command:

```bash
python fm_nm_cost_comarison.py
```



## Figure5

Figure 5 gives all the practical cost of the variants of Polalrd's P-1 algorithm for factoring a randomly generated $P-1$ with small prime factors. We 've stored them in the log file `pollardpm1_test.log`. One can also get the test result after running the following command:

```bash
g++ -O3 -funroll-loops -std=c++14 algorithm_test.cpp generate_primes.cpp pollard_pm1.cpp utils.cpp -lgmp -lfplll -pthread -o algorithm_test
./algorithm_test 1 
```

Since the value of  $N = PQ$ is generated randomly, the output of the  result maybe quit different from our `pollardpm1_test.log`, but the conclusion is same, i.e. our variant is most efficient than others.











