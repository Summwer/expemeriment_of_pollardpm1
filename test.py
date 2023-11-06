
from math import log, log2, sqrt, floor, ceil
from scipy.special import lambertw
from copy import deepcopy
from bisect import bisect

import warnings
warnings.filterwarnings('ignore')

#cost of gcd
def GCD_cost(n):
    return n**2

#cost of a multiplication 
def M_cost(n):
    return n*log2(n)

#return: cost of fast multiplication, and the bit size of the product.
def FMP_cost(n):
    # $print(n)
    # print(n,(2**(n-1)) * (n**2), (2**n)*(n-1)/(n*log(2)))
    return (2**(n-1)) * (n**2), (2**n)*(n-1)/(n*log(2))


def prime_prob(n):
    return max(0.,min(1./(n*log(2)),1.))

def prob_1(_):
    return 1

def FMP_cost_n1_n2(n1,n2,prob = prime_prob):
    """_summary_

    Args:
        n1 (int): lower bound bit size of list L
        n2 (int): upper bound bit size of list L
        prob (float): probability in number occurance, for all consecutive numbers, p = 1; for prime numbers, prob is a function prob = 1/(nln2)

    Returns:
        _type_: _description_
    """
    
    T = 0
    
    if(n1 < 1 or n1 != round(n1) or n2 != round(n2) or n1>n2):
        raise "n1 is not a valid value."
    if(n1 == 1):
        bitsizes_num = dict([(i, round(prob(i) * 2**(i-1))) for i in range(2,n2+1)])#bitsizes_num[0] is the amount of 1, since 1*a = a, so let the bitsize of 1 become 0 to avoid the increasement of bit size after product
    else:
        bitsizes_num = dict([(i, 2**(i-1)) for i in range(n1,n2+1)]) 
    amount = sum(bitsizes_num.values())
    bitsizes = sorted(bitsizes_num.keys())
    indices = [sum([bitsizes_num[i] for i in bitsizes[:j]]) for j in range(len(bitsizes_num))] + [amount]
    
    while(amount > 1):
        print(amount)
        if(amount % 2 != 0):
            if(0 not in bitsizes_num):
                bitsizes_num[0] = 1
            else:
                bitsizes_num[0] += 1
            amount += 1
            bitsizes = sorted(bitsizes_num.keys())
            indices = [sum([bitsizes_num[i] for i in bitsizes[:j]]) for j in range(len(bitsizes_num))] + [amount]
            
        tmp_bitsizes = {}
        for i in range(amount//2, amount):
            print("\r%10d/%10d" %(i+1 - amount//2,amount//2), end='')
            #return the bitsize of the i-th number
            l1 = bitsizes[bisect(indices,amount-i-1)-1]
            l2 = bitsizes[bisect(indices,i)-1]
            # print(bisect(indices,amount-i-1), bisect(indices,i))
            # print(l1,l2)
            if(l2 > 0):
                T+=  l2 * log2(l2)
            if(l2+l1 not in tmp_bitsizes):
                tmp_bitsizes[l2+l1] = 1
            else:
                tmp_bitsizes[l2+l1] +=1
        print()
        bitsizes_num = deepcopy(tmp_bitsizes)
        amount //=2
        bitsizes = sorted(bitsizes_num.keys())
        indices = [sum([bitsizes_num[i] for i in bitsizes[:j]]) for j in range(len(bitsizes_num))] + [amount]

    return T , list(bitsizes_num.keys())[0]


# input: bit size of the exponent
# return: Modular Powering cost
def MP_cost(n):
    w = 2*lambertw(sqrt((n-1)*log(2)))/log(2)-1
    return 2**(w-1) - 1 + (n-1)/(w+1)

def IPP1v2_complexity(bit_N, pt, umax):
    bit_pt = floor(log2(pt))+1
    print(bit_pt)
    T_FMP, bitsize_FMP = FMP_cost_n1_n2(1, bit_pt, prime_prob)
    return  T_FMP + umax * GCD_cost(bit_N) + umax * MP_cost(bitsize_FMP)*M_cost(bit_N)

def dynamic_scaling_pollard_pm1_complexity(bit_N, pt, px, ux):
    T = 0
    k0 = ceil(log2(ux*log2(px)))
    bit_Pj1 = 1
    bit_pt = floor(log2(pt))+1
    
    for j in range(1,k0+1):
        if(pt >= 2**(2**j)):
            Tpj2, bit_Pj2 = FMP_cost_n1_n2(2**(j-1), 2**j, prime_prob)
        elif(pt > 2**(2**(j-1)) and pt < 2**(2**j)):
            Tpj2, bit_Pj2 = FMP_cost_n1_n2(2**(j-1), bit_pt, prime_prob)
        else:
            bit_Pj2, Tpj2 = 0, 0
        T += Tpj2
        bit_Pj = bit_Pj1 + bit_Pj2
        bit_Pj1 += bit_Pj
        T += M_cost(bit_Pj1) + M_cost(bit_Pj) + MP_cost(bit_Pj)*M_cost(bit_N) + GCD_cost(bit_N)
    return T

def gen_bit_Pj2(j):
    bit_Pj2 = [1, FMP_cost(2)[1] - 1]   
    for i in range(2,j):
        # print(FMP_cost(2**j)[1])
        bit_Pj2.append(FMP_cost(2**i)[1] - FMP_cost(2**(i-1))[1])
    return bit_Pj2

def gen_bit_Pj1(j,bit_Pj2):
    Sum = 0
    Pj1 = [None]
    for i in range(j):
        Sum += bit_Pj2[i] * 2**(j-i-1)
        Pj1.append(Sum)
    return Pj1


bit_N = 510
pt = 2**31
umax = 110
px = 3
ux = 110
T1 = IPP1v2_complexity(bit_N, pt, umax)
# T2 = dynamic_scaling_pollard_pm1_complexity(bit_N, pt, px, ux)

# print(float(T1)/float(T2))
# print(FMP_cost(0,4))
# print(FMP_cost_n1_n2(1,4,prob_1))