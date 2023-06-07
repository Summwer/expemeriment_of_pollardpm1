from scipy.special import lambertw
from math import sqrt,log
import matplotlib.pyplot as plt 
import numpy
from scipy import optimize


def f_nm(x, A, B):
    return A * 4**(x) + B #*(x**2) 

def f_fm(x, A, B):
    return A * 2**(x-1) * (x**2) + B #*(x**2) 

def compute_R2(actual,predict):
    corr_matrix = numpy.corrcoef(actual, predict)
    corr = corr_matrix[0,1]  #相关系数
    R_sq = corr**2
    return R_sq
    
    
ls = [21, 22, 23, 24, 25, 26]
nms = [2.3, 10.22, 42.07, 142.46, 605.49, 2378.83]
fms = [0.07, 0.19, 0.6, 0.96, 2.55, 4.27]
    
with plt.style.context(['science','ieee']):
    fig, ax = plt.subplots(figsize=(4,3),dpi=600)   
    
  
    ax.scatter(ls, nms, label='normal_mul cost',color = "orange", marker='.', zorder=4)
    ax.scatter(ls, fms, label='fast_mul cost',color = "red", marker='*', zorder=3)
    
    
    A1, B1= optimize.curve_fit(f_nm, ls, nms)[0] 

    x1 = numpy.array(ls)
    y1 = A1 * 4**(x1)  + B1 #*(x1**2)
    R = compute_R2(nms,y1)
    
    gap = 0.05
    xs = [ls[0]+ gap*i for i in range(int(len(ls)/gap)) if ls[0]+ gap*i <=ls[-1]]
    x1 = numpy.array(xs)
    y1 = A1 * 4**(x1)  + B1 #*(x1**2)
    
    plt.plot(x1, y1, color="gray", linestyle="--", label=" y = %1.3e $4^n $ + %3.4f, $R^2$ = %3.4f" %(A1,B1,R)) #n^2
    
    
    A2, B2= optimize.curve_fit(f_fm, ls, fms)[0] 

    x2 = numpy.array(ls)
    y2 =  A2 * 2**(x2) *(x2**2)  + B2 #*(x1**2)
    R = compute_R2(fms,y2)
    
    gap = 0.05
    x2 = numpy.array(xs)
    y2 = A2 * 2**(x2) *(x2**2)  + B2 #*(x1**2)
    
    plt.plot(x2, y2, color="brown", linestyle="--", label=" y = %1.3e $2^{n-1} n^2 $ + %3.4f, $R^2$ = %3.4f" %(A2,B2,R)) #n^2
    
    

    plt.title("Cost Comparison Between normal_mul and fast_mul in Prime List")
    
    ax.legend()
    ax.set(xlabel='bit size of the maximal prime multiplier')
    # ax.set(ylabel='$f(w(l))_{\min}$')
    ax.set(ylabel='cost/s')
    
    fig.savefig(r'mul_cost.png')