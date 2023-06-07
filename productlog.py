from scipy.special import lambertw
from math import sqrt,log
import matplotlib.pyplot as plt 
import numpy
from scipy import optimize


def f_1(x, A, B, C):
    return A * x**2 + B * x + C

def compute_R2(actual,predict):
    corr_matrix = numpy.corrcoef(actual, predict)
    corr = corr_matrix[0,1]  #相关系数
    R_sq = corr**2
    return R_sq



def asym_w(x):
    L1 = log(x)
    L2 = log(log(x))
    return L1 - L2 + L2/L1 #+ L2*(-2+L2)/(2*L1**2) + L2*(6-9*L2+2*L2**2)/(6*L1**3)+L2*(-12 + 36*L2 - 22 * L2 **2 + 3* L2**3)/(12*L1**4) 

def fw(w,l):
    floor_w = int(w)
    fw1 = 2**(floor_w-1)-1+(l-1)/(floor_w+1)
    fw2 = 2**(floor_w+1-1)-1+(l-1)/(floor_w+1+1)
    if fw1 <= fw2:
        return floor_w,fw1
    else:
        return floor_w+1,fw2



#picture of w

max_num = int(1.2e9)
start = int(1e7)
gap = int(1e7)
ls = [_ for _ in range(start,max_num,gap)]
ws = []
intws = []
asym_ws =[]
fws = []
w_dic ={}

for l in ls:
    w = 2*lambertw(sqrt((l-1)*log(2)))/log(2)-1
    ws.append(w)
    intw,fw_ = fw(w,l)
    intws.append(intw)
    fws.append(fw_)
    asym_ws.append(2*asym_w(sqrt((l-1)*log(2)))/(log(2))-1)
    if intw not in w_dic:
        w_dic[intw] = [l]
    else:
        w_dic[intw] += [l]
        
for w in w_dic:
    print(w, min(w_dic[w]),max(w_dic[w]))
    
    

with plt.style.context(['science','ieee']):
    fig, ax = plt.subplots(figsize=(4,3),dpi=600)   
    R2 = compute_R2(ws,asym_ws)
    ax.plot(ls, ws, label=r"$w = \frac{2W(\sqrt{(l-1)\cdot \ln 2})}{\ln 2} - 1$")
    ax.plot(ls, intws, label=r"$\lceil w \rfloor$")
    ax.plot(ls, asym_ws, label=r"Asymptotic $w = L_1 - L_2 +\frac{L_2}{L_1}$, $R^2 = %.2f$" %R2, color = "orange")
    # plt.title(r'$w = \frac{2W(\sqrt{(l-1)\cdot \ln 2})}{\ln 2} - 1$')
    ax.legend()
    ax.set(xlabel=r'$l$')
    ax.set(ylabel='$w$')
    fig.savefig(r'w.png')
    
    
with plt.style.context(['science','ieee']):
    fig, ax = plt.subplots(figsize=(4,3),dpi=600)   
    
    ax.plot(ls, fws, label='$f(w(l))_{\min}$')
    
    
    # plt.title(r'$w = \frac{2W(\sqrt{(l-1)\cdot \ln 2})}{\ln 2} - 1$')
    
    #linear fitting
    A1, B1, C1= optimize.curve_fit(f_1, ls, fws)[0] 
    # x1 = numpy.arange(start,max_num,gap)
    x1 = numpy.array(ls)
    y1 = A1 * x1**2 + B1 * x1 + C1
    
    
    R = compute_R2(fws,y1)
    plt.plot(x1, y1, color="blue", linestyle="--", label=" y = %1.3e $l^2$ + %1.3e$l$ + %3.3f, $R^2$ = %3.4f" %(A1,B1,C1,R))
    
    
    plt.title("Modular Powering Cost (in theory) -- bitsize")
    
    ax.legend()
    ax.set(xlabel=r'$l$')
    ax.set(ylabel='$f(w(l))_{\min}$')
    
    fig.savefig(r'fw_theo.png')
    
    
ls = [1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824 ]
fws = [1.3889, 2.73596, 5.51345, 11.7506, 22.8856, 46.698, 90.972, 183.055, 350.514, 702.463, 1412.65]
    
with plt.style.context(['science','ieee']):
    fig, ax = plt.subplots(figsize=(4,3),dpi=600)   
    
    # ax.plot(ls, fws, label='$f(w(l))_{\min}$')
    ax.scatter(ls, fws, label='Experiment in bit(N) = bit(a) = 2048',color = "red", marker='.')
    
    # plt.title(r'$w = \frac{2W(\sqrt{(l-1)\cdot \ln 2})}{\ln 2} - 1$')
    
    #linear fitting
    A1, B1, C1= optimize.curve_fit(f_1, ls, fws)[0] 
    # x1 = numpy.arange(start,max_num,gap)
    x1 = numpy.array(ls)
    y1 = A1 * x1**2 + B1 * x1 + C1
    
    
    R = compute_R2(fws,y1)
    plt.plot(x1, y1, color="blue", linestyle="--", label=" y = %1.3e $l^2$ + %1.3e$l$ + %3.3f, $R^2$ = %3.4f" %(A1,B1,C1,R))
    
    
    
    
    plt.title("Modular Powering Cost (in experiment) -- bitsize")
    
    ax.legend()
    ax.set(xlabel='bit size of the exponent')
    # ax.set(ylabel='$f(w(l))_{\min}$')
    ax.set(ylabel='Modular Exponent Cost/s')
    
    fig.savefig(r'fw_exp.png')