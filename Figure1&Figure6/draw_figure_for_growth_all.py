
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log, log2, ceil,sqrt

#ax.figure(figsize=(9, 7),dpi=1000)

#plt.style.use('fivethirtyeight')
#mpl.rcParams['text.usetex'] = True

fig, ax = plt.subplots(figsize=(10, 6),dpi=1000)
#ax.gcf().subplots_adjust(left = 0.04, right=0.05,top=0.05, bottom = 0.004)

#ax.subplots(constrained_layout=True)


def theo_DSP_cost(bitN, px, ux, pt, C1, C2):
    return log2(ux * log2(px))* C2 * bitN**2 + C1*pt*(1/2. * (log2(pt)**2) - log2(pt)+ux*log2(px)+2 )+C1*pt*(ux*log2(px))/log2(pt)*bitN*log2(bitN)


def theo_som22_cost(bitN, px, ux, pt, umax, C1, C2):
#    print(C1*pt*(log2(pt)**2) + umax * C2 * (bitN**2) + umax * C1 * pt * bitN * log2(bitN))
    return C1*pt*(log2(pt)**2) + umax * C2 * (bitN**2) + umax * C1 * pt * bitN * log2(bitN)





#bit_N = 1024
som22_cost_in_estimator = [5.34356e+11, 7.22133e+12, 1.42603e+13, 2.12992e+13, 2.83381e+13, 3.5377e+13 ]

som22_cost_in_estimator = som22_cost_in_estimator[:6]


DSP_cost_in_estimator = [6.29006e+11, 6.29006e+11, 1.12597e+12, 2.10275e+12, 2.10275e+12, 4.02609e+12 ]

DSP_cost_in_estimator = DSP_cost_in_estimator[:6]
DSP_cost_in_experiment = [404.33,749.87,1450.15,1454.84,2833.08,2891.40 ]
som22_cost_in_experiment = [394.85,6810.99,13579.94,20344.55,26452.08,34563.31 ]

#cost in cycles
DSP_cost_in_experiment = [_*1E9 for _ in DSP_cost_in_experiment]
som22_cost_in_experiment = [_*1E9 for _ in som22_cost_in_experiment]

bitN = 1024
pt = 1073741789 #3175639079
C1 = 0.5
C2 = 0.5
C3 = 6.41E-02

theo_DSP_costs = []
theo_som22_costs = []
for u in range(0,min(120,ceil(512*log(2)/log(3))),20):
    if(pt > 3**u):
        px = pt
        ux = 1
        umax = max(1,u)
    else:
        px = 3
        ux = u
        umax = u
    theo_DSP_costs.append(C3*theo_DSP_cost(bitN, px, ux, pt, C1, C2))
    theo_som22_costs.append(C3*theo_som22_cost(bitN, px, ux, pt, umax, C1, C2))
    
us = list(range(0,min(120,ceil(512*log(2)/log(3))),20))


#bit_N = 2048
#DSP_cost_in_estimator = [1.476291e+14, 1.476291e+14, 2.938020e+14, 5.857555e+14,
#                         5.857555e+14, 1.168941e+15, 1.168941e+15, 1.168941e+15,
#                         1.168941e+15, 2.333950e+15, 2.333950e+15, 2.333950e+15,
#                         2.333950e+15, 2.333950e+15, 2.333950e+15, 2.333950e+15,
#                         2.333950e+15, 4.661310e+15, 4.661310e+15, 4.661310e+15,
#                         4.661310e+15, 4.661310e+15, 4.661310e+15, 4.661310e+15,
#                         4.661310e+15, 4.661310e+15, 4.661310e+15, 4.661310e+15,
#                         4.661310e+15, 4.661310e+15, 4.661310e+15, 4.661310e+15,
#                         4.661310e+15]
#som22_cost_in_estimator = [1.472225e+14, 2.912719e+15, 5.823769e+15, 8.734818e+15,
#                           1.164587e+16, 1.455692e+16, 1.746797e+16, 2.037902e+16,
#                           2.329007e+16, 2.620112e+16, 2.911216e+16, 3.202321e+16,
#                           3.493426e+16, 3.784531e+16, 4.075636e+16, 4.366741e+16,
#                           4.657846e+16, 4.948951e+16, 5.240056e+16, 5.531161e+16,
#                           5.822266e+16, 6.113371e+16, 6.404476e+16, 6.695581e+16,
#                           6.986686e+16, 7.277791e+16, 7.568896e+16, 7.860001e+16,
#                           8.151105e+16, 8.442210e+16, 8.733315e+16, 9.024420e+16,
#                           9.315525e+16]



ax.plot(us,theo_DSP_costs,label="theo cost of our variant", color = "black", linewidth =2.0, linestyle="--", zorder = 1)
ax.plot(us,theo_som22_costs,label="theo IPP1 cost", color = "black", linewidth =2.0, linestyle="-.", zorder = 2)
ax.scatter(us,DSP_cost_in_estimator,label="cost of our variant in estimator", color = "orange", marker ='^', edgecolors = "black", zorder = 3)
ax.scatter(us,som22_cost_in_estimator,label="IPP1 cost in estimator", color = "red", edgecolors = "black", zorder = 4)
ax.scatter(us,DSP_cost_in_experiment,label="cost of our variant in experiment", color = "green", marker ='^', edgecolors = "black", zorder = 3)
ax.scatter(us,som22_cost_in_experiment,label="IPP1 cost in experiment", color = "blue", edgecolors = "black", zorder = 4)



ax.legend(fontsize=15)
plt.xlabel(r"$u_x$",fontsize=17)
plt.xticks(fontsize=17)
plt.ylabel(r"operation",fontsize=17)
plt.yticks(fontsize=17)
#ax.set_title(r"bit(N) = %d, P-1 = 2 $\times$ %d $\times$ $3^u$" %(bitN,pt))
#, $C_1$ = %d, $C_2$ = %d" %(bitN,pt,C1,C2))
plt.savefig("growth-fitness.png", bbox_inches='tight')

