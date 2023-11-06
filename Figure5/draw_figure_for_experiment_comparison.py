
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log, log2, ceil,sqrt
from scipy.interpolate import make_interp_spline
import numpy as np



#plt.style.use("ggplot")

fig, ax = plt.subplots(figsize=(10, 6),dpi=1000)

plt.grid(zorder = -1)
width=0.2

index=np.arange(1,9)
#print(index)
dsp_costs = [412.18, 748.37, 402.70, 401.68, 401.86, 407.08, 404.31, 402.60 ]
IPP1_costs = [5157.74,5110.87,4439.44,8827.70,8512.96,3483.75,4108.09,9498.75 ]
pol74_costs = [5573.02,1816.96,2666.92,6349.34,398.66,6682.96,5786.20,5164.79]
bis03_costs = [11912.68,3523.75,5364.30,13722.64,661.86,14506.06,12432.47,11038.81 ]
plt.bar(index,dsp_costs,width,tick_label = index, label='Our Variant with FM',zorder = 4)#color='orange',
plt.bar(index+width,IPP1_costs,width,label='IPP1 with FM',zorder = 3)
plt.bar(index+2*width,pol74_costs,width, label='Original Pollard\'s P-1',zorder = 2)
plt.bar(index+3*width,bis03_costs,width,label='Trivial Pollard\'s P-1',zorder = 2)


ax.legend(fontsize=15)
plt.xlabel(r"instance index",fontsize=19)
#plt.xlim(100,180)
plt.xticks(fontsize=19)
plt.ylabel(r"cost/s",fontsize=19)
plt.yticks(fontsize=19)
plt.ylim(None,18000)
#ax.set_title(r"bit(N) = %d, P-1 = 2 $\times$ %d $\times$ $3^u$" %(bitN,pt))
#, $C_1$ = %d, $C_2$ = %d" %(bitN,pt,C1,C2))
plt.savefig("pollardpm1-comparison.png",bbox_inches='tight')
plt.close()



