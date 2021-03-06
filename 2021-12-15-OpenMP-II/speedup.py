import numpy as np
import matplotlib.pyplot  as plt

l1=[1,16]
l2=[1,16]
plt.plot(l1,l2,label=f"Teórico",linestyle="dashed")

x=np.arange(1,17)
yM=np.loadtxt("Man.txt")[:]
yA=np.loadtxt("Auto.txt")[:]

TM1=yM[0]
TA1=yA[0]

for i in x:
    yM[i-1]=(TM1/yM[i-1])
    yA[i-1]=(TA1/yA[i-1])

plt.scatter(x,yM,marker='+',color="k",label=f"Partición manual")
plt.scatter(x,yA,marker='x',color="r",label=f"Patición automática")

plt.xlabel("# Threads")
plt.ylabel("speedup")
plt.title("Speedup en función del número de threds (6 cores/12 threads)")
plt.legend()
plt.savefig("speedup.png",dpi=200)
