import numpy as np
import matplotlib.pyplot  as plt

x=np.loadtxt("scaling.txt")[:,0]
y=np.loadtxt("scaling.txt")[:,1]
plt.scatter(x,y)
plt.xlabel("# Threads")
plt.ylabel("speedup")
plt.title("Speedup (in a computer with 12 threads)")
plt.savefig("speedup.png",dpi=200)
