import numpy as np
import matplotlib.pyplot  as plt

x=np.loadtxt("scaling.txt")[:,0]
y=np.loadtxt("scaling.txt")[:,2]
plt.scatter(x,y)
plt.xlabel("# Threads")
plt.ylabel("parallel eficiency")
plt.title("Parallel eficiency (in a computer with 12 threads)")
plt.savefig("parallel.png",dpi=200)
