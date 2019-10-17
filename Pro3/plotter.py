import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from math import *

# Reading from file
f=open("int_data_2", "r")
lines =f.readlines()[0:]

n = 7
f.close()
Polar = np.zeros(n)
Kart =np.zeros(n)
N =np.zeros(n)
Analy =np.zeros(n)


# Reading out values
i = 0
for line in lines[0:]:
    words = line.split()
    N[i] = float(words[0])
    Kart[i] = float(words[1])
    Polar[i]=float(words[2])
    Analy[i]=0.192766
    i=i+1

# Plotting
# plt.hold(True)
plt.plot(N,Polar,'-*r')
plt.plot(N,Kart,'-*b')
plt.plot(N,Analy,'-g')
plt.xlabel('Polynomial degrees n',fontsize = 14)
plt.ylabel('Integral value',fontsize = 14)
plt.legend(['Polar', 'Cartesian','Analytical'])
plt.title("Numerical solution (Read and Blue) and exact analytical \n solution (Green)", fontsize = 14)
plt.savefig('Int.png')

"""
    plt.figure()
    ax = plt.subplot(111)
    #plt.subplot(3,1,1)
    plt.plot(xx, u, 'r')
    plt.plot(x,v,'b')
    #plt.axis([0, max(xx), min(u) - 0.5, max(u) + 0.5])   # [ xx_min, xx_max, u_cordmin, u_max]
    plt.xlabel('x',fontsize = 14)
    plt.ylabel('u(x)',fontsize = 14)

    ax.tick_params(labelsize = 14)
    plt.show()​
"""
