import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from math import *




ZERO=1e-1
def f(r):





    f1 =exp(-2.0*(abs(r)))

    return f1


x=np.linspace(-10.0,10.0,1000)


y1=np.zeros(len(x))

for i in range(len(x)):
        y1[i]=f(x[i])







plt.subplot(2,1,1)

plt.plot(x,y1,'r')
plt.title(r'$\psi=e^{-2r}$')

plt.axis([-10,10,0,1])
plt.xlabel("r")
plt.ylabel(r'$\psi$')
plt.legend([r'$r \epsilon[-10,10]$'])

plt.subplot(2,1,2)
plt.xlabel("r")
plt.ylabel(r'$\psi$')
plt.plot(x,y1,'r')
plt.legend([r'$r \epsilon[-2,2]$'])
plt.axis([-2,2,0,1])
plt.savefig('Pro3a.png')
