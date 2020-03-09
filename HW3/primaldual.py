import numpy as np
from matplotlib import pyplot as plt



f0 = lambda x: 1 + x**2
L = lambda x,v: (1 + v)*x**2 - 6*v*x + (1 + 8*v) #f0(x) + v*(x-2)*(x-4)

plt.figure(1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

x = np.linspace(2,4,1000)
plt.plot(x, f0(x), lw = 2, label = 'Objective $f_0(x)$')

plt.plot(2, 5, marker = '*', markersize = 15, label = 'Optimal point')

# The Lagrange multipliers to plot.
vv = [1, 2, 3]
for v in vv:
    plt.plot(x, L(x,v), label = 'Lagrangian $L(x,\lambda)$, $\lambda = {:0.1f}$'.format(v))

plt.xlabel('$x$')
plt.ylabel('Function value')
plt.title('Objective function and the Lagrangian over the feasible set')
plt.legend()


plt.figure(2)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
v = np.linspace(0,10,1000)
g = -9*v**2/(1+v) + 8*v + 1

plt.plot(v, g, lw = 2)
plt.xlabel('Dual variable $\lambda$')
plt.ylabel('Dual function $g(\lambda)$')
plt.title('Lagrange dual function')


p = lambda u: u - 6*np.sqrt(u + 1) + 11

plt.figure(3)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

uu = np.linspace(-1, 10, 1000)
plt.plot(uu, p(uu), lw = 2)
plt.xlabel('Parameter $u$')
plt.ylabel('Optimal value $p^*(u)$')
plt.title('Optimal value of the parametric problem')

plt.show()
