"""
Solve linear system using LU decomposition and Gaussian elimination
"""
from pylab import *

### Input values
rhod1 = 1.0e-3; rhod2 = 1.0e-3; rhog=1.0e3;
alp   = 0.005
dt    = 1e-1

### Matrix
A = np.array([[1+ dt*(rhod1/rhog*alp + rhod2/rhog*alp), -dt*rhod1/rhog*alp, -dt*rhod2/rhog*alp],
              [-dt*alp                                , 1+dt*alp          , 0              ],
              [-dt*alp                                , 0                 , 1+dt*alp       ]])

### Vy
b1 = np.array([[2], [-3], [-6]])
b2 = np.array([[2], [-3+4e-16], [-6]])

### Initial condition for solvers
x1 = b1
x2 = b2
t = 0

save_diff = [ abs(x1[1][0] - x2[1][0])/abs(x1[1][0]),]
time      = [0,]

while(t<=2000):

    ##### Solvers
    x1n = linalg.solve(A, x1)
    x2n = linalg.solve(A, x2)

    x1 = x1n; x2 = x2n;

    save_diff.append( abs(x1[1][0]-x2[1][0])/abs(x1[1][0]) )

    t+=dt
    
    time.append(t)

plot(time,  save_diff  , 'x-') 

xscale('log')
yscale('log')
ylabel('Delta')
xlabel('Time')
print()
print()
print( 'Max difference of Gauss-Elimination  vs LinAlg %1e18e'% max(save_diff) )
show()
