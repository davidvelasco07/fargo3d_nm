from pylab import *

def stream(x,A,B,C):
    return sqrt(-C/A+B/(4.0*A)*arctan(x**2))

def compute_C(A, B, ymax, n):
    C_max = A*abs(ymax)**2-B*pi/(8.0*A)
    return np.linspace(-C_max,C_max,n)
    
A = 1.0
B = 1.0
xmin = -5
xmax = 5
ymin = -2
ymax = 2

x = np.linspace(-5,5,10000)
c = compute_C(A,B,ymax,50)

figure(figsize=(15,5))

separ = stream(x,A,B,0)

for c0 in c:
    s = stream(x,A,B,c0)
    plot(x,s,'k',linewidth=2)
    plot(x,-s,'k',linewidth=2)
    plot(x,separ,'r',linewidth=3)
    plot(x,-separ,'r',linewidth=3)
    
xlim(xmin,xmax)    
ylim(ymin,ymax)

show()
