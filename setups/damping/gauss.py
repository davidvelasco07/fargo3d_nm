"""
Solve linear system using LU decomposition and Gaussian elimination
"""
from pylab import *
import numpy as np
from scipy.linalg import lu, inv

def gausselim(A,B):
    """
    Solve Ax = B using Gaussian elimination and LU decomposition.
    A = LU   decompose A into lower and upper triangular matrices
    LUx = B  substitute into original equation for A
    Let y = Ux and solve:
    Ly = B --> y = (L^-1)B  solve for y using "forward" substitution
    Ux = y --> x = (U^-1)y  solve for x using "backward" substitution
    :param A: coefficients in Ax = B
    :type A: numpy.ndarray of size (m, n)
    :param B: dependent variable in Ax = B
    :type B: numpy.ndarray of size (m, 1)
    """
    # LU decomposition with pivot
    pl, u = lu(A, permute_l=True)
    # forward substitution to solve for Ly = B
    y = np.zeros(B.size)
    for m, b in enumerate(B.flatten()):
        y[m] = b
        # skip for loop if m == 0
        if m:
            for n in range(m):
                y[m] -= y[n] * pl[m,n]
        y[m] /= pl[m, m]

    # backward substitution to solve for y = Ux
    x = np.zeros(B.size)
    lastidx = B.size - 1  # last index
    for midx in range(B.size):
        m = B.size - 1 - midx  # backwards index
        x[m] = y[m]
        if midx:
            for nidx in range(midx):
                n = B.size - 1  - nidx
                x[m] -= x[n] * u[m,n]
        x[m] /= u[m, m]
    return x


if __name__ == '__main__':

    ### Input values
    rhod1 = 1.0e-3; rhod2 = 1.0e-3; rhog=1.0e3;
    alp   = 0.005
    dt    = 1e-1

    ### Matrix
    A = np.array([[1+ dt*(rhod1/rhog*alp + rhod2/rhog*alp), -dt*rhod1/rhog*alp, -dt*rhod2/rhog*alp],
                  [-dt*alp                                , 1+dt*alp          , 0              ],
                  [-dt*alp                                , 0                 , 1+dt*alp       ]])

    ### Vy
    b = np.array([[2], [-3], [-6]])

    ### Initial condition for solvers
    x = x2 = xs = b
    t = 0

    save_diff   = [ log10(abs(x[0] - xs[0][0])/x[0]),]
    save_diff_2 = [ log10(abs(x2[0][0] - xs[0][0])/x2[0][0]),]
    save_diff_3 = [ log10(abs(x[0] - x2[0][0])/x[0]),]
    time        = [0,]
    xsol0       = [xs[0][0],]
    xsol1       = [xs[1][0],]
    xsol2       = [xs[2][0],]
    
    while(t<2000):

        ##### Solvers
        xn  = gausselim(A, x)
        x2n = np.linalg.solve(A, x2)

        s   = dt*alp/(1+dt*alp)
        As  = xs[0][0] + s*rhod1/rhog*xs[1][0] + s*rhod2/rhog*xs[2][0]
        Bs  = rhod1/rhog*s  + rhod2/rhog*s
        xsn = np.array( [ [As/(1+Bs)] ,
                        [ s*As/(1+Bs) + xs[1][0]/(1+dt*alp)],
                        [ s*As/(1+Bs) + x2[2][0]/(1+dt*alp)] ] )  
    
        x = xn; x2 = x2n; xs = xsn

        
        print("First component")        
        print("Gausselim vs Analitical: %1.16e"% (x[0]-xs[0][0])     )
        print("LinAlg    vs Analitical: %1.16e"% (x2[0][0]-xs[0][0]) )

        print("Second component")        
        print("Gausselim vs Analitical: %1.16e"% (x[1]-xs[1][0])     )
        print("LinAlg    vs Analitical: %1.16e"% (x2[1][0]-xs[1][0]) )

        print("Third component")        
        print("Gausselim vs Analitical: %1.16e"% (x[2]-xs[2][0])     )
        print("LinAlg    vs Analitical: %1.16e"% (x2[2][0]-xs[2][0]) )

        save_diff.append(   (abs((x[0]-xs[0][0]))/x[0]) )
        save_diff_2.append( (abs((x2[0][0]-xs[0][0]))/x2[0][0]) )
        save_diff_3.append( (abs((x[0]-x2[0][0]))/x[0]) )

        t+=dt
        
        time.append(t)

        ### Solutions
        xsol0.append( xs[0][0] )
        xsol1.append( xs[1][0] )
        xsol2.append( xs[2][0] )

    plot(time,  save_diff  , 'x-', label='Gauss-Elimination  vs Analitical') 
    plot(time,  save_diff_2, '-', label='LinAlg  vs Analitical' ) 

    save("python_data_vy", array(save_diff))
    save("python_data_time", array(time))
    
    #plot( save_diff_3, '-', label='Gauss-Elimination  vs LinAlg' )
    xscale('log')
    yscale('log')
    ylabel('Delta')
    xlabel('Time')
    print()
    print()
    print( 'Max difference of Gauss-Elimination  vs LinAlg %1e18e'% max(save_diff_3) )
    legend()
    show()

    clf()
    plot(time, xsol0)
    plot(time, xsol1)
    plot(time, xsol2)
    show()
