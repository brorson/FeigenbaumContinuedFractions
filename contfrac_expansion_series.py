import gmpy2
from gmpy2 import mpfr
import numpy as np

import matplotlib.pyplot as plt


def compute_expansion_coeffs(x, N):
    """
    This computes the continued fraction coeffs
    of the real number x.  It computes the first
    N coeffs.
    """
    a = np.zeros((N,), dtype=object)

    a[0] = np.floor(x)
    xi = x-a[0]
    b = mpfr("1.0")/xi

    for i in range(1,N):
        a[i] = np.floor(b)
        xi = b - a[i]
        b = mpfr("1.0")/xi
    return a


def print_expansion_coeffs(a, N):
    for i in range(len(a)):
        print(a[i])

        
def plot_expansion_coeffs(a, N):
    """
    This plots the expansion coeffs
    console.
    """
    plt.plot(a, 'r.')
    plt.yscale('log')
    plt.show()
    

def compute_x(a, N):
    """
    Given expansion coeffs a and the number of coeffs N,
    this computes the real number x from the continued
    fraction expansion.  It starts at the end of the 
    expansion and works backwards.
    """
    x = mpfr(a[N-1])
    for i in range((N-2), -1, -1):
        x = a[i] + mpfr("1.0")/x
    return x

def write_test_number(x, N, filename):
    """
    Use this to write test numbers like sqrt(2) to
    a file to verify the program.  Use it like this:
    delta_expansion_series.write_test_number(5,100,'sqrt2.txt').
    """
    gmpy2.get_context().precision=3*N
    FID = open(filename, 'w')
    y = "%s" % gmpy2.sqrt(x)
    for i in range(len(y)):
        FID.write("%s" % y[i])
    FID.close()


def main():
    """
    Main fcn.  This reads in the number to analyze, then
    loops over number of terms in continued fraction 
    expansion and computes the Khinchin coeff for each.
    """
    gmpy2.get_context().precision=10001

    FID = open('alpha.txt', 'r')
    str = FID.read().replace('\n','').replace('-','')
    FID.close()
    
    #Nmax = len(str)
    Nmax = 9999
    #print(str[0:10+2])
    
    K = np.zeros([Nmax,1])
    ldiff = np.zeros([Nmax,1])
    for n in range(1, Nmax):
        x = mpfr(str[0:n+2]) # Add 2 to compensate for "4."
        a = compute_expansion_coeffs(x,n)

        # Make sure everything adds up
        xcomp = compute_x(a, n)
        diff = abs(x-xcomp)
        ldiff[n] = gmpy2.log10(diff)
        # Change this back to 2.0 eventually
        tol = 1000000.0*np.power(mpfr(10.0), mpfr(-n))

        # print_expansion_coeffs(a, n)
        K[n] = np.power(np.prod(a),(mpfr("1.0")/mpfr(n))) 
        print("n = %d, K = %f, diff = %s, tol = %s" % (n, K[n], '{0:.5e}'.format(diff),'{0:.5e}'.format(tol) ))
        
        if (n > 20 and diff > tol):
            print("!!!!!Failed round trip!!!!!")
            print("diff = ")
            print(diff)
            print("tol = ")
            print(tol)
            #break
        
    # Plot convergence to Khinchin constant
    plt.figure(1)
    plt.plot(K, 'r.', markersize=1)
    plt.plot([1,Nmax],[2.6854520010, 2.6854520010], 'b', linewidth=1 )
    plt.title('Khinchin constant vs. number of convergents')
    plt.ion()
    plt.show()

    plt.figure(2)
    plt.plot(a, 'r.')
    plt.yscale('log')
    plt.title('Partial fraction coeffs')
    plt.ion()
    plt.show()

    plt.figure(3)
    plt.plot(ldiff, 'r.')
    plt.title('log diff')
    plt.ion()
    plt.show()


#======================================
if (__name__ == "__main__"):
    main()

    
