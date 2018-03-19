from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor, trim_zeros, append
from numpy.linalg import inv
from time import time
from sys import argv

def block(x):
    # preliminaries
    n = len(x); d = int(log2(n)); s, gamma = zeros(d), zeros(d);
    mu = mean(x); t0 = time()

    # estimate the auto-covariance and variances 
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])
   
    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print "Warning: Use more data"
    ans = s[k]/2**(d-k)
    print "Runtime: %g sec" % (time()-t0); print "Blocking Statistics :"
    print "average            iterations      std. error"
    print "%8g %20g %15g" % (mu, k, ans**.5)
    return ans

# input data must be a power of two
#x = loadtxt(argv[1])
#x = x[int(0.75*len(x)):len(x)]
eq = 0.25


x = loadtxt("/home/mathisre/Dropbox/Programs/VMC-bosons2/with_MPI/ground_state_0.dat")

x0 = loadtxt("/home/mathisre/Dropbox/Programs/VMC-bosons2/with_MPI/ground_state_0.dat")
x1 = loadtxt("/home/mathisre/Dropbox/Programs/VMC-bosons2/with_MPI/ground_state_1.dat")
x2 = loadtxt("/home/mathisre/Dropbox/Programs/VMC-bosons2/with_MPI/ground_state_2.dat")
x3 = loadtxt("/home/mathisre/Dropbox/Programs/VMC-bosons2/with_MPI/ground_state_3.dat")

x0 = x0[int((1-eq)*len(x)):len(x)]
x1 = x1[int((1-eq)*len(x)):len(x)]
x2 = x2[int((1-eq)*len(x)):len(x)]
x3 = x3[int((1-eq)*len(x)):len(x)]
print len(x)

x5 = append(x0,x1)
x4 = append(x2,x3)
x = append(x5,x4)

x = loadtxt("/home/mathisre/Dropbox/Programs/VMC-bosons2/build-no-matr-Desktop_Qt_5_10_1_GCC_64bit-Release/energy.dat")
x = x[2*10**6-2**20:2*10**6]
ans = block(x)



