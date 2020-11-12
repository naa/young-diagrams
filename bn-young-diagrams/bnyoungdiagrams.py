"""This module can be used to find most probable rotated generalized
Young diagrams for the decomposition of the N-th tensor power of
so(2n+1) spinor representation, to plot the rotated generalized Young
diagrams and to plot the limit shape.

"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from scipy import integrate
from matplotlib import transforms
import math
from numpy.lib import scimath 

def rho(x,c):
    """Charge density or density of steps of the diagram in the limit
    N,n\to\infty with c=(N+2n-1)/n

    """
#    return (1-np.sign(x-c/2))/2*((1-np.sign(c-4))/4-np.sign(c-4)/(4*np.pi)*
#                                 (np.imag(scimath.log(-8-c*(x-4)+np.abs(c-4)*scimath.sqrt(x*x-2*c+4))+
#                                  scimath.log(-8+c*(x+4)+np.abs(c-4)*scimath.sqrt(x*x-2*c+4)))-np.pi))
#    
    res= (1-np.sign(np.abs(x)-c/2))/2*np.real((1-np.sign(c-4))/4+np.sign(c-4)/(4*np.pi)*
                                              np.where(np.abs(x)<np.sqrt(2*c-4),
                                                       (np.arctan((-8-c*(x-4))/(np.abs(c-4)*scimath.sqrt(-x*x+2*c-4)))+
                                                        np.arctan((-8+c*(x+4))/(np.abs(c-4)*scimath.sqrt(-x*x+2*c-4)))),0))
    return res


def plot_limit_shape(n,N):
    """
    Plot limit shape for so(2n+1) rank n and tensor power N, c=(N+2n-1)/n
    """
    cc=(N+2*n-1)/n
    x=np.arange(0, cc/2+1/2,(cc+1)/400) 
    y=rho(x,cc)
    y = np.where(np.isnan(y), 0, y)
    ft=1+integrate.cumtrapz(1-4*y, x,initial=0)
    plt.plot(x,ft)


def plot_rotated_diagram(al):
    """
    Plot rotated generalized Young diagram from the list of coordinates [a_1,a_2,...,a_n]
    """
    n=len(al)
    for i in range(n):
        for j in range((al[i]-2*(n-i)+1+al[-1]-1)//2):
            if j<al[-1]-1:
                plt.plot([(2*n-2*i+j)/(2*n),
                          (2*n-2*i+j+1)/(2*n),
                          (2*n-2*i+j-1)/(2*n),
                          (2*n-2*i+j-2)/(2*n),
                          (2*n-2*i+j)/(2*n)],
                         [(2*i+j)/(2*n),
                          (2*i+j+1)/(2*n),
                          (2*i+j+3)/(2*n),
                          (2*i+j+2)/(2*n),
                          (2*i+j)/(2*n)],'gray')
            else:
                plt.plot([(2*n-2*i+2*j-al[-1]+1)/(2*n),
                          (2*n-2*i+2*j-al[-1]+3)/(2*n),
                          (2*n-2*i+2*j-al[-1]+1)/(2*n),
                          (2*n-2*i+2*j-al[-1]-1)/(2*n),
                          (2*n-2*i+2*j-al[-1]+1)/(2*n)],
                         [(2*i+2*j-al[-1]+1)/(2*n),
                          (2*i+2*j-al[-1]+3)/(2*n),
                          (2*i+2*j-al[-1]+5)/(2*n),
                          (2*i+2*j-al[-1]+3)/(2*n),
                          (2*i+2*j-al[-1]+1)/(2*n)],'gray')




def log(x):
    """
    Natural logarithm that returns -\infty for x<0
    """
    return math.log(x) if x > 0 else -math.inf

def logstirling(n):
    """
    Log of stirling approximation, uses lgamma but returns \infty for n<0
    """
    if n<0:
        return math.inf
    return math.lgamma(n+1)    

def logbnmultiplicity(a,  n,  N):
    """
    Computes the logarithm of the multiplicity in the tensor power decomposition
    """
    res=0
    for i in range(len(a)):
        res+=logstirling(N+2*i)-2.0*i*log(2)-logstirling((N+a[i]+2*n-1)/2)-logstirling((N-a[i]+2*n-1)/2)+log(a[i])
        for j in range(i):
            res+=log(a[i]*a[i]-a[j]*a[j])
    return res

def logbndimension(a, n):
    """
    Computes the logarithm of the dimension of an irrep with 
    the highest weight with the coordinates [a_1,...,a_n]
    """
    res=(-n*n+2*n)*log(2)+logstirling(n)
    for i in range(n):
        res+=log(a[i])-logstirling(2*n-2*i)
        for j in range(i):
            res+=log(a[i]*a[i]-a[j]*a[j])
    return res


def logthetameasure(a, k, N, theta):
  return -theta+k*log(theta)-logstirling(k)+logbnmultiplicity(a,k,N)+logbndimension(a,k)

def logmun(a, n, N):
    """
    Computes the logarithm of the probability
    """
    res=logbnmultiplicity(a,n,N)+logbndimension(a,n)-n*N*log(2)
    return res


def findmax(n, N):
    """
    Finds most probable configuration for the rank n and tensor power N, returns
    coordinates [a_1,...,a_n]
    """
    a=[]
    for i in range(n):
        a.append(2*i+1+N%2)
    cont=True
    while cont:
        cont=False
        for i in range(n):
            logp=logmun(a,n,N)
            a[i]+=2
            newlogp=logmun(a,n,N)
            while newlogp>logp:
                cont=True
                logp=newlogp
                a[i]+=2
                newlogp=logmun(a,n,N)
            a[i]-=2
    a.reverse()
    return a


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="When used as a standalone program, this module computes the most probable diagram for the given values of algebra rank and tensor power. It prints the coordinates of {a} of the rotated diagram to the standard output and saves the plot of the diagram and the limit shape to a pdf file.")
    parser.add_argument("-n", "--rank", help="Rank n of the algebra so(2n+1), default n=30",type=int)
    parser.add_argument("-N", "--power", help="Tensor power of spinor representation of the algebra so(2n+1), default N=60",type=int)
    parser.add_argument("-f", "--file", help="Diagram output filename (pdf), default name is fig.pdf",
                        type=str)
    args = parser.parse_args()
    rank=30
    power=60
    filename="fig.pdf"
    if args.rank:
        if args.rank>100:
            print("Rank value too big, it will take too long. Modify the code if you really mean it.")
            exit(1)
        if args.rank<=0:
            print("Rank should be positive")
            exit(1)
        rank=args.rank
    if args.power:
        if args.power>10000:
            print("Tensor power value is too big, it will take too long. Modify the code if you really mean it.")
            exit(1)
        if args.power<=0:
            print("Tensor power should be positive")
            exit(1)
        power=args.power
    if args.file:
        filename=args.file
    maxconf=findmax(rank,power)
    print(maxconf)
    plot_rotated_diagram(maxconf)
    plot_limit_shape(rank,power)
    plt.savefig(filename)
