from math import*
from numpy import linspace
from numpy.linalg import*
from pylab import *
import scipy.integrate


L = 5*10**-10
a = 10*1.6*10**-19
hbar = 1.05*10**-34
mass = 9.1*10**-31    




def simpson(f, a, b, h):

    n = int(abs(b - a)/h)
    if n % 2 != 0:
        n -= 1

    sum1 = 0
    s3 = 0
    for i in range(1, n, 2):
        s3 += f(a + (i*h))
    s3 *= 4
    sum1 += s3
    sum2 = 0
    s4 = 0
    for i in range(2, n-1, 2):
        s4 += f(a + i*h)
    s4 *= 2
    sum2 = s4
    approx = h/(3.0)*(f(a) + f(b) + sum1 + sum2)
    return approx


def V(x):
    return 0 if 0<x<1 else 100000


def Hmn(m, n, V):
    return (2/L)*simpson(lambda x: sin(m*pi*x/L)*((n**2*pi**2*hbar**2/(2*mass*L**2))*sin(n*pi*x/L)+V(x)*sin(n*pi*x/L)), 0, L, .01*L)







def phi(n, x):
    return sqrt(2/L)*sin(n*pi*x/L)


def generate_psi(phi, col):
    return lambda x: sum(col[sam] * phi(sam, x) for sam in range(len(col)))

def An(x, n):
    if 0 <= x < L/2:
        return sqrt(12/L**3) * simpson(lambda x: x * psis[n](x), 0, L/2, .01*L)
    elif L/2 <= x <= L:
        return sqrt(12/L**3) * simpson(lambda x: (L-x) * psis[n](x), L/2, L, .01*L)
    else:
        return


def An(x, n):
    return sqrt(12/L**3) * simpson(lambda x: x * psis[n](x), 0, L/2, .01*L) + sqrt(12/L**3) * simpson(lambda x: (L-x) * psis[n](x), L/2, L, .01*L)


def Psi(x, t):
    return sum([An(x, n) * estuff(t, n) * psis[n](x) for n in range(len(psis))])


def estuff(t, n):
    return e ** (-(0+1j)*energy[n] * t*(1/hbar))




def generate(t):
    i = t[0]
    t = t[1]
    print(f"Starting Frame {i}")
    ans = [abs(Psi(x, t))**2 for x in linspace(0, L, 100)]
    # img.append(plot(linspace(0, L, 100), ans, label="$\Psi(x, t=%s)$" % t))
    # integ = scipy.integrate.simps(ans, linspace(0, L, 100))
    # print(integ)

    print(f"Finished Frame {i}")
    plt.plot(linspace(0, L, 100), ans)
    plt.savefig(
        "psis_frames_new/frame{}.png".format(str(i+1).zfill(7)))
    clf()


H = zeros([20, 20], float)
for m in range(0, 20):
    for n in range(0, 20):
        H[m][n] = Hmn(m+1, n+1, V)
energy, Anp = eigh(H)
energy *= 6.2*10**18

psis = []
for vec in Anp:
    psis.append(generate_psi(phi, vec))


energy, Anp = eigh(H)
energy *= 6.2*10**18


if __name__ == '__main__':
    import multiprocessing
    pool = multiprocessing.Pool(12)


    NUM_FRAMES = 1000
    
    
    t = linspace(0, 5, NUM_FRAMES)
    pool.map(generate,  [(i, t[i]) for i in range(len(t))])

