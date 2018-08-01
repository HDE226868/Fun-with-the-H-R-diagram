import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

m_min = 0.01
m_max = 100

N = float(input("Total number of stars: ")) # Number of stars in sample

def kroupa_index(m):
    """Returns the power law index in a given range, with m in solar masses"""
    if m < .08:
        a = -0.3
    elif 0.08 <= m and m < 0.5:
        a = -1.3
    else:
        a = -2.3
    return a

# normalization constant
xi0 = 1/(integrate.quad(lambda x: x**(kroupa_index(x)), m_min, m_max)[0])

def kroupa(m1,m2):
    """Integrates IMF from m1 to m2, with masses in solar masses"""
    return N*xi0*(integrate.quad(lambda x: x**(kroupa_index(x)), m1, m2)[0])


def Mmin(t):
    """Determines the lowest-mass star to enter main sequence at time t"""
    return (t/10**(7))**(-2/5)

def Mmax(t):
    """Determines the highest-mass star to enter main sequence at time t"""
    return (t/10**(10))**(-2/5)

def NumPMS(t):
    """Calculates the number of pre-main sequence stars"""
    return kroupa(m_min,Mmin(t))

def NumMS(t):
    """Calculates the number of main sequence stars"""
    return kroupa(Mmin(t),Mmax(t))

def NumPoMS(t):
    """Calculates the number of post-main sequence stars and stellar remannts"""
    return N - NumPMS(t) - NumMS(t)

M_range = np.linspace(m_min,m_max)

plt.loglog(M_range,[xi0*M**(kroupa_index(M)) for M in M_range],'k')
plt.xlabel("Mass (solar masses)")
plt.ylabel(r"$\xi(m)$")
plt.title("Kroupa IMF")
plt.show()

Tot = [kroupa(m_min,M) for M in M_range]
plt.loglog(M_range,Tot,'k')
plt.xlabel("Mass (solar masses)")
plt.ylabel("Number of stars of mass less than M")
plt.title("Cumulative mass function of Kroupa IMF")
plt.show()

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

T = np.linspace(10**5,10**(11),1000)

NumsPMS = [NumPMS(t) for t in T]
ax1.semilogx(T,NumsPMS,'k')
ax1.set_ylabel("Number of stars")
ax1.set_title("Pre-main sequence evolution")

NumsMS = [NumMS(t) for t in T]
ax2.semilogx(T,NumsMS,'k')
ax2.set_ylabel("Number of stars")
ax2.set_title("Main sequence evolution")

NumsPoMS = [NumPoMS(t) for t in T]
ax3.semilogx(T,NumsPoMS,'k')
ax3.set_xlabel("Time (years)")
ax3.set_ylabel("Number of stars")
ax3.set_title("Post-main sequence evolution")

plt.show()