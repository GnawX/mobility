import numpy as np
from math import pi,sqrt
import matplotlib.pyplot as plt

#input
m   = 0.25
eps = 6.5
epsr = 30
z   = 1  # charge
omega = 0.0126 # eV

# parameters SI units
e    =  1.602176634e-19
m0   =  9.10938356e-31
kb   =  1.38064852e-23
kbe  =  8.617333262145e-5
hbar =  1.054571817e-34
eps0 =  8.8541878128e-12
ha   =  27.211
#
alpha = sqrt(m/2/(omega/ha))*(1/eps-1/epsr)
eps = eps*eps0
m = m*m0

#Brooks-Herring mobility 
def bk_mu(T,ni,n):
    b = 24*m*eps/n*(kb*T/e/hbar)**2
    a = np.log(1+b) - b/(1+b)
    return 8*(4*pi*eps)**2*kb*T*sqrt(2*kb*T/pi/m)/ \
           (pi*e**3*z**2*ni*a)

#polar phonon limited mobility (ACS Energy Lett. 2019, 4, 456−463)
def ph_mu(T):
    return (0.52*(omega/kbe/T)**3.3+0.34)*e*hbar/alpha/m/kb/T

n = 1e+10*1e+6
  
n1 = 1e+14*1e+6
n2 = 1e+15*1e+6
n3 = 1e+16*1e+6
n4 = 1e+17*1e+6
n5 = 1e+18*1e+6



T = np.arange(5,501,5)
mu1 = np.zeros(len(T))
mu2 = np.zeros(len(T))
mu3 = np.zeros(len(T))
mu4 = np.zeros(len(T))
mu5 = np.zeros(len(T))
mup = np.zeros(len(T))

for i in range(len(T)):
    mu1[i] = bk_mu(T[i],n1,n)*1e+4
    mu2[i] = bk_mu(T[i],n2,n)*1e+4
    mu3[i] = bk_mu(T[i],n3,n)*1e+4
    mu4[i] = bk_mu(T[i],n4,n)*1e+4
    mu5[i] = bk_mu(T[i],n5,n)*1e+4
    mup[i] = ph_mu(T[i])*1e+4
    mu1[i] = 1/(1/mu1[i]+1/mup[i])
    mu2[i] = 1/(1/mu2[i]+1/mup[i])
    mu3[i] = 1/(1/mu3[i]+1/mup[i])
    mu4[i] = 1/(1/mu4[i]+1/mup[i])
    mu5[i] = 1/(1/mu5[i]+1/mup[i])

plt.plot(T,mu1,label=r'$10^{14} cm^{-3}$')
plt.plot(T,mu2,label=r'$10^{15} cm^{-3}$')
plt.plot(T,mu3,label=r'$10^{16} cm^{-3}$')
plt.plot(T,mu4,label=r'$10^{17} cm^{-3}$')
plt.plot(T,mu5,label=r'$10^{18} cm^{-3}$')
plt.plot(T,mup,label='phonon')

plt.yscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel(r'Mobility ($cm^2V^{-1}s^{-1}$)')
plt.legend()
plt.show()
