#---------------------
# CubeSat thermal model
#
#Author: Nicolas Beaudoin
#Course: MCG43221
#
#Reference: https://s3vi.ndc.nasa.gov/ssri-kb/static/resources/Preliminary_Thermal_Analysis_of_Small_Satellites.pdf
#---------------------


import numpy as np
import matplotlib.pyplot as plt


# Variables -------
cp = 903            # sepcific heat of satellite body (aluminium) [J/(kg*K)]
m = 4               # mass of satellite [kg]
sigma = 5.67*10**(-8) # stephen-boltzman constant [W * m^-2 * K^-4]
T = 293.15          # (initial) temperature of satellite [K]
r_earth = 6873000   # radius of the earth [m]
h = 400000          # altitude of satellite [m]
tau = 5300          # revolution period [s]
alpha = 0.96        # absorbtivity of the satellite
eps = 0.90          # emissivity of satellite
summer = True       # bool, choses which q_sun to use
q_sun_hot = 1414    # specific heat rate from sun, summer dist [W * m^2]
q_sun_cold = 1322   # specific heat rate from sun, winter dist [W * m^2]
Qgen = 36           # heat generated onboard the satellite [W]

beta_crit = np.rad2deg( np.asin(r_earth/(r_earth+h)) )      # critical angle for eclipse [degrees]


# Utility --------

# albedo factor
def a(beta):        
    if beta < 30:
        return 0.14
    else:
        return 0.19

# earth IR specific heat rate [W]    
def q_IR(beta):     
    if beta < 30:
        return 228
    else:
        return 218

# fraction of orbit time spent in eclipse
def frac_E(beta):       
    if abs(beta) < beta_crit:
        # brace yourself, this is a fun one
        return 1/180 * np.rad2deg( np.acos(np.sqrt(h**2 + 2*r_earth*h)/((r_earth+h)*np.cos(np.deg2rad(beta)))) )
    else:
        return 0

# eclipsed factor (basically bool)
def s(t, beta):

    if ( tau/2 * (1-frac_E(beta)) < t%tau < tau/2 * (1+frac_E(beta)) ): 
        return 1
    else:
        return 0

# area of satellite facing sun (const for now) [m^2]
def A_sat():        
    return 0.3142 

# sun heat rate
def q_sun():
    if summer:
        return q_sun_hot
    else:
        return q_sun_cold

# Q dot funtion []
def Q(beta, t, temp):
    return q_IR(beta)* A_sat() + (1+a(beta))*q_sun()*A_sat()*s(t, beta)*alpha + Qgen - A_sat()*sigma*temp**4


# MAIN -----

#chose beta (or to iterate, simply add loop)
beta = 31

#choose time step
delt = 1    # [s]

#choose numbger of orbit cycles
n=10

#step
t = 0
Tmax = T
Tmin = T
Ts = []
while t < np.ceil(n*tau/delt):
    T += Q(beta, t, T)*delt/(cp*m)
    if T > Tmax:
        Tmax = T
    if T < Tmin:
        Tmin = T
    Ts.append(T)
    t += delt

#result
print("Tmax: "+ str(Tmax))
print("Tmin: " + str(Tmin))
print("final T: " + str(T))

#plot
plt.plot(Ts)
plt.ylabel("T[k]")
plt.xlabel("t [s]")
plt.show()











