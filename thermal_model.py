#---------------------
# CubeSat thermal model
#
# Author: Nicolas Beaudoin
# Course: MCG43221
#
# Reference: https://s3vi.ndc.nasa.gov/ssri-kb/static/resources/Preliminary_Thermal_Analysis_of_Small_Satellites.pdf
#---------------------

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


# Parameters -------
cp = 903            # sepcific heat of satellite body (~aluminium) [J/(kg*K)]
m = 4               # mass of satellite [kg]
sigma = 5.67*10**-8 # stefan-boltzman constant [W * m^-2 * K^-4]
T = 273.15          # (initial) temperature of satellite [K]
r_earth = 6873000   # radius of the earth [m]
h = 1200000         # altitude of satellite [m]
tau = 6550          # revolution period [s]
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
        return 1/np.pi * np.acos(np.sqrt(h**2 + 2*r_earth*h)/((r_earth+h)*np.cos(np.deg2rad(beta))))
    else:
        return 0

# eclipsed factor (basically bool)
def s(t, beta):

    if ( tau/2 * (1-frac_E(beta)) < t%tau < tau/2 * (1+frac_E(beta)) ): 
        return 0
    else:
        return 1

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
    return q_IR(beta)* A_sat() + (1+a(beta))*q_sun()*A_sat()*s(t, beta)*alpha + s(t, beta)*Qgen - A_sat()*sigma*temp**4



# MAIN -----

#choose time step
delt = 10    # [s]

#choose numbger of orbit cycles
n=2

# datapoints for plot
betas = []
ts = []
Ts = []

#step
for beta in range(0,90,5):
    t = 0
    while t < n*tau:
        T += Q(beta, t, T)*delt/(cp*m)
        Ts.append(T)
        ts.append(t)
        t += delt
        betas.append(beta)

#results
print("Tmax: "+ str(max(Ts)))
print("Tmin: " + str(min(Ts)))
print("final T: " + str(T))

# 2D plot
plt.plot(Ts)
plt.ylabel("T [K]")
plt.xlabel("t [s]")
plt.show()

# 3D plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_trisurf(ts,betas,Ts, alpha = 0.1)
ax.scatter(ts,betas,Ts, )
ax.set_xlabel("t [s]")
ax.set_ylabel("$\\beta$ [degrees]")
ax.set_zlabel("T [K]")
plt.show()












