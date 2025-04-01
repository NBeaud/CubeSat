#---------------------
# CubeSat 6 node thermal model
#
# Author: Nicolas Beaudoin
# Course: MCG43221
#
# Reference: https://s3vi.ndc.nasa.gov/ssri-kb/static/resources/Preliminary_Thermal_Analysis_of_Small_Satellites.pdf
#---------------------

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt


# Environment Parameters -------
sigma = 5.67*10**-8 # stefan-boltzman constant [W * m^-2 * K^-4]
T = 293.15          # (initial) temperature of satellite [K]
r_earth = 6873000   # radius of the earth [m]
h = 1200000         # altitude of satellite [m]
tau = 6550          # revolution period [s]
summer = True       # bool, choses which q_sun to use
q_sun_hot = 1414    # specific heat rate from sun, summer dist [W * m^2]
q_sun_cold = 1322   # specific heat rate from sun, winter dist [W * m^2]

beta_crit = np.rad2deg( np.asin(r_earth/(r_earth+h)) )      # critical angle for eclipse [degrees]

# Satellite Parameters -------
k = 400             # conductance of body, conservative est. [W/(m^2 K)]
Qgen = 36           # heat generated onboard the satellite [W] 
cp = 903            # sepcific heat of satellite body (~aluminium) [J/(kg*K)]
m = 4               # mass of satellite [kg]

# index: {zenith, -v, +v, nadir, north, south}
areas = [0.03, 0.01, 0.01, 0.03, 0.03, 0.03]            # area of faces [m^2]
alphas = [0.96, 0.96, 0.96, 0.96, 0.96, 0.96]           # absorbtivity of faces 
epss = [0.90, 0.90, 0.90, 0.90, 0.90, 0.90]             # emissivity of faces


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
        return 1/np.pi * np.acos(np.sqrt(h**2 + 2*r_earth*h)/((r_earth+h)*np.cos(np.deg2rad(beta))))
    else:
        return 0

# eclipsed factor (basically bool)
def s(t, beta):

    if ( tau/2 * (1-frac_E(beta)) < t%tau < tau/2 * (1+frac_E(beta)) ): 
        return 0
    else:
        return 1

# sun heat rate
def q_sun():
    if summer:
        return q_sun_hot
    else:
        return q_sun_cold


# View Factors -------

def Fzen(t, beta):
    if tau/4 > t > 3*tau/4:
        return np.cos(np.deg2rad(2*np.pi*t/tau))*np.cos(np.deg2rad(beta)) 
    return 0

def Fvneg(t, beta):
    if t < tau/2 *(1 - frac_E(beta)):
        return np.sin(np.deg2rad(2*np.pi*t/tau))*np.cos(np.deg2rad(beta))
    return 0

def Fvpos(t, beta):
    if t > tau/2 *(1 + frac_E(beta)):
        return -np.sin(np.deg2rad(2*np.pi*t/tau))*np.cos(np.deg2rad(beta))
    return 0

def Fnad(t, beta):
    if tau/4 < t < tau/2 *(1 - frac_E(beta)) or tau/2 *(1 + frac_E(beta)) < t < 3*tau/4:
        return -np.cos(np.deg2rad(2*np.pi*t/tau))*np.cos(np.deg2rad(beta))
    return 0

def Fns(t, beta):
    if tau/2 *(1 - frac_E(beta)) > t > tau/2 *(1 + frac_E(beta)):
        return np.sin(np.deg2rad(beta))
    return 0


# Heat Transfer Q dots -------

def Qzen(Ts, t, beta):
    return ( Fzen(t, beta)*areas[0]*q_sun()*alphas[0] 
         - sigma*epss[0]*areas[0]*Ts[0]**4 
         + k*areas[1]*(Ts[1]-Ts[0])
         + k*areas[2]*(Ts[2]-Ts[0])
         + k*areas[4]*(Ts[4]-Ts[0])
         + k*areas[5]*(Ts[5]-Ts[0]) )
    
def Qvneg(Ts, t, beta):
    return ( Fvneg(t, beta)*areas[1]*q_sun()*alphas[1] 
         - sigma*epss[1]*areas[1]*Ts[1]**4 
         + k*areas[0]*(Ts[0]-Ts[1])
         + k*areas[3]*(Ts[3]-Ts[1])
         + k*areas[4]*(Ts[4]-Ts[1])
         + k*areas[5]*(Ts[5]-Ts[1]) )

def Qvpos(Ts, t, beta):
    return ( Fvpos(t, beta)*areas[2]*q_sun()*alphas[2] 
         - sigma*epss[2]*areas[2]*Ts[2]**4 
         + k*areas[0]*(Ts[0]-Ts[2])
         + k*areas[3]*(Ts[3]-Ts[2])
         + k*areas[4]*(Ts[4]-Ts[2])
         + k*areas[5]*(Ts[5]-Ts[2]) )

def Qnad(Ts, t, beta):
    return ( (Fnad(t,beta) + a(beta))*areas[3]*q_sun()*alphas[3]    #which alpha? is Azen a typo? probably
            + q_IR(beta)*areas[3]
            - sigma*epss[3]*areas[3]*Ts[3]**4
            + k*areas[1]*(Ts[1]-Ts[3])
            + k*areas[2]*(Ts[2]-Ts[3])
            + k*areas[4]*(Ts[4]-Ts[3])
            + k*areas[5]*(Ts[5]-Ts[3]) )

def Qn(Ts, t, beta):
    return ( Fns(t, beta)*areas[4]*q_sun()*alphas[4] 
         - sigma*epss[4]*areas[4]*Ts[4]**4 
         + k*areas[0]*(Ts[0]-Ts[4])
         + k*areas[1]*(Ts[1]-Ts[4])
         + k*areas[2]*(Ts[2]-Ts[4])
         + k*areas[3]*(Ts[3]-Ts[4]) )

def Qs(Ts, t, beta):
    return ( - sigma*epss[5]*areas[5]*Ts[5]**4 
         + k*areas[0]*(Ts[0]-Ts[4])
         + k*areas[1]*(Ts[1]-Ts[4])
         + k*areas[2]*(Ts[2]-Ts[4])
         + k*areas[3]*(Ts[3]-Ts[4]) )     # another Tzen typo?

def Q(Ts,t,beta):
    return ( Qzen(Ts,t,beta) 
            + Qvneg(Ts,t,beta)
            + Qvpos(Ts,t,beta)
            + Qnad(Ts,t,beta)
            + Qn(Ts,t,beta)
            + Qs(Ts,t,beta) 
            + Qgen*s(t,beta))


# MAIN -----

#choose time step
delt = 10    # [s]

#choose numbger of orbit cycles
n = 2

# datapoints for plot
betas = [0]
ts = [0]
Ts = [6*[T]]

#step
for beta in range(0,90,5):
    t = 0
    while t < n*tau:
        temp = Ts[-1].copy()
        temp[0] += Qzen(temp, t, beta)*delt/(cp*m)   # m as fraction with area weight     
        temp[1] += Qvneg(temp, t, beta)*delt/(cp*m)
        temp[2] += Qvpos(temp, t, beta)*delt/(cp*m)
        temp[3] += Qnad(temp, t, beta)*delt/(cp*m)
        temp[4] += Qn(temp, t, beta)*delt/(cp*m)
        temp[5] += Qs(temp, t, beta)*delt/(cp*m)
        Ts.append(temp)
        ts.append(t)
        t += delt
        betas.append(beta)

print(Ts[1])


Tavg = []
for i in range(len(Ts)):
    Tavg.append(sum(Ts[i])/len(Ts[i]))

#results
print("Tmax: "+ str(max(Tavg)))
print("Tmin: " + str(min(Tavg)))
print("final T: " + str(T))

# 2D plot
plt.plot(ts,Tavg,".")
plt.ylabel("T [K]")
plt.xlabel("t [s]")
plt.show()

# 3D plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_trisurf(ts,betas,Tavg, alpha = 0.1)
ax.scatter(ts,betas,Tavg, )
ax.set_xlabel("t [s]")
ax.set_ylabel("$\\beta$ [degrees]")
ax.set_zlabel("T [K]")
ax.set_title("Temperature Evolution According to Orbit Beta Angle")
plt.show()