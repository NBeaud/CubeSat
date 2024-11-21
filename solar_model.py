#---------------------
# CubeSat solar power
#
#Author: Nicolas Beaudoin
#Course: MCG43221
#
#Reference: https://publications-cnrc.canada.ca/eng/view/accepted/?id=d590e5a1-d443-4e56-bb45-a3f91d69443d
#---------------------

import numpy as np

'''
Coords:
x: cubesat coords for "width", origin on center of face
y: cubesat coords along 3U length, central skewer
z: cubesat coord for "height", origin on center of face

X,Y,Z : world coordinates, start same as x,y,z
    sun points along X

alpha: angle of satellite about x
beta: angle of satellite about z

'''



def calc_power(A,cos_phi):
    '''
    params:
    A: Area of photovoltaic cell [m^2]
    cos_phi: cosine of incident angle
    '''
    
    eta = 0.30       #efficiency of cell
    I = 1370        #incident light, ~ solar constant [W/m^2]
    
    return A*eta*I*cos_phi


#####Configurations#####
'''
panel = [A,a,b, x,y,z], where A is the area, and x,y,z is the vector to the center

where
A: area of the panel
a: angular position of panel about x axis, 0 is xz plane, positive normal along y
b: angular position of panel about z axis
x: vector component to center of panel
y: vector component to center of panel
z: vector component to center of panel

not currently using x,y,z vector
       
'''

conf1 = [[0.03, 0, 0],   #0, 0, 0.20],    #folds
         [0.03, 0, 0],   #0, 0, -0.20],   
         [0.03, 0, 0],   #0.20, 0, 0],
         [0.03, 0, 0],  #-0.20, 0, 0],
         [0.03, np.pi, 0],   #0, 0, 0.20],    #rev folds
         [0.03, np.pi, 0],   #0, 0, -0.20],   
         [0.03, 0, np.pi],   #0.20, 0, 0],
         [0.03, 0, np.pi],  #-0.20, 0, 0],
         [0.03, -np.pi/2, 0],   #0, 0.15, 0.05], #body long
         [0.03, np.pi/2, 0],   #0, 0.15, -0.05],
         [0.03, 0, -np.pi/2],   #0.05, 0.15, 0],
         [0.03, 0, np.pi/2],   #-0.05, 0.15, 0],
         [0.01, np.pi, 0],    #0, 0, 0], #caps
         [0.01, 0, 0]]    #0, 0.3, 0]]

conf2 =[[0.09, np.pi/4, np.pi/2],  #0.2, 0.15, 0],      #top up diag panel
        [0.09, -3*np.pi/4, np.pi/2],  #-0.2, 0.15, 0],       #top down
        [0.09, np.pi/4, np.pi/2],   #0.2, 0.15, 0],        #bottom up
        [0.09, -3*np.pi/4, np.pi/2],#,   -0.2, 0.15, 0]]         #bottom down
        [0.03, np.pi/2, 0],                                 #body long
        [0.03, -np.pi/2, 0],
        [0.03, 0, np.pi/2],     
        [0.03, 0, -np.pi/2],
        [0.01, np.pi, 0],     #caps
        [0.01, 0, 0]]


conf3 =[[0.06, np.pi/4, np.pi/2],  #0.2, 0.15, 0],      #top up diag panel
        [0.06, -3*np.pi/4, np.pi/2],  #-0.2, 0.15, 0],       #top down
        [0.06, np.pi/4, np.pi/2],   #0.2, 0.15, 0],        #bottom up
        [0.06, -3*np.pi/4, np.pi/2],#,   -0.2, 0.15, 0]]         #bottom down
        [0.03, np.pi/2, 0],                                 #body long
        [0.03, -np.pi/2, 0],
        [0.03, 0, np.pi/2],     
        [0.03, 0, -np.pi/2],
        [0.01, np.pi, 0],     #caps
        [0.01, 0, 0]]

conf4 = [[0.06, -np.pi/2, 0],    #body long
         [0.06, np.pi/2, 0],   
         [0.06, 0, -np.pi/2],   
         [0.06, 0, np.pi/2],   
         [0.01, 0, 0]]  #cap 

#####--------------#####


def calc_total_power(conf, alpha, beta):
    '''
    params:
    conf : configuration matrix
    alpha: alpha angle of satellite about X
    beta: beta angle of satellite about Z

    Doesn't account for shading of satellite body on itself
    and other sources of light (reflection off earth, of sat itself)

    '''
    tp = 0
    for i in conf:
        a_new = i[1]+alpha
        b_new = i[2]+beta
        cos_phi = np.dot([np.sin(a_new)*np.cos(b_new),np.sin(a_new)*np.sin(b_new) ,np.cos(a_new)],[1,0,0]) # calc incident angle

        if cos_phi<0:
            cos_phi = 0     # if value is neg, opposite direction of sun, 0 power  
        
        tp += calc_power(i[0], cos_phi)
        
    return tp

#-----------------------
# MAIN
#-----------------------

#choose conf
conf = conf4

#test a bunch of angles, 1 deg step, specific orbit remains TBD
vals = []
for alpha in range(360): 
    for beta in range(360):
        #add shaded region case, also not all angles of second dim are pertinent?
        vals.append(calc_total_power(conf,np.pi*alpha/180, np.pi*beta/180))


print("Maximal power: "+ str(max(vals)))                # max power  
print("Average* power: " + str(sum(vals)/len(vals)))      # avg power (*kinda)


