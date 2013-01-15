# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 13:39:14 2012

@Case: Simple 2D Cylinder
@Methodology: Solving the Cylinder using Panel Method
@author: lento

"""
############################################################### 
'''            Panel Method: Simple Cylinder                '''
###############################################################


############################################################### 
'''                 Initialization                          '''
# Pollute the namespace but also provide MATLAB-like experience:
from pylab import *  #analysis:ignore
# Enable Matplotlib's interactive mode:
ion()
from numpy import *
from scipy import *
from pylab import *

#importing custom modules
from vort2D import*
from unitVec import *

import streamlines as sl
  
############################################################### 
'''                 Control Parameters                      '''

# Flow Conditions
Uinf    = 10.   # Free stream velocity : x-direction
Winf    = 0.    # Free stream velocity : y-direction

# Basic Geometry
R       = 1.    # radius of cylinder
n       = 100   # no. of  points
m       = n     # no. of panels (same)

meshDensity = 100  # Mesh [for flow field]

###############################################################
'''                Panel + Control Points                   '''
    
# Generating Panel joints 
# Cylindrical Coordinates:
r       = repeat(R,n+1)  # radius for n+1 point of n panels
theta   = linspace(pi,-pi,(n+1))
# In cartesian coordinates
xp      = r*cos(theta) # Panel Joints: x-coordinates
yp      = r*sin(theta) # Panel Joints: y-coordinates
# Grouping the coorinates in 
panelCoor = concatenate(([xp],[yp]))

# Genrating Control points
xcp     = (xp[1:]+xp[:-1])/2 # x-coordinates
ycp     = (yp[1:]+yp[:-1])/2 # y-coordinates
# Grouping the coorinates in 
controlCoor = zeros([2,n])

norm = normVec(panelCoor[:,:-1],panelCoor[:,1:])
tang = tangVec(panelCoor[:,:-1],panelCoor[:,1:])

controlCoor[0,:] = xcp[:] + 100*norm[0]*finfo(float).eps # pushing control point outwards
controlCoor[1,:] = ycp[:] + 100*norm[1]*finfo(float).eps # to erase numerical inaccuracy       
    
####################################################################
'''                 Solving the problem                         '''

# Initial Variables 

# Generating meshgrid of points to evalute (shape n[points] by m[panels])
sigma = ones([n,m])

controlPointsx = tile(array([controlCoor[0]]).transpose(),[1,m])
controlPointsy = tile(array([controlCoor[1]]).transpose(),[1,m])

panelOriginx = tile(panelCoor[0,:-1], [n,1])
panelOriginy = tile(panelCoor[1,:-1], [n,1])

panelEndx = tile(panelCoor[0,1:], [n,1])
panelEndy = tile(panelCoor[1,1:], [n,1])

normx = tile(array([norm[0]]).transpose(),[1,m])
normy = tile(array([norm[1]]).transpose(),[1,m])

# Calculating the A matrix
u,w     = sor2D(sigma, (controlPointsx, controlPointsy), (panelOriginx, panelOriginy), (panelEndx, panelEndy))
A       = u*normx + w*normy
RHS     = -(Uinf*norm[0] + Winf*norm[1])
# Solving for the Sigma Matrix
Sigma = solve(A,RHS)

####################################################################
'''               Calculating the results                        '''

#Induced Velocity at control point
u,w = sor2D(Sigma, (controlPointsx, controlPointsy), (panelOriginx, panelOriginy), (panelEndx, panelEndy))
Q     = concatenate(([sum(u,axis=1)],[sum(w,axis=1)])) + array([[Uinf],[Winf]]) 
Qt    = Q[0,:]*tang[0] + Q[1,:]*tang[1] # tangential velocity
Qn    = Q[0,:]*norm[0] + Q[1,:]*norm[1] # normal velocity
Qres  = sqrt(pow(Q[0,:],2) + pow(Q[1,:],2)) # resultant velocity

    
# Calculating the flow field 
meshx,meshy = meshgrid(linspace(-2*R,2*R,meshDensity),linspace(-2*R,2*R,meshDensity))

# Solving the problem using 3D array (no. panel,meshx,meshy)
# Generating points
controlPointsx = tile(meshx,[m,1,1])
controlPointsy = tile(meshy,[m,1,1])

panelOriginx = tile(reshape(panelCoor[0,:-1],(m,1,1)),[1,meshDensity,meshDensity])
panelOriginy = tile(reshape(panelCoor[1,:-1],(m,1,1)),[1,meshDensity,meshDensity])

panelEndx = tile(reshape(panelCoor[0,1:],(m,1,1)),[1,meshDensity,meshDensity])
panelEndy = tile(reshape(panelCoor[1,1:],(m,1,1)),[1,meshDensity,meshDensity])

sourceSigma = tile(reshape(Sigma,(m,1,1)),[1,meshDensity,meshDensity])

u,w = sor2D(sourceSigma, (controlPointsx, controlPointsy), (panelOriginx, panelOriginy), (panelEndx, panelEndy))

meshVelx = sum(u,axis=0) + Uinf
meshVely = sum(w,axis=0) + Winf
meshVelres = sqrt(pow(meshVelx,2) + pow(meshVely,2))

Cp = 1-(Qt**2)/(Uinf**2) # Calculating the pressure coefficients

streamlns  = sl.Streamlines(meshx,meshy,meshVelx,meshVely)
     
####################################################################
'''                     Post-Processing                         '''

#Figure 1
figure(1)
title('Single 2D Cylinder [%d panels]'%(n))
plot(panelCoor[0,:],panelCoor[1,:],'b')
plot(controlCoor[0,:],controlCoor[1,:],'k.',lw=5)
quiver(controlCoor[0,:],controlCoor[1,:],Q[0,:],Q[1,:],Qres)
quiver(meshx,meshy,meshVelx,meshVely,meshVelres)
colorbar()
xlabel('x-coordinate [-]')
ylabel('y-coordinate [-]')
axis([-2, 2, -2, 2])
axis('equal')
grid()

#Figure 2
figure(2)
title('Pressure Coefficient')
plot(controlCoor[0,:],Cp)
xlabel('x-coordinate [-]')
ylabel('Pressure Coefficient $C_{p}$ [-]')
axis([-2, 2, 2, -2])
axis('equal')
grid()

figure(3)
title('Flow Field')
contourf(meshx,meshy,meshVelres)
plot(panelCoor[0,:],panelCoor[1,:],'k')
xlabel('x-coordinate [-]')
ylabel('y-coordinate [-]')
axis([-2,2,-2,2])
axis('equal')
colorbar()
streamlns.plot()