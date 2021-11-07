#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 18:05:08 2021

@author: arefe
"""


import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import matplotlib.animation
from matplotlib.animation import FuncAnimation
import ffmpeg

def polar2euc (r , theta ):
    x = r* np.cos(theta)
    y = r* np.sin(theta)

    return x , y




#constants
N_init = 10000
R_max = 16
R_min = 1
G = 6.67e-11 * ((3e16)**-3) * (2e30)
Const = G*1e6
#Const = 1
T = 3*10**13
delta_t = 10**11
num_per = int(0.05*N_init)
num_spir = 4
N = N_init+num_per*num_spir

close_abgle_interval = 10*np.pi/180
Cx = 0.1*Const
Cy = 0.1*Const
num_avg = 5 

r_init = R_min+ np.random.rand(N_init)*R_max
theta_init = np.random.rand(N_init)*2*np.pi
x_init  , y_init= polar2euc (r_init , theta_init)

#plt.figure(figsize=(5,5))
#plt.scatter(x_init,y_init)

r = r_init
theta = theta_init

for i in range (num_spir ):
    
    theta_per = (0 + i * int (360/num_spir))*np.pi/180
    theta_add = theta_per + np.random.rand(num_per)*0.04
    theta = np.append (theta , theta_add)



r_add =  R_min+ np.random.rand(num_per*num_spir)*R_max
r = np.append (r  ,r_add)


x , y = polar2euc (r , theta)

#plt.figure(figsize=(5,5))
#plt.scatter(x,y ,  0.5)

def velocity  (C , r):
    return (np.sqrt (C/r))

def velocity_close (G , r , theta , theta_spir):
    
    return (np.sqrt(G/(r*(np.abs(theta-theta_spir)))))
    

def dist (x1,y1,x2, y2):
    return np.sqrt((x2-x1)**2+(y2-y1)**2)


t_arr = np.arange (T/delta_t)
vel = velocity (Const , r)
vel_close = 0.001
delta_theta = delta_t * vel/r 
delta_theta_prime = delta_t*vel_close/r
t_steps = int (T/delta_t)
theta_arr = np.zeros ((N,t_steps)) 
x_arr = np.zeros ((N,t_steps)) 
y_arr = np.zeros ((N,t_steps)) 

#theta_arr [0 , :] = theta    
for i in range (t_steps):
    #print (theta_arr.shape)
    theta_arr[: , i ] = theta + i * delta_theta

    x_arr [: ,i ] , y_arr[: , i] = polar2euc (r, theta_arr[: , i])
    #for n in range (N):
        #x_close_arr = np.abs (x_arr[: , i] - x_arr[n,i])
        #y_close_arr = np.abs (y_arr[: , i] - y_arr[n,i])
        #x_close_avg = np.mean(np.sort (x_close_arr)[0:num_avg])
        #y_close_avg = np.mean(np.sort (y_close_arr)[0:num_avg])
        
        #dist =  dist(x_arr[n,i] , y_arr[n,i] , x_arr[: , i] , y_arr[: , i])
        #closest_neighbors = np.argsort (dist)[1:num_avg]
        #print (closest_neighbors)
      
       
        #x_close_avg = np.mean(x_arr[closest_neighbors , i])
        #y_close_avg = np.mean(y_arr[closest_neighbors , i])
        #x_arr[n ,i] += np.sign (x_arr[n ,i ] - x_close_avg) * vel_close
        #y_arr[n ,i] += np.sign (y_arr[n ,i ] - y_close_avg) * vel_close
        #print (x_close_avg , y_close_avg)
        
        

def animate(i):
    
    theta_f = theta_arr [: , i  ]
    #print (theta_f)
    x_f = x_arr[: , i]
    y_f = y_arr[: , i]

    
    ax.clear()
    ax.scatter(x_f, y_f , 0.1)
    #ax.set_xlim([0,20])
    #ax.set_ylim([0,10])
    
fig, ax = plt.subplots()    
ani = FuncAnimation(fig, animate, frames=t_steps, interval=1, repeat=True , save_count=200)

plt.show()


Writer = matplotlib.animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
ani.save('im.mp4', writer=writer)
#f = r"c://Users/Arefe/Desktop/animation.mp4" 
#writervideo = matplotlib.animation.FFMpegWriter(fps=30) 
#ani.save(f, writer=writer)
#ani.save ("ani.fps")

