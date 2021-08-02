import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import os.path
from os import path
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from statistics import mode
"""
This script do:
-open one by one Profil_Vx.txt file
-Draw graphs resutls :
+3D time evolution of Profil_Vx intantaneous
+3D time evolution of Profil_Vx time averaged
-Time evolution of Flow rate indtantaneous and time-averaged

"""
#Indicate what simulation to exploit
num_of_grains = 5000


#simulation attibute:
grain_size       = 0.001
ground_thickness = 10 * grain_size
ground_y_shift   = 0.09
y0               = - ground_thickness - ground_y_shift

hstep = 1e-4

# ordre : e-mu-angle											
e_init = 0.0
e_end  = 0.5
d_e    = 0.25
e = np.arange(e_init,e_end+d_e,d_e) 

mu_init = 0.3
mu_end  = 0.7
d_mu 	= 0.2
mu = np.arange(mu_init,mu_end+d_mu,d_mu)

ang_init = 16
ang_end  = 32
d_ang 	 = 2
ang = np.arange(ang_init,ang_end+d_ang,d_ang)


print('e = ', e)
print('mu = ',mu)
print('ang = ',ang)


# we check if exist and open all files .txt one by one and take out graphs results
for i_e in range(len(e)):
	for i_mu in range(len(mu)):
		for i_ang in range(len(ang)):
			#filename_input_input = "flow_height_N-3000-e-0.5-mu-0.5-angle-30.txt"  #use for 1 file hdf5
			#use for all file hdf5
			filename_input = 'Profil_Vxy_N-{}-e-{}-mu-{}-angle-{}.txt'.format(num_of_grains,str(e[i_e]),str(mu[i_mu]),str(ang[i_ang]))
			if (path.exists(filename_input)==True):
				#open Profil_Vx.txt
				f    = open(filename_input,'r')
				data =f.readlines()
				
				#filter data (attention line in txt begin at 1 while data begin at 0)
				t  = np.array(data[9][1:-2].split(', '),dtype=np.float64)
				dt = t[1]-t[0]

				x_position_from_gate = float(data[7])
			
				h1_time_averaged_ratio      = np.array(data[17][1:-2].split(', '),dtype=np.float64)
				h1_ratio = np.array(data[19][1:-2].split(', '),dtype=np.float64)

				Vx_profil_100steps_averaged = data[24][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				Vx_profil_100steps_averaged = [p.split(', ') for p in Vx_profil_100steps_averaged] #split string profil_com to list of string
				Vx_profil_100steps_averaged = [[np.float64(x) for x in y] for y in Vx_profil_100steps_averaged]# convert each string to float value
				
				Vx_profil_100steps_averaged_std = data[27][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				Vx_profil_100steps_averaged_std = [p.split(', ') for p in Vx_profil_100steps_averaged_std] #split string profil_com to list of string
				Vx_profil_100steps_averaged_std = [[np.float64(x) for x in y] for y in Vx_profil_100steps_averaged_std]# convert each string to float value
				
				Vx_profil_instantaneous = data[36][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				Vx_profil_instantaneous = [p.split(', ') for p in Vx_profil_instantaneous] #split string profil_com to list of string
				Vx_profil_instantaneous = [[np.float64(x) for x in y] for y in Vx_profil_instantaneous]# convert each string to float value
				
				Vy_profil_100steps_averaged = data[47][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				Vy_profil_100steps_averaged = [p.split(', ') for p in Vy_profil_100steps_averaged] #split string profil_com to list of string
				Vy_profil_100steps_averaged = [[np.float64(x) for x in y] for y in Vy_profil_100steps_averaged]# convert each string to float value
				
				Vy_profil_100steps_averaged_std = data[50][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				Vy_profil_100steps_averaged_std = [p.split(', ') for p in Vy_profil_100steps_averaged_std] #split string profil_com to list of string
				Vy_profil_100steps_averaged_std = [[np.float64(x) for x in y] for y in Vy_profil_100steps_averaged_std]# convert each string to float value
				
				Vy_profil_instantaneous = data[53][2:-2].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				Vy_profil_instantaneous = [p.split(', ') for p in Vy_profil_instantaneous] #split string profil_com to list of string
				Vy_profil_instantaneous = [[float(x) for x in y] for y in Vy_profil_instantaneous]# convert each string to float value


				flow_rate_time_averaged = np.array(data[30][1:-2].split(','),dtype=np.float64)
				flow_rate_instantaneous = np.array(data[39][1:-2].split(','),dtype=np.float64)

				Vx_integration_h_time_averaged = np.array(data[33][1:-2].split(','),dtype=np.float64)
				Vx_integration_h_instantaneous = np.array(data[42][1:-2].split(','),dtype=np.float64)
				
				#3D time evolution of profil_Vx_time averaged
     
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')

				z1    = h1_time_averaged_ratio
				y1    = t
				Z1,Y1 = np.meshgrid(z1,y1)
				X1    = np.zeros((len(y1),len(z1)))

				for i in range(len(y1)):
					X1[i]=np.array(Vx_profil_100steps_averaged)[i,:]

				ax.scatter3D(X1,Y1,Z1, color = "green")
				#ax.view_init(23,-100)
				plt.suptitle(' Profil_Vx time-averaged versus time at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
				plt.title(filename_input, fontsize=10)
				plt.xlabel('Vx_time-averaged (m/s)')
				plt.ylabel('t(s)')
				ax.set_zlabel("y/d")



				#3D time evolution of profil_Vx instantaneous
				fig = plt.figure()
				ax  = fig.add_subplot(111, projection='3d')

				z2    = h1_ratio
				y2    = t
				Z2,Y2 = np.meshgrid(z2,y2)
				X2    =np.zeros((len(y2),len(z2)))

				#f_out.write(len(y2))
				for i in range(len(y2)):
					X2[i]=np.array(Vx_profil_instantaneous)[i,:]

				ax.scatter3D(X2,Y2,Z2, color = "green")
				#ax.set_xlim(0, 23)
				#ax.set_zlim(0, 2)
				#ax.view_init(1,-100)
				plt.suptitle('Profil_Vx_instantaneou versus time at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
				plt.title(filename_input, fontsize=10)
				plt.xlabel('Vx (m/s)')
				plt.ylabel('t(s)')
				ax.set_zlabel("y/d")

				#3D time evolution of profil_Vy_time averaged
     
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')

				z1    = h1_time_averaged_ratio
				y1    = t
				Z1,Y1 = np.meshgrid(z1,y1)
				X1    = np.zeros((len(y1),len(z1)))

				for i in range(len(y1)):
					X1[i]=np.array(Vy_profil_100steps_averaged)[i,:]

				ax.scatter3D(X1,Y1,Z1, color = "green")
				#ax.view_init(23,-100)
				plt.suptitle(' Profil_Vy time-averaged versus time at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
				plt.title(filename_input, fontsize=10)
				plt.xlabel('Vy_time-averaged (m/s)')
				plt.ylabel('t(s)')
				ax.set_zlabel("y/d")

				#3D time evolution of profil_Vy instantaneous
				fig = plt.figure()
				ax  = fig.add_subplot(111, projection='3d')

				z2    = h1_ratio
				y2    = t
				Z2,Y2 = np.meshgrid(z2,y2)
				X2    =np.zeros((len(y2),len(z2)))

				#f_out.write(len(y2))
				for i in range(len(y2)):
					X2[i]=np.array(Vy_profil_instantaneous)[i,:]

				ax.scatter3D(X2,Y2,Z2, color = "green")
				#ax.set_xlim(0, 23)
				#ax.set_zlim(0, 2)
				#ax.view_init(1,-100)
				plt.suptitle('Profil_Vy_instantaneou versus time at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
				plt.title(filename_input, fontsize=10)
				plt.xlabel('Vy (m/s)')
				plt.ylabel('t(s)')
				ax.set_zlabel("y/d")


				#graph results of time evolution of flow_rate
				width = dt/4
				plt.figure()
				plt.plot(t,flow_rate_time_averaged,'-o',label='flow_rate_time-average')
				plt.plot(t,flow_rate_instantaneous,label='flow_rate')
				plt.legend()
				plt.suptitle('Flow_rate versus time at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
				plt.title(filename_input, fontsize=10)
				plt.xlabel('time t(s)')
				plt.ylabel('m3.s-1')

				#graph results of time evolution of integral of Vx on flow_height
				width = dt/4
				plt.figure()
				plt.plot(t,Vx_integration_h_time_averaged,'-o',label='integral of Vx_time averaged')
				plt.plot(t,Vx_integration_h_instantaneous,label='integral of Vx_time instantaneous')
				plt.legend()
				plt.suptitle('integral of Vx versus time at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
				plt.title(filename_input, fontsize=10)
				plt.xlabel('time t(s)')
				plt.ylabel('s-1')
				


				plt.show();
				f.close()