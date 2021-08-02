import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
import os.path
from os import path
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from statistics import mode

#Indicate what simulation to exploit
num_of_grains = 3000

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
			#filename_input = "flow_height_N-3000-e-0.5-mu-0.5-angle-30.txt"  #use for 1 file hdf5
			#use for all file hdf5
			filename_input = 'flow_height_N-{}-e-{}-mu-{}-angle-{}.txt'.format(num_of_grains,str(e[i_e]),str(mu[i_mu]),str(ang[i_ang]))
			if (path.exists(filename_input)==True):

					data 	 = open(filename_input,'r').readlines()

					#print(data)

					grain_size  = float(data[3])
					hstep		= float(data[5])
					porosity	= float(data[7])
					num_grains 	= float(data[9])
					tank_with	= float(data[11])
					x_position_from_gate = float(data[13])

					#filter data (attention line in txt begin at 1 while data begin at 0)
					t  = np.array(data[15][1:-2].split(),dtype=np.float64)
					dt = t[1]-t[0]

					h  = np.array(data[34][1:-2].split(','),dtype=np.float64)

					h_time_averaged      = np.array(data[39][1:-2].split(','),dtype=np.float64)
					h_time_average_std   = np.array(data[42][1:-2].split(','),dtype=np.float64)
					h_compacity_averaged = np.array(data[26][1:-2].split(','),dtype=np.float64)

					profil_compacity_100steps_averaged = data[20][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
					profil_compacity_100steps_averaged = [p.split(', ') for p in profil_compacity_100steps_averaged] #split string profil_com to list of string
					profil_compacity_100steps_averaged = [[np.float64(x) for x in y] for y in profil_compacity_100steps_averaged]# convert each string to float value

					profil_compacity_100steps_std_all = data[23][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
					profil_compacity_100steps_std_all = [p.split(', ') for p in profil_compacity_100steps_std_all] #split string profil_com to list of string
					profil_compacity_100steps_std_all = [[np.float64(x) for x in y] for y in profil_compacity_100steps_std_all]# convert each string to float value

					profil_density_instantaneous = data[31][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
					profil_density_instantaneous = [p.split(', ') for p in profil_density_instantaneous] #split string profil_com to list of string
					profil_density_instantaneous = [[np.float64(x) for x in y] for y in profil_density_instantaneous]# convert each string to float value


					#Grapsh results Time evolution of flow_heigh
					#Flow_height time evolution graphs results
					width = dt/4
					plt.figure(1)
					plt.plot(t,h,'o-',label='h_instantaneous')
					plt.errorbar(t,h_time_averaged,h_time_average_std,label='h_time_average')
					plt.plot(t,h_compacity_averaged,'o--',label='h_compacity_averaged')
					plt.legend()
					plt.suptitle('Evolution of Flow_height at x ={x_pos} dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
					plt.title(filename_input,fontsize=10)
					plt.xlabel('time t(s)')
					plt.ylabel('y/d')

					#3D Time evolution of Profil_compacity_100steps_averaged
					h_compacity_averaged_max = max(h_compacity_averaged)
					for i in range(len(profil_compacity_100steps_averaged)):
					    while(len(profil_compacity_100steps_averaged[i])<h_compacity_averaged_max):
					    	profil_compacity_100steps_averaged[i].append(0.)

					fig = plt.figure(2)
					ax = fig.add_subplot(111, projection='3d')

					z1 = np.arange(0,max(h_compacity_averaged),1)
					y1 = t
					Z1,Y1 = np.meshgrid(z1,y1)
					X1=np.zeros((len(y1),len(z1)))

					for i in range(len(y1)):
					    X1[i]=profil_compacity_100steps_averaged[i][:]


					ax.scatter3D(X1,Y1,Z1, color = "green")
					#ax.view_init(23,-100)
					plt.suptitle('Time evolution of Profil_compacity 100steps averaged at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
					plt.title(filename_input,fontsize=10)
					plt.xlabel('Compacity')
					plt.ylabel('t(s)')
					ax.set_zlabel("y/d")
					h_max = max(h)

					#3D Time evolution of density instantaneous (fraction volume)
					for i in range(len(profil_density_instantaneous)):
					    while(len(profil_density_instantaneous[i])<h_max):
					        profil_density_instantaneous[i].append(0.)

					fig = plt.figure(3)
					ax = fig.add_subplot(111, projection='3d')

					z1 = np.arange(0,h_max,1)
					y1 = t
					Z1,Y1 = np.meshgrid(z1,y1)
					X1=np.zeros((len(y1),len(z1)))

					for i in range(len(y1)):
					    X1[i]=profil_density_instantaneous[i][:]


					ax.scatter3D(X1,Y1,Z1, color = "green")
					#ax.view_init(23,-100)
					plt.suptitle('Evolution of Profil_compacity instantaneous at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
					plt.title(filename_input,fontsize=10)
					plt.xlabel('Compacity')
					plt.ylabel('t(s)')
					ax.set_zlabel("y/d")
					plt.show()