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
This script will open one by one:
-Flow_height.txt and draw comparaison graphs of flow_height and profil_compacity
-profil_Vx .txt and draw comparaison graphs of profil Vx and flow_rate

"""

#simulation attibute:
grain_size       = 0.001
ground_thickness = 10 * grain_size
ground_y_shift   = 0.09
y0               = - ground_thickness - ground_y_shift

hstep = 1e-4
"""
num_of_grains = input('Please enter num_of_grains  = ')
angle    = input('Please enter inclination [16°,32°](step=2°)  = ')
ref      = input('Please choose the fixed attribute to take out graphs results (selon e or mu ?) : ')
t_graphs = np.float64(input('Please enter time to profil_compacity and profil_Vx graphs ([0.,3.],dt=0.5): t = '))
"""
num_of_grains = 30000
angle    = 20
ref      = 'e'
t_graphs = 3
coef 	 = 1 #V_thresold = coef * V_X

# To choose position x_wall from gate to the calculation
x_position_from_gate = 700

# ordre : e-mu-angle											
e_init = 0.0
e_end  = 0.5
d_e    = 0.25
e = np.arange(e_init,e_end+d_e,d_e) 

mu_init = 0.3
mu_end  = 0.7
d_mu 	= 0.2
mu = np.arange(mu_init,mu_end+d_mu,d_mu)


print('e = ', e)
print('mu = ',mu)
print('angle = ',angle)


if (ref=='e'):

	L_dead_zone_fig, L_dead_zone_graphs  	 = plt.subplots(2,len(e),figsize=(14, 6), sharey=True)
	S_dead_zone_fig, S_dead_zone_graphs  	 = plt.subplots(2,len(e),figsize=(14, 6), sharey=True)
	num_fig, num_graphs  	 			     = plt.subplots(2,len(e),figsize=(14, 6), sharey=True)
	com_dead_zone_fig, com_dead_zone_graphs  = plt.subplots(2,len(e),figsize=(14, 6), sharey=True)
	
	force_fig_inst, force_graphs_inst   	 = plt.subplots(3,len(e),figsize=(14, 6), sharey=True)
	force_fig_time_ave,force_graphs_time_ave = plt.subplots(3,len(e),figsize=(14, 6), sharey=True)
	'''	
	wbound = -0.1
	ebound = 3.1
	sbound = -0.1
	nbound = 150.1
	#plt.xlim(wbound, ebound)
	#plt.ylim(sbound, nbound)
	'''
	# we check if exist and open all files .txt one by one and take out graphs results
	for i_e in range(len(e)):
		for i_mu in range(len(mu)):
			
			"""
			open Dead_zone.txt and take out graphs:
			- compare L_dead_zone by Vx and Vnorme
			- compare S_dead_zone by Vx and Vnorme
			"""
			dead_zone_input = 'Dead_zone_wall_{}_N-{}-e-{}-mu-{}-angle-{}.txt'.format(x_position_from_gate,num_of_grains,e[i_e],mu[i_mu],angle)
			if (path.exists(dead_zone_input)==True):
				f_d 	   	= open(dead_zone_input,'r')
				data_d 		= f_d.readlines()

				num_grains 	= float(data_d[3])
				tank_with	= float(data_d[5])
				#x_position_from_gate = float(data_d[7])
				t  = np.array(data_d[9][1:-2].split(', '),dtype=np.float64)
				dt = t[1]-t[0]

				#we take results at t = 3.0 (line 1119)
				row_ecart  = coef*(40-11)

				index_L_vx    = int(18+row_ecart)
				index_S_vx	  = int(20+row_ecart)
				index_num_vx  = int(22+row_ecart)
				index_frac_vx = int(24+row_ecart)

				index_L_vn    = int(32+row_ecart)
				index_S_vn	  = int(34+row_ecart)
				index_num_vn  = int(36+row_ecart)
				index_frac_vn = int(38+row_ecart)

				L_dead_zone_all_steps_Vx = np.array(data_d[index_L_vx][1:-2].split(','),dtype=np.float64)
				S_dead_zone_all_steps_Vx = np.array(data_d[index_S_vx][1:-2].split(','),dtype=np.float64)
				num_all_steps_Vx 	 	 = np.array(data_d[index_num_vx][1:-2].split(','),dtype=np.float64)
				frac_all_steps_Vx 	 	 = np.array(data_d[index_frac_vx][1:-2].split(','),dtype=np.float64)
				

				L_dead_zone_all_steps_Vn = np.array(data_d[index_L_vn][1:-2].split(','),dtype=np.float64)
				S_dead_zone_all_steps_Vn = np.array(data_d[index_S_vn][1:-2].split(','),dtype=np.float64)
				num_all_steps_Vn 	 	 = np.array(data_d[index_num_vn][1:-2].split(','),dtype=np.float64)
				frac_all_steps_Vn 	 	 = np.array(data_d[index_frac_vn][1:-2].split(','),dtype=np.float64)
				
				#Graphs L_dead_zone Vx and V_norme
				L_dead_zone_fig.suptitle('L_dead_zone versus time for Vx and V_norme at x ={x_pos}, coef_Vthresold = {co} , dt = {dt}, ang={a}'.format(co = coef,x_pos= x_position_from_gate,dt=dt,a=angle),fontsize=15)
				L_dead_zone_graphs[0,i_e].plot(t,L_dead_zone_all_steps_Vx,'o-',label='L_Vx-mu-{}'.format(mu[i_mu]))
				L_dead_zone_graphs[0,i_e].legend()
				L_dead_zone_graphs[0,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				L_dead_zone_graphs[0,i_e].set_xlabel('time t(s)')
				L_dead_zone_graphs[0,i_e].set_ylabel('x/d')
				
				L_dead_zone_graphs[1,i_e].plot(t,L_dead_zone_all_steps_Vn,'o-',label='L_Vnorme-mu-{}'.format(mu[i_mu]))
				L_dead_zone_graphs[1,i_e].legend()
				L_dead_zone_graphs[1,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				L_dead_zone_graphs[1,i_e].set_xlabel('time t(s)')
				L_dead_zone_graphs[1,i_e].set_ylabel('x/d')

				#Graphs S_dead_zone Vx and V_norme
				S_dead_zone_fig.suptitle('S_dead_zone versus time for Vx and V_norme at x ={x_pos}, coef_Vthresold={co} , dt = {dt}, ang={a}'.format(co= coef,x_pos= x_position_from_gate,dt=dt,a=angle),fontsize=15)
				S_dead_zone_graphs[0,i_e].plot(t,S_dead_zone_all_steps_Vx,'o-',label='L_Vx-mu-{}'.format(mu[i_mu]))
				S_dead_zone_graphs[0,i_e].legend()
				S_dead_zone_graphs[0,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				S_dead_zone_graphs[0,i_e].set_xlabel('time t(s)')
				S_dead_zone_graphs[0,i_e].set_ylabel('x/d')
				
				S_dead_zone_graphs[1,i_e].plot(t,S_dead_zone_all_steps_Vn,'o-',label='L_Vnorme-mu-{}'.format(mu[i_mu]))
				S_dead_zone_graphs[1,i_e].legend()
				S_dead_zone_graphs[1,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				S_dead_zone_graphs[1,i_e].set_xlabel('time t(s)')
				S_dead_zone_graphs[1,i_e].set_ylabel('x/d')

				#graphs num_grains in dead_zone with Vx and V_norme
				num_fig.suptitle('Quantity of grains in dead_zone versus time for Vx and V_norme at x ={x_pos}, coef_Vthresold = {co} , dt = {dt}, ang={a}'.format(co = coef,x_pos= x_position_from_gate,dt=dt,a=angle),fontsize=15)
				num_graphs[0,i_e].plot(t,num_all_steps_Vx,'o-',label='num-grains_Vx-mu-{}'.format(mu[i_mu]))
				num_graphs[0,i_e].legend()
				num_graphs[0,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				num_graphs[0,i_e].set_xlabel('time t(s)')
				num_graphs[0,i_e].set_ylabel('grains')
				
				num_graphs[1,i_e].plot(t,num_all_steps_Vn,'o-',label='num-grains_Vnorme-mu-{}'.format(mu[i_mu]))
				num_graphs[1,i_e].legend()
				num_graphs[1,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				num_graphs[1,i_e].set_xlabel('time t(s)')
				num_graphs[1,i_e].set_ylabel('grains')

				#graphs compacity of dead_zone with Vx and V_norme
				com_dead_zone_fig.suptitle('Compacity of dead_zone versus time for Vx and V_norme at x ={x_pos}, coef_Vthresold = {co} , dt = {dt}, ang={a}'.format(co = coef,x_pos= x_position_from_gate,dt=dt,a=angle),fontsize=15)
				com_dead_zone_graphs[0,i_e].plot(t,frac_all_steps_Vx,'o-',label='com_Vx-mu-{}'.format(mu[i_mu]))
				com_dead_zone_graphs[0,i_e].legend()
				com_dead_zone_graphs[0,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				com_dead_zone_graphs[0,i_e].set_xlabel('time t(s)')
				com_dead_zone_graphs[0,i_e].set_ylabel('fraction_volume')
				
				com_dead_zone_graphs[1,i_e].plot(t,frac_all_steps_Vn,'o-',label='com_Vnorme-mu-{}'.format(mu[i_mu]))
				com_dead_zone_graphs[1,i_e].legend()
				com_dead_zone_graphs[1,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				com_dead_zone_graphs[1,i_e].set_xlabel('time t(s)')
				com_dead_zone_graphs[1,i_e].set_ylabel('fraction_volume')


			else:
 				print('file',dead_zone_input,'does not exists')


			"""
			open Force_on_wall.txt and take out graphs:
			- compare normal force , tangantial force and norme force
			"""
			force_wall_input = 'force_on_wall_{}d_N-{}-e-{}-mu-{}-angle-{}.txt'.format(x_position_from_gate,num_of_grains,e[i_e],mu[i_mu],angle)
			if (path.exists(force_wall_input)==True):
				f_w 	   	 = open(force_wall_input,'r')
				data_w 		 = f_w.readlines()

				normal_force_inst     = np.array(data_w[11][1:-2].split(', '),dtype=np.float64)
				tangential_force_inst = np.array(data_w[13][1:-2].split(', '),dtype=np.float64)
				norm_force_inst       = np.array(data_w[15][1:-2].split(', '),dtype=np.float64)

				normal_force_ta     = np.array(data_w[20][1:-2].split(', '),dtype=np.float64)
				tangential_force_ta = np.array(data_w[24][1:-2].split(', '),dtype=np.float64)
				norm_force_ta       = np.array(data_w[28][1:-2].split(', '),dtype=np.float64)
				
				normal_force_ta_std     = np.array(data_w[22][1:-2].split(', '),dtype=np.float64)
				tangential_force_ta_std = np.array(data_w[26][1:-2].split(', '),dtype=np.float64)
				norm_force_ta_std       = np.array(data_w[30][1:-2].split(', '),dtype=np.float64)
				
				#Graphs forc_on_wall
				# instantaneous
				force_fig_inst.suptitle('Force (normal,tangential and norm) versus time at x ={x_pos}, dt = {d_t}, ang={a}'.format(co = coef,x_pos= x_position_from_gate,d_t=dt,a=angle),fontsize=15)
				force_graphs_inst[0,i_e].plot(t,normal_force_inst,'o-',label='F_normal-mu-{}'.format(mu[i_mu]))
				force_graphs_inst[0,i_e].legend()
				force_graphs_inst[0,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				force_graphs_inst[0,i_e].set_xlabel('time t(s)')
				force_graphs_inst[0,i_e].set_ylabel('F (N.m-1)')
				
				force_graphs_inst[1,i_e].plot(t,tangential_force_inst,'o-',label='F_tangential-mu-{}'.format(mu[i_mu]))
				force_graphs_inst[1,i_e].legend()
				force_graphs_inst[1,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				force_graphs_inst[1,i_e].set_xlabel('time t(s)')
				force_graphs_inst[1,i_e].set_ylabel('F (N.m-1)')

				force_graphs_inst[2,i_e].plot(t,norm_force_inst,'o-',label='F_norm-mu-{}'.format(mu[i_mu]))
				force_graphs_inst[2,i_e].legend()
				force_graphs_inst[2,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				force_graphs_inst[2,i_e].set_xlabel('time t(s)')
				force_graphs_inst[2,i_e].set_ylabel('F (N.m-1)')
				
				#time_averaged
				force_fig_time_ave.suptitle('Force (normal,tangential and norm) time_averaged versus time at x ={x_pos}, dt = {dt}, ang={a}'.format(co = coef,x_pos= x_position_from_gate,dt=dt,a=angle),fontsize=15)
				force_graphs_time_ave[0,i_e].errorbar(t,normal_force_ta,normal_force_ta_std,label='F_normal-mu-{}'.format(mu[i_mu]))
				force_graphs_time_ave[0,i_e].legend()
				force_graphs_time_ave[0,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				force_graphs_time_ave[0,i_e].set_xlabel('time t(s)')
				force_graphs_time_ave[0,i_e].set_ylabel('F (N.m-1)')
				
				force_graphs_time_ave[1,i_e].errorbar(t,tangential_force_ta,tangential_force_ta_std,label='F_tangential-mu-{}'.format(mu[i_mu]))
				force_graphs_time_ave[1,i_e].legend()
				force_graphs_time_ave[1,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				force_graphs_time_ave[1,i_e].set_xlabel('time t(s)')
				force_graphs_time_ave[1,i_e].set_ylabel('F (N.m-1)')

				force_graphs_time_ave[2,i_e].errorbar(t,norm_force_ta,norm_force_ta_std,label='F_norm-mu-{}'.format(mu[i_mu]))
				force_graphs_time_ave[2,i_e].legend()
				force_graphs_time_ave[2,i_e].set_title('e = '+str(e[i_e]),fontsize=10)
				force_graphs_time_ave[2,i_e].set_xlabel('time t(s)')
				force_graphs_time_ave[2,i_e].set_ylabel('F (N.m-1)')

plt.show()