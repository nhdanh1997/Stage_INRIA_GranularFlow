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
This script will do :
- open one by one dead_zone.txt
- take out results and draw graphs: 
+ evolution of dead_zone
+ evolution of S_dead_zone and L_dead_zone.
"""
#Indicate what simulation to exploit
num_of_grains = 30000

#enter wall x_position_from_gate
x_position_from_gate = 700

#enter coefficient de thresold
coef = 1

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

ang_init = 12
ang_end  = 40
d_ang 	 = 2
ang = np.arange(ang_init,ang_end+d_ang,d_ang)


print('e = ', e)
print('mu = ',mu)
print('ang = ',ang)

# we check if exist and open all files .txt one by one and take out graphs results
for i_e in range(len(e)):
	for i_mu in range(len(mu)):
		for i_ang in range(len(ang)):
			#use for all file hdf5
			filename_input = 'Dead_zone_wall_{}_N-{}-e-{}-mu-{}-angle-{}.txt'.format(x_position_from_gate,num_of_grains,str(e[i_e]),str(mu[i_mu]),str(ang[i_ang]))
			if (path.exists(filename_input)==True):
				#open Profil_Vx.txt
				f    = open(filename_input,'r')
				data =f.readlines()
				
				#filter data (attention line in txt begin at 1 while data begin at 0)
				t  = np.array(data[9][1:-2].split(', '),dtype=np.float64)
				dt = t[1]-t[0]
				tank_width = float(data[5])

				#we take results at t = 3.0 (line 1119)
				row_ecart  = coef*(40-11)

				# Dead_zone with V_thresold = coef*Vx
				index_dyn_vx  = int(14+row_ecart)
				#index_line_vx = int(16+row_ecart)
				index_L_vx    = int(18+row_ecart)
				index_S_vx	  = int(20+row_ecart)
				index_num_vx  = int(22+row_ecart)
				index_frac_vx = int(24+row_ecart)

				grains_dead_zone_dyn_vx = data[index_dyn_vx][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				grains_dead_zone_dyn_vx = [p.split(', ') for p in grains_dead_zone_dyn_vx] #split string profil_com to list of string
				grains_dead_zone_dyn_vx = [[np.float64(x) for x in y] for y in grains_dead_zone_dyn_vx]# convert each string to float value
				#print(grains_dead_zone_dyn_vx)
				"""
				grains_dead_zone_line_vx = data[index_line_vx][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				grains_dead_zone_line_vx = [p.split(', ') for p in grains_dead_zone_line_vx] #split string profil_com to list of string
				grains_dead_zone_line_vx = [[np.float64(x) for x in y] for y in grains_dead_zone_line_vx]# convert each string to float value
				"""

				L_dead_zone_vx = np.array(data[index_L_vx][1:-2].split(', '),dtype=np.float64)
				S_dead_zone_vx = np.array(data[index_S_vx][1:-2].split(', '),dtype=np.float64)
				ratio_S_vx 	   = np.array(data[index_frac_vx][1:-2].split(', '),dtype=np.float64)
				
				# Dead_zone with V_thresold = coef*Vnorme (Vnorme = sqrt(Vx² + Vy²))
				index_dyn_vn  = int(28+row_ecart)
				#index_line_vn = int(28+row_ecart)
				index_L_vn    = int(32+row_ecart)
				index_S_vn	  = int(34+row_ecart)
				index_num_vx  = int(36+row_ecart)
				index_frac_vn = int(38+row_ecart)

				grains_dead_zone_dyn_vn = data[index_dyn_vn][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				grains_dead_zone_dyn_vn = [p.split(', ') for p in grains_dead_zone_dyn_vn] #split string profil_com to list of string
				grains_dead_zone_dyn_vn = [[np.float64(x) for x in y] for y in grains_dead_zone_dyn_vn]# convert each string to float value
				"""
				grains_dead_zone_line_vn = data[index_line_vn][2:-3].split('], [') #enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_compacity at each 0.5s)
				grains_dead_zone_line_vn = [p.split(', ') for p in grains_dead_zone_line_vn] #split string profil_com to list of string
				grains_dead_zone_line_vn = [[np.float64(x) for x in y] for y in grains_dead_zone_line_vn]# convert each string to float value
				"""
				
				L_dead_zone_vn = np.array(data[index_L_vn][1:-2].split(', '),dtype=np.float64)
				S_dead_zone_vn = np.array(data[index_S_vn][1:-2].split(', '),dtype=np.float64)
				ratio_S_vn 	   = np.array(data[index_frac_vn][1:-2].split(', '),dtype=np.float64)

				"""
				#We grab 3d string list stocks in one line and slipt it to 3d array 
				grains_dead_zone_dyn_all_steps = data[13][3:-4].split(']], [[') 
				grains_dead_zone_dyn_all_steps = [a.split('], [') for a in grains_dead_zone_dyn_all_steps] 
				grains_dead_zone_dyn_all_steps = [[b.split(', ') for b in a] for a in grains_dead_zone_dyn_all_steps] 
				grains_dead_zone_dyn_all_steps = [[[np.float64(x) for x in y] for y in z] for z in grains_dead_zone_dyn_all_steps]

				grains_dead_zone_line_all_steps = data[15][3:-4]
				grains_dead_zone_line_all_steps = grains_dead_zone_line_all_steps.split(']], [[')
				#print(grains_dead_zone_line_all_steps)
				grains_dead_zone_line_all_steps = [a.split('], [') for a in grains_dead_zone_line_all_steps] 
				grains_dead_zone_line_all_steps = [[b.split(', ') for b in a] for a in grains_dead_zone_line_all_steps] 
				grains_dead_zone_line_all_steps = [[[np.float64(x) for x in y] for y in z] for z in grains_dead_zone_line_all_steps]
				"""
				
				# graph for Vx_thresold                      
				x_grains_vx = (np.array(grains_dead_zone_dyn_vx)[:,2] - tank_width*grain_size)/grain_size
				y_grains_vx = (np.array(grains_dead_zone_dyn_vx)[:,3] -y0)/grain_size
				
				x_grains_vn = (np.array(grains_dead_zone_dyn_vn)[:,2] - tank_width*grain_size)/grain_size
				y_grains_vn = (np.array(grains_dead_zone_dyn_vn)[:,3] -y0)/grain_size

				#Dead_zone graph
				#Vx_thresold
				Dead_zone_fig, dead_zone_graph  = plt.subplots(1,2,figsize=(14, 6))
				Dead_zone_fig.suptitle('Dead_zone-N-{}-e-{}-mu-{}-angle-{} at x = {}, t = {}, coef ={}'.format(num_of_grains,str(e[i_e]),str(mu[i_mu]),str(ang[i_ang]),x_position_from_gate,3,coef))
				dead_zone_graph[0].plot(x_grains_vx, y_grains_vx, 'o', color='green')
				dead_zone_graph[0].set_title('With VX_thresold ')
				dead_zone_graph[0].set_xlabel('x/d')
				dead_zone_graph[0].set_ylabel('y/d')
				dead_zone_graph[0].axis('equal')
				#Vn_thresold
				dead_zone_graph[1].plot(x_grains_vn, y_grains_vn, 'o', color='green')
				dead_zone_graph[1].set_title('With Vnorm_thresold ')
				dead_zone_graph[1].set_xlabel('x/d')
				dead_zone_graph[1].set_ylabel('y/d')
				dead_zone_graph[1].axis('equal')

				#fraction_volume in dead_zone
				ratios_S_fig, ratio_S_graph = plt.subplots(1,2,figsize=(14, 6))
				#Vx_thresold
				ratios_S_fig.suptitle('Compacity-dead_zone-N-{}-e-{}-mu-{}-angle-{} at x = {}, t = {}, coef ={}'.format(num_of_grains,str(e[i_e]),str(mu[i_mu]),str(ang[i_ang]),x_position_from_gate,3,coef))
				ratio_S_graph[0].plot(t, ratio_S_vx, '-o', color='green')
				ratio_S_graph[0].set_title('With VX_thresold ')
				ratio_S_graph[0].set_xlabel('time (s)')
				ratio_S_graph[0].set_ylabel('fraction-volume')
				#Vn_thresold
				ratio_S_graph[1].plot(t, ratio_S_vn, '-o', color='green')
				ratio_S_graph[1].set_title('With Vnorm_thresold ')
				ratio_S_graph[1].set_xlabel('time (s)')
				ratio_S_graph[1].set_ylabel('fraction-volume')


				#L_dead_zone_graph
				dead_zone_size_fig, dead_zone_size_graph  = plt.subplots(2,2,figsize=(14, 6))
				dead_zone_size_fig.suptitle('Dead_zone_size-N-{}-e-{}-mu-{}-angle-{} at x = {}, t = {}, coef ={}'.format(num_of_grains,str(e[i_e]),str(mu[i_mu]),str(ang[i_ang]),x_position_from_gate,3,coef))
				dead_zone_size_graph[0,0].plot(t, L_dead_zone_vx, 'o-', color='green')
				dead_zone_size_graph[0,0].set_title('L_dead_zone with VX_thresold')
				dead_zone_size_graph[0,0].set_xlabel('t (s)')
				dead_zone_size_graph[0,0].set_ylabel('L_dead_zone (x/d)')
			                      
				dead_zone_size_graph[0,1].plot(t, L_dead_zone_vn, 'o-', color='green')
				dead_zone_size_graph[0,1].set_title('L_dead_zone With Vnorm_thresold')
				dead_zone_size_graph[0,1].set_xlabel('t (s)')
				dead_zone_size_graph[0,1].set_ylabel('L_dead_zone (x/d)')

				#S_dead_zone_graph
				dead_zone_size_graph[1,0].plot(t, S_dead_zone_vx, 'o-', color='green')
				dead_zone_size_graph[1,0].set_title('S_dead_zone with VX_thresold')
				dead_zone_size_graph[1,0].set_xlabel('t (s)')
				dead_zone_size_graph[1,0].set_ylabel('S_dead_zone ')

				dead_zone_size_graph[1,1].plot(t, S_dead_zone_vn, 'o-', color='green')
				dead_zone_size_graph[1,1].set_title('S_dead_zone With Vnorm_thresold')
				dead_zone_size_graph[1,1].set_xlabel('t (s)')
				dead_zone_size_graph[1,1].set_ylabel('S_dead_zone ')


				plt.show()