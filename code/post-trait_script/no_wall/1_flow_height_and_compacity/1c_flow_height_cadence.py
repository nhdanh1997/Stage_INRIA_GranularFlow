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
x1_ratio : position of x1 scale with grain_size (x1 = grain_size*x1_ratio)
(x1_ratio,y1_ratio) : point le plus bas à gauche of box test (here we use box with area=5x1)
dx = larger of the box
dy = hauteur of box

idea algorithm:
-calculer (somme des surface /dx*dy) par rectangle dx_dy de bas en haut jusqu’à quand <0.15 => trouver un hauteur h_t for a moment t

This srcipt calculate:
(1) profil of 100step time-averaged compacity for each 0.5s (we calculate mean of 100 value ratio_S(compacity) => if ratio_S<0.15 so lift up y1_ratio)
(2) h_compacity_averaged ( calculate by (1))
(3) profil_compacity (density) instantaneous for each  0.5s 
(4) h_instantaneous at each 0.5s
(5) h_time_averaged (we calculate 100 profil_density => 100 values of h => 1 value averaged with ecarte type for each 0.5s)

"""
# Set attibute to calculate results
#To choose quality of grains in simulation
num_of_grains = 30000

#To choose position from gate to the calculation 
x_position_from_gate = 300

#To choose time simulation to generate results
t_init  = 1.
t_end   = 3.
dt      = 0.5
t       = np.arange(t_init,t_end+dt,dt)

#Defaults attribute of hdf5 simulation 
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

n_row = num_of_grains/100
n_col = 100
distribution      = ('uniform', 0.15) 
tank_width_number = n_col+2
tank_width        = tank_width_number*grain_size*(1.0+distribution[1])

#position of x studied and size of box_dx_dy
x_position     = tank_width/grain_size + x_position_from_gate 
y1_ratio       = 0 # start at 0
dx             = 5
dy             = 1 	

porosity = 0.15

print('e = ', e)
print('mu = ',mu)
print('ang = ',ang)

def fraction_volume_dx_dy(x1_ratio=200,dx=5,y1_ratio=0,dy=1,t=0.):
    
    #return ratio_S = (sum of surface of all grains / dx*dy) at time t 
    #fraction_volume
    
    ratio_S = 0.

    x1 = x1_ratio*grain_size
    x2 = (x1_ratio + dx)*grain_size
    #f_out.write(x1,"====",x2)

    y1 = y0 + y1_ratio*grain_size
    y2 = y0 + (y1_ratio + dy)*grain_size
    #f_out.write(y1,"====",y2)
    
    raw_times_dyn =[] # stock all time steps  t[] in dynt[]
    dynt = [] #sotck data of dynamics a time = t

    iden = [] #iden of grains satisfy the requirement in box dx_dy
    rayon_spheres_dxdy =[] # stocks rayon of grains in box dx_dy

    #dynt[]
    raw_times_dyn = dyn[:,0]
    test_indice_dyn_t = np.logical_and(raw_times_dyn>t-hstep, raw_times_dyn<t+hstep)
    dynt = dyn[test_indice_dyn_t]
    
    #iden[]
    iden_raw = dynt[:,1]
    x_dynt   = dynt[:,2]
    y_dynt   = dynt[:,3]

    bool_x_dynt_filtered = (x_dynt>x1)&(x_dynt<x2)
    bool_y_dynt_filtered = (y_dynt>y1)&(y_dynt<y2)

    iden = iden_raw[bool_x_dynt_filtered&bool_y_dynt_filtered]
    #f_out.write("num of grain in box",len(iden))
    #dynt_good = dynt[bool_x_dynt_filtered&bool_y_dynt_filtered]
    #f_out.write(dynt_good)

    for i in range(len(iden)):
        for j in range(len(id_spheres)): 
            if(iden[i] == id_spheres[j]):
                rayon_spheres_dxdy.append(float(rayon_spheres[j][0]))
    
    # Sum of the squares of radius
    for i in range(len(rayon_spheres_dxdy)):
        ratio_S += rayon_spheres_dxdy[i]*rayon_spheres_dxdy[i] 

    
    ratio_S = math.pi*ratio_S / (dx*grain_size*dy*grain_size)

    return ratio_S


# we check if exist and open all files hdf5 one by one, 
# run the post-trait computation and write to a output file .txt for each simulation 
for i_e in range(len(e)):
	for i_mu in range(len(mu)):
		for i_ang in range(len(ang)):
			#filename_input = 2d_sphere_flow_N-30000-e-0.5-mu-0.5-angle-30.hdf5 #use for 1 file hdf5
			#use for all file hdf5
			filename_input = '2d_sphere_flow_N-{}-e-{}-mu-{}-angle-{}.hdf5'.format(num_of_grains,str(e[i_e]),str(mu[i_mu]),str(ang[i_ang]))
			if (path.exists(filename_input)==True):
				with h5py.File(filename_input, "r") as f:

				    ls    = list(f.keys())
				    #f_out.write('list of group : \n',ls,'\n')
				    group =  np.array(f.get('data'))
				    #f_out.write('list of datasets : \n', group,'\n')
				    
				    base_items = list(f.items())
				    #f_out.write('infos of items in the base director  : \n', base_items,'\n')
				    data       = f.get('data')
				    data_items = list(data.items())
				    #f_out.write('infos of data in "data" :\n ',data_items,'\n')
				    
				    dyn = np.array(data.get('dynamic'))
				    v   = np.array(data.get('velocities'))
				    cf  = np.array(data.get('cf'))


				    ref          = data.get('/data/ref')
				    ref_items    = np.array(list(ref.items()))
				    #f_out.write('infos in "data/ref" :',len(ref_items),'\n')
				    #f_out.write(ref_items[1][0],"====")

				    rayon_spheres=[] # stock radius of all spheres
				    name_spheres =[] # stock nams of all spheres
				    id_spheres   =[] # stock id of all spheres

				    for i in range(len(ref_items)):
				        

				        rayon = np.array(ref.get(ref_items[i][0]))
				        rayon_spheres.append(rayon)
				        
				        name   = ref_items[i][0]
				        name_spheres.append(name)
				        
				        director =  '/data/ref/'+name
				        iden = f[director].attrs['id']
				        id_spheres.append(iden) 

				    f.close()
				
				#create file results .txt
				filename_ouput = "flow_height_"+str(int(x_position_from_gate))+"d"+filename_input[14:-4]+"txt"
				f_out = open(filename_ouput,'w')
				f_out.write(filename_input)

				f_out.write('\n\ngrain_size\n')
				f_out.write(str(grain_size))
				f_out.write('\nhstep\n')
				f_out.write(str(hstep))

				f_out.write('\nporosity\n')
				f_out.write(str(porosity))

				h                    =[] # stock all flow height instantaneous at every t (each 0.5s)

				h_time_averaged      =[] # stock all flow_height_time_average by 100steps at every t (100 profil_compacity => 100 h => h_time_averaged)
				h_time_average_std   =[] # stocks standard deviation (ecart-type) for each h_time_averaged

				h_compacity_averaged = [] # lift up when ratio_S_100steps_averaged < 0.15 (1 profil_compacity by averaging ratio_S_100steps for lifting-up each box dx_dy form bottom)

				profil_density_instantaneous =[] #data of all profil_density instantaneous each 0.5s
				profil_compacity_100steps_averaged =[] #data of all profil of compacity 100step time-averaged for each 0.5s (ratio_S_100steps_averaged)
				profil_compacity_100steps_std_all  =[]# data of standard deviation of 100steps averaged compacity

				#tank_size identify with 2d_sphere_recirculation.py
				f_out.write('\nnum_grains\n')
				f_out.write(str(num_of_grains))

				f_out.write("\ntank_width\n")
				f_out.write(str(tank_width/grain_size))

				f_out.write('\nx_position_from_gate\n')
				f_out.write(str(x_position_from_gate))

				f_out.write("\ntime simulation\n")
				f_out.write(str(t))

				f_out.write("\n\n=============Profil_compacity_100steps_averaged =============\n")
				# (1) and (2) profil of compacity time-averaged
				# for each 0,5s, we calculate 100step time_averaged of porosity each box dx_dy from bottom to top , if ratio_S_100steps_averaged <0.15 => lift up the flow_height

				for i in range(len(t)):
				    ratio_S = 1  # >porosity to loop while  
				    profil_compacity_100steps     =[] # each 0.5 t, stocks profil of 100steps averaged value of compacity form bottom to flow_height
				    profil_compacity_100steps_std =[] # each 0.5 t, stocks profil of standard deviation of 100steps averaged value of compacity form bottom to flow_height
				 
				    while(ratio_S >= porosity):
				        compacity_100steps        =[] # stocks 100 compacitys for 100steps
				        t_averaged = np.arange(t[i]-100*hstep,t[i]+hstep,hstep)
				        for j in range(len(t_averaged)):
				            ratio_S= fraction_volume_dx_dy(x_position,dx,y1_ratio,dy,t_averaged[j])
				            compacity_100steps.append(ratio_S)

				        profil_compacity_100steps.append(np.mean(compacity_100steps)) #from bottom to top, stocks all mean values compacities
				        profil_compacity_100steps_std.append(np.std(compacity_100steps))#stocks standard deviation (ecart-type) of each mean value compacity
				        y1_ratio += 1
				    
				    profil_compacity_100steps_averaged.append(profil_compacity_100steps)
				    profil_compacity_100steps_std_all.append(profil_compacity_100steps_std)
				    h_compacity_averaged.append(y1_ratio)

				    y1_ratio = 1 # reset y1 to next compacity computation

				f_out.write("\nprofil_compacity_100steps_averaged\n")    
				f_out.write(str(profil_compacity_100steps_averaged))
				f_out.write("\n\nprofil_compacity_100steps_std_all\n") #ecart-type of compacity averaged
				f_out.write(str(profil_compacity_100steps_std_all))
				f_out.write("\n\nh_compacity_averaged\n")
				f_out.write(str(h_compacity_averaged))

				f_out.write("\n\n====================== h_instantaneous ======================\n")
				# (3) and (4) we calculate flow_height and profil_compacity instantaneous for each 0.5s
				for i in range(len(t)):# calculate flow_height at every t[]
				    ratio_S = 1 # >porosity to loop while  
				    profil_density =[]
				    #f_out.write("\n ===== at time t =  ", t[i] )
				    while(ratio_S >= porosity):
				        ratio_S = fraction_volume_dx_dy(x_position,dx,y1_ratio,dy,t[i])
				        profil_density.append(ratio_S)
				        #f_out.write('y/d = ',y1_ratio,', ratio_s = ',ratio_S)
				        y1_ratio+=1 #level up box_dx_dy until reach flow_heigt (ratio_s < 0.15)
				    
				    profil_density_instantaneous.append(profil_density)
				    h.append(y1_ratio)

				    y1_ratio = 1  

				h_average = sum(h)/len(h)
				h_max     = max(h)

				f_out.write("\nprofil_density_instantaneous\n")
				f_out.write(str(profil_density_instantaneous))
				f_out.write("\n\nh_instantaneous\n")
				f_out.write(str(h))


				f_out.write("\n\n========================== h_time-averaged ==================\n")
				# (5) we calculate flow_height for each step in 100step before each 0.5s and calculate finally h_100steps_averaged 
				ratio_S  = 1 # >porosity to loop while  
				y1_ratio = 1
				#time moyenne each pas
				for i in range(len(t)):
				    t_averaged = np.arange(t[i]-100*hstep,t[i]+hstep,hstep)
				    h_100steps = []
				    for j in range(len(t_averaged)):# calculate flow_height at every t[]
				        #f_out.write("\n ===== at time t =  ", t_averaged[j] )
				        while(ratio_S >= porosity):
				            ratio_S = fraction_volume_dx_dy(x_position,dx,y1_ratio,dy,t_averaged[j])
				            #f_out.write('y/d = ',y1_ratio,', ratio_s = ',ratio_S)
				            y1_ratio+=1 #level up box_dx_dy until reach flow_heigt (when ratio_s < 0.15)
				        
				        h_100steps.append(y1_ratio)
				        ratio_S  = 1  # reset value to next time t[i] computation
				        y1_ratio = 1
				    h_time_averaged.append(np.mean(h_100steps)) #stock all h_averaged of 100 values h
				    h_time_average_std.append(np.std(h_100steps)) #stock ecart-type of each h_averaged 

				f_out.write("\nh_time_averaged\n")
				f_out.write(str(h_time_averaged))
				f_out.write("\n\nh_time_average_std\n")
				f_out.write(str(h_time_average_std))

				f_out.close()
			else :
				print('File ',filename_input,' does not exists')
