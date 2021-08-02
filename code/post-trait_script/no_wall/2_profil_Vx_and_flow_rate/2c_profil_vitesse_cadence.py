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
-open hdf5 to take value
-open file flow_height.txt to take value of h_time-averaged and h_instantaneous
-calculate Profil_Vx 100steps time averaged and Profil_Vx instantaneous
-write ou value to Profil_vx.txt
"""
#Set attibute to calculate results

#To choose quality of grains in simulation
num_of_grains = 30000

#To choose position from gate to the calculation 
x_position_from_gate = 300

#To choose time simulation to generate results
t_init  = 1.
t_end   = 3.
dt      = 0.5
t       = np.arange(t_init,t_end+dt,dt)



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

#Defaults attribute of hdf5 simulation 
grain_size       = 0.001
ground_thickness = 10 * grain_size
ground_y_shift   = 0.09
y0               = - ground_thickness - ground_y_shift
hstep = 1e-4

#tank_size identify with 2d_sphere_recirculation.py
n_row = num_of_grains/100
n_col = 100
distribution      = ('uniform', 0.15) 
tank_width_number = n_col+2
tank_width        = tank_width_number*grain_size*(1.0+distribution[1])

#position of x studied and size of box_dx_dy
dx             = 5 #change dx,dy > change size of box computation
dy             = 1 	


print('e = ', e)
print('mu = ',mu)
print('ang = ',ang)
def vitesse_int_moy_dx_dy(x1_ratio=50, dx=5, y1_ratio=0, dy=1, t=0.25):
    """
    #calcul velocity averaged of all grains in the box dx_dy at time t	
    
    return : Velocity averaged following X,Y and Momentum Z
    """

    x1 = x1_ratio*grain_size
    x2 = (x1_ratio + dx)*grain_size

    y1 = y0 + y1_ratio*grain_size
    y2 = y0 + (y1_ratio + dy)*grain_size

    dynt = [] # stock data of dynamics at time t
    vt   = []   # stock data of velocities at time t
    iden = [] # stock iden of billes satisfy the requirement(in the box dx_dy)
    sumVx = sumVy = sumMz = 0

    #dynt[]
    raw_times_dyn     = dyn[:,0]
    test_indice_dyn_t = np.logical_and(raw_times_dyn>t-hstep, raw_times_dyn<t+hstep)
    dynt              = dyn[test_indice_dyn_t]

    #vt[]
    raw_times_v     = v[:,0]
    test_indice_v_t = np.logical_and(raw_times_v>t-hstep, raw_times_v<t+hstep)
    vt              = v[test_indice_v_t]   

    #iden
    iden_raw = dynt[:,1]
    x_dynt   = dynt[:,2]
    y_dynt   = dynt[:,3]

    bool_x_dynt_filtered = (x_dynt>x1)&(x_dynt<x2)
    bool_y_dynt_filtered = (y_dynt>y1)&(y_dynt<y2)

    iden = iden_raw[bool_x_dynt_filtered&bool_y_dynt_filtered]

    if(len(iden) == 0):
        moyenne_Vx = 0
        moyenne_Vy = 0
        moyenne_Mz = 0
    else:
        for i in range(len(iden)):
            # take the grains in vt[] with iden similar to iden[] and calculate the average
            for j in range(len(vt)):
                if(vt[j][1] == iden[i]):
                    sumVx += vt[j][2]
                    sumVy += vt[j][3]
                    sumMz += vt[j][7]
        moyenne_Vx = sumVx/len(iden)
        moyenne_Vy = sumVy/len(iden)
        moyenne_Mz = sumMz/len(iden)

    return moyenne_Vx, moyenne_Vy, moyenne_Mz

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
                
                #open the results flow_height.txt to take value of h_instantaneous and h_timeaveraged
                flow_height_input = "flow_height_"+str(int(x_position_from_gate))+"d"+filename_input[14:-4]+"txt"
                if (path.exists(flow_height_input)==True):    
                    f_h_input         = open(flow_height_input,'r')
                    data_flow_height  = f_h_input.readlines()
                    
                    #flow_height time-averaged and instantaneous calculated by flow_height.py
                    h_time_averaged  = np.array(data_flow_height[39][1:-2].split(','),dtype=np.float64)
                    h_instantaneous  = np.array(data_flow_height[34][1:-2].split(','),dtype=np.float64)
                    f_h_input.close()

                    x_position_from_gate = float(data_flow_height[13])
                    x_position           = tank_width/grain_size + x_position_from_gate 

                    #we calculate all profil_Vx to value h_max 
                    h1_time_averaged_ratio = np.arange(0,max(h_time_averaged.astype(int))+1,1) #make sure these h_ratio are int values 
                    h1_ratio   = np.arange(0,max(h_instantaneous.astype(int))+1,1)
                    
                    filename_ouput = "Profil_Vxy_"+str(int(x_position_from_gate))+"d"+filename_input[14:-4]+"txt"
                    f_out = open(filename_ouput,'w')
                    f_out.write(filename_input)

                    #Write out defaults attributes of simulation
                    f_out.write('\n\nnum_grains\n')
                    f_out.write(str(num_of_grains))
                    f_out.write("\ntank_width\n")
                    f_out.write(str(tank_width/grain_size))
                    f_out.write('\nx_position_from_gate\n')
                    f_out.write(str(x_position_from_gate))
                    f_out.write("\ntime simulation\n")
                    f_out.write(str(list(t)))
                    
                    f_out.write("\n\nh_time_averaged\n")
                    f_out.write(str(list(h_time_averaged)))
                    f_out.write("\nh_instantaneous\n")
                    f_out.write(str(list(h_instantaneous)))                    
                                     

                    f_out.write("\n\n h1_time_averaged_ratio\n")
                    f_out.write(str(list(h1_time_averaged_ratio)))
                    f_out.write("\nh1_ratio\n")
                    f_out.write(str(list(h1_ratio))) 
                    f_out.write("\n\n=====  Profil_vitesse_Vx time-averaged ========\n")

                    Vx_profil_100steps_averaged     =[] #data of all profil of velocity 100steps-averaged
                    Vx_profil_100steps_averaged_std =[] #data of standard deviation
                    Vx_profil_instantanenous        =[]#data of all profil of velocity instantaneous
                    
                    Vy_profil_100steps_averaged     =[] #data of all profil of velocity 100steps-averaged
                    Vy_profil_100steps_averaged_std =[] #data of standard deviation
                    Vy_profil_instantanenous        =[]#data of all profil of velocity instantaneous

                    #caculate Vitesse_int_moy in dx_dy at t[i] and after that, vitesse_moyenne following all t[]
                    for j in range(len(t)):
                        Vx_profil_100steps     = []                 
                        Vx_profil_100steps_std = []
                        
                        Vy_profil_100steps     = []                 
                        Vy_profil_100steps_std = []

                        t_averaged = np.arange(t[j]-100*hstep,t[j]+hstep,hstep) # calculate average in 10 next steps of simulation
                        for y1_ratio in h1_time_averaged_ratio:
                            Vx_100steps =[] # at each y1_ration stocks 100 values of Vx
                            Vy_100steps =[] # at each y1_ration stocks 100 values of Vy
                            for i in range(len(t_averaged)):
                                mVx, mVy, mMz = vitesse_int_moy_dx_dy(x_position, dx, y1_ratio, dy, t_averaged[i])
                                Vx_100steps.append(mVx)
                                Vy_100steps.append(mVy)
        
                            Vx_profil_100steps.append(np.mean(Vx_100steps))
                            Vx_profil_100steps_std.append(np.std(Vx_100steps))
                            
                            Vy_profil_100steps.append(np.mean(Vy_100steps)) 
                            Vy_profil_100steps_std.append(np.std(Vy_100steps))
                        
                        Vx_profil_100steps_averaged.append(Vx_profil_100steps)
                        Vx_profil_100steps_averaged_std.append(Vx_profil_100steps_std)
                        
                        Vy_profil_100steps_averaged.append(Vy_profil_100steps)
                        Vy_profil_100steps_averaged_std.append(Vy_profil_100steps_std)

                    f_out.write("\nVx_profil_100steps_averaged\n") 
                    f_out.write(str(Vx_profil_100steps_averaged)) # matrix 2d (n_row = len(t),n_col = len(h1_ratio)) 
                    f_out.write("\n\nVx_profil_100steps_averaged_std\n") #ecart-type of compacity averaged
                    f_out.write(str(Vx_profil_100steps_averaged_std)) 
                    Vx_profil_100steps_averaged=np.float64(Vx_profil_100steps_averaged)
                

                    #flow_rate_time_averaged
                    flow_rate_time_averaged = []
                    for i in range(len(Vx_profil_100steps_averaged)):
                        flow_rate = np.mean(Vx_profil_100steps_averaged[i])*h_time_averaged[i]*grain_size
                        flow_rate_time_averaged.append(flow_rate)

                    f_out.write("\n\nflow_rate time averaged\n")
                    f_out.write(str(flow_rate_time_averaged))

                    #calculate integration of Vx/h_time_averaged
                    Vx_integration_h_time_averaged = []
                    for i in range(len(Vx_profil_100steps_averaged)):
                        inter = np.mean(Vx_profil_100steps_averaged[i])/(h_time_averaged[i]*grain_size)
                        Vx_integration_h_time_averaged.append(inter)
                    f_out.write("\n\nintegration of Vx/h_time_averaged\n")
                    f_out.write(str(Vx_integration_h_time_averaged))
                    
                    """
                    # Find t when flow_rate become constant    
                    counter_flow_rate_time_averaged    = Counter(flow_rate_time_averaged)
                    const_flow_rate_time_averaged      = mode(flow_rate_time_averaged) # value at this flow_rate become constant

                    iden_flow_rate_10steps_const = np.searchsorted(flow_rate_time_averaged,const_flow_rate_time_averaged)
                    #f_out.write(iden_flow_rate_10steps_const)
                    print(len(t),iden_flow_rate_10steps_const,len(flow_rate_time_averaged))
                    T_const_flow_rate_time_averaged    = t[iden_flow_rate_10steps_const]

                    f_out.write("\n\nflow_rate_10steps_constant \n")
                    f_out.write(str(const_flow_rate_time_averaged))
                    f_out.write("\nnum of value flow_rate_time_averaged_constant \n")
                    f_out.write(str(counter_flow_rate_time_averaged[const_flow_rate_time_averaged]))
                    f_out.write("\nT_const_time_averaged \n")
                    f_out.write(str(T_const_flow_rate_time_averaged))
                    """
                    
                    f_out.write("\n\n=======  Profil_vitesse_Vx instantaneous =======\n")
                    
                    for i in range(len(t)):
                        #f_out.write("\n value at time t = ", t[i])
                        Vx_profil_t =[]
                        Vy_profil_t =[]
                        for y1_ratio in h1_ratio:
                            moyenne_Vx_t,moyenne_Vy_t,moyenne_Mz_t = vitesse_int_moy_dx_dy(x_position, dx, y1_ratio, dy, t[i])
                            #if moyenne_Vx_t!=0 :
                                #f_out.write("value at t = ",t[i],", h = ",y1_ratio)
                            #f_out.write("----------velocity_average by time at y1_ratio =",y1_ratio, " ------------------------")
                            #f_out.write("Vx_instant_average = ", moyenne_Vx_t)
                            Vx_profil_t.append(moyenne_Vx_t)
                            Vy_profil_t.append(moyenne_Vy_t)
                        #f_out.write(Vx_profil_t)
                        Vx_profil_instantanenous.append(Vx_profil_t)    
                        Vy_profil_instantanenous.append(Vy_profil_t)    

                    f_out.write(str(Vx_profil_instantanenous)) # matrix 2d (n_row = len(t),n_col = len(h1_ratio) 
                    Vx_profil_instantanenous=np.float64(Vx_profil_instantanenous)
                    
                    #caclutate time evolution of flow_rate each t
                    flow_rate_instantaneous = []
                    for i in range(len(Vx_profil_instantanenous)):
                        flow_rate = np.mean(Vx_profil_instantanenous[i])*h_instantaneous[i]*grain_size
                        flow_rate_instantaneous.append(flow_rate)
                    f_out.write("\n\nflow_rate_instantaneous\n")
                    f_out.write(str(flow_rate_instantaneous))

                    #calculate integration of Vx/h_instantaeous
                    Vx_integration_h_instantaneous = []
                    for i in range(len(Vx_profil_instantanenous)):
                        inter = np.mean(Vx_profil_instantanenous[i])/(h_instantaneous[i]*grain_size)
                        Vx_integration_h_instantaneous.append(inter)
                    f_out.write("\n\nintegration of Vx/h_instantaeous\n")
                    f_out.write(str(Vx_integration_h_instantaneous))

                    f_out.write("\n\n=====  Profil_vitesse_Vy time-averaged ========\n")

                    f_out.write("\nVy_profil_100steps_averaged\n") 
                    f_out.write(str(Vy_profil_100steps_averaged)) # matrix 2d (n_row = len(t),n_col = len(h1_ratio)) 
                    f_out.write("\n\nVy_profil_100steps_averaged_std\n") #ecart-type of compacity averaged
                    f_out.write(str(Vy_profil_100steps_averaged_std)) 

                    f_out.write("\n\n=======  Profil_vitesse_Vy instantaneous =======\n")
                    f_out.write(str(Vy_profil_instantanenous)) # matrix 2d (n_row = len(t),n_col = len(h1_ratio) 

                    """  
                    # Find t when flow_rate become constant  
                    counter_flow_rate    = Counter(flow_rate_instantaneous)
                    const_flow_rate      = mode(flow_rate_instantaneous) # value at this flow_rate become constant
                    iden_flow_rate_const = np.searchsorted(flow_rate_instantaneous,const_flow_rate)
                    #f_out.write(iden_flow_rate_const)
                    T_const_flow_rate    = t[iden_flow_rate_const]

                    f_out.write("\n\nflow_rate_instantaneous_constant \n")
                    f_out.write(str(flow_rate_instantaneous))
                    f_out.write("\nnum of value flow_rate_instantaneous_constant \n")
                    f_out.write(str(counter_flow_rate[const_flow_rate]))
                    f_out.write("\nT_const_instantaneous \n")
                    f_out.write(str(T_const_flow_rate))
                    """  
                    f_out.close()
                else:
                    print('File ',flow_height_input,' does not exists') 

            else:
                print('File ',filename_input,' does not exists')                    
