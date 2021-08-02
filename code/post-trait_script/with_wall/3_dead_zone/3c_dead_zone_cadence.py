import h5py
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import os.path
from os import path
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from statistics import mode
"""
This script will do :
- open one by one hdf5 file
- open one by one Profil_vx. txt file, grab the Profil_vx timle averaged
- caculate :
+ dead_zone forme time evolution
+ L_dead_zone and S_dead_zone time evolution

-stock in dead_zone.txt

"""
# Set attibute to calculate results
# To choose quality of grains in simulation
num_of_grains = 30000

# To choose time simulation to generate results
t_init = 1.
t_end  = 3.
dt     = 0.5
t = np.arange(t_init, t_end+dt, dt)

# To choose position x_wall from gate to the calculation
x_position_from_gate = 700

# Defaults attribute of hdf5 simulation
grain_size       = 0.001
ground_thickness = 10 * grain_size
ground_y_shift   = 0.09
y0 = - ground_thickness - ground_y_shift
hstep = 1e-4

# ordre : e-mu-angle
e_init  = 0.0
e_end   = 0.5
d_e     = 0.25
e = np.arange(e_init, e_end+d_e, d_e)

mu_init = 0.3
mu_end  = 0.7
d_mu    = 0.2
mu = np.arange(mu_init, mu_end+d_mu, d_mu)

ang_init = 20
ang_end  = 20
d_ang    = 2
ang = np.arange(ang_init, ang_end+d_ang, d_ang)

n_row = num_of_grains/100
n_col = 100
distribution      = ('uniform', 0.15)
tank_width_number = n_col+2
tank_thickness    = 10*grain_size
tank_width        = tank_width_number*grain_size*(1.0+distribution[1])

# position of x studied and size of box_dx_dy
#x_position = tank_width/grain_size + x_position_from_gate
y1_ratio   = 0  # start at 0
dx = 1

wall_thickness = 10

print('e = ', e)
print('mu = ', mu)
print('ang = ', ang)


def dead_zone(wall_x_position,V_thresold, t, dx):
    """
    This function detect the grains in dead_zone
    t : time
    dx: step for detection dead_zone_line

    ideal:
    if v < V_thresold => grain quasi-statique situate in dead_zone
       v > V_thresold => grain out of range of dead zone
    V_thresold = 0.05*u with u : the depth-averaged velocity (Profil_Vx averaged following h)

    return  grains_dead_zone_dyn : dynamics of grains in dead zone
            grains_dead_zone_line : coordination of grains in dead zone line
            L_dead_zone: length of dead zone
            S_dead_zone: area of dead zone

    """
    dynt          = []  # sotck data of dynamics a time = t
    raw_times_dyn = []  # stock all time steps  t[] in dynt[]

    vt          = []  # sotck data of velocities a time = t
    raw_times_v = []  # stock all time steps  t[] in dynt[]

    grains_dead_zone_dyn    = [] # stocks dynamics grains detected in dead_zone
    rayon_spheres_dead_zone = [] # stocks rayon of grains in 

    ratio_S     = 0 # fraction_volume of dead_zone
    L_dead_zone = 0
    S_dead_zone = 0

    # Build dynt[]
    raw_times_dyn     = dyn[:, 0]
    test_indice_dyn_t = np.logical_and(raw_times_dyn > t-hstep, raw_times_dyn < t+hstep)
    dynt              = dyn[test_indice_dyn_t]

    # Build vt[]
    raw_times_v     = v[:, 0]
    test_indice_v_t = np.logical_and(raw_times_v > t-hstep, raw_times_v < t+hstep)
    vt              = v[test_indice_v_t]

    # Build grains_dead_zone_dyn
    for i in range(len(vt)):
        if((vt[i, 2] <V_thresold) and (dynt[i, 2] >(tank_width*2)) and (dynt[i, 2] < wall_x_position*grain_size)):
            #(i,'====', vt[i,2])
            grains_dead_zone_dyn.append(dynt[i, :])

    grains_dead_zone_dyn = np.array(grains_dead_zone_dyn)
    #print(len(grains_dead_zone_dyn))
 
    """
    """
    # Calculate L_dead_zone
    x_grains_dead_zone = []

    for i in range(len(grains_dead_zone_dyn)):
        x_grains_dead_zone.append(grains_dead_zone_dyn[i][2])
    x_grains_dead_zone = np.array(x_grains_dead_zone)

    if len(x_grains_dead_zone)!=0:
        L_dead_zone = max(x_grains_dead_zone) - min(x_grains_dead_zone)
        L_dead_zone_grain_size_scaled = L_dead_zone/grain_size
    else:
        L_dead_zone_grain_size_scaled = 0
    x_ratio_dead_zone = np.arange(0, L_dead_zone_grain_size_scaled+1, 1)
    #print(x_ratio_dead_zone)

    # Caculate S_dead_zone by building grains_dead_zone_line
    # Build grains_dead_zone_line
    grains_dead_zone_line = []  # stocks all grains at tallest position in each pas dx

    for j in range(len(x_ratio_dead_zone)-dx):
        x1 = (wall_x_position -x_ratio_dead_zone[j]-wall_thickness*0.5)*grain_size
        x2 = (wall_x_position -x_ratio_dead_zone[j+dx]-wall_thickness*0.5)*grain_size

        # grain_in_dx : we find all grain in each vertical bar in dead_zone
        # to detecte the tallest grains => dead_zone line
        grains_in_dx = []  # stocks grains in each pas dx to find the grain at tallest position
        """
        k = 0
        for i in range(len(grains_dead_zone_dyn)):
            if(grains_dead_zone_dyn[i, 2] < x1 and grains_dead_zone_dyn[i, 2] > x2):
                grains_in_dx.append(grains_dead_zone_dyn[i, :])
                k += 1
                #print(j, "====", k)
        grains_in_dx = np.array(grains_in_dx)
        #print(len(grains_in_dx))
        """
        x_dead_zone_dyn = grains_dead_zone_dyn[:, 2]
        good_x_in_dx    = np.logical_and(x_dead_zone_dyn < x1, x_dead_zone_dyn > x2)
        grains_in_dx    = grains_dead_zone_dyn[good_x_in_dx]

        if(len(grains_in_dx) != 0):
            y_grains_in_dx   = grains_in_dx[:, 3]
            y_tallest_grain  = max(y_grains_in_dx)
            grains_dead_zone_line.append([(x1+x2)*0.5, y_tallest_grain])
    #grains_dead_zone_line = np.array(grains_dead_zone_line)
    #print(grains_dead_zone_line)

    # Caculate S_dead_zone
    for i in range(len(grains_dead_zone_line)-1):
        h = grains_dead_zone_line[i+1][0] - grains_dead_zone_line[i][0]
        # print(h)
        S_trapeze = (grains_dead_zone_line[i+1][1]+grains_dead_zone_line[i][1])*h
        # print(S_parallepipede)
        S_dead_zone += S_trapeze
    S_dead_zone = 0.5*S_dead_zone
    #print('S_dead_zone =',S_dead_zone )

    #calculate fraction volume in dead_zone
    iden = grains_dead_zone_dyn[:, 1] # identity of grains in dead_zone

    for i in range(len(iden)):
        for j in range(len(id_spheres)): 
            if(iden[i] == id_spheres[j]):
                rayon_spheres_dead_zone.append(float(rayon_spheres[j][0]))
    
    # Sum of the squares of radius
    for i in range(len(rayon_spheres_dead_zone)):
        ratio_S += rayon_spheres_dead_zone[i]*rayon_spheres_dead_zone[i] 

    
    ratio_S = mt.pi*ratio_S / S_dead_zone

    return grains_dead_zone_dyn, grains_dead_zone_line, L_dead_zone_grain_size_scaled, S_dead_zone, ratio_S


# we check if exist and open all files hdf5 one by one,
# run the post-trait computation and write to a output file .txt for each simulation
for i_e in range(len(e)):
    for i_mu in range(len(mu)):
        for i_ang in range(len(ang)):
            # use for all file hdf5
            filename_input = '2d_sphere_flow_wall_{}_N-{}-e-{}-mu-{}-angle-{}.hdf5'.format(str(int(x_position_from_gate)), str(num_of_grains), str(e[i_e]), str(mu[i_mu]), str(ang[i_ang]))
            if (path.exists(filename_input) == True):
                with h5py.File(filename_input, "r") as f:
                    ls     = list(f.keys())
                    #f_out.write('list of group : \n',ls,'\n')
                    group  = np.array(f.get('data'))
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
                # read file Profil_Vxy    
                Profil_Vx_input = "Profil_Vxy_"+str(int(x_position_from_gate))+"d"+filename_input[23:-4]+"txt"
                if (path.exists(Profil_Vx_input) == True):
                    f_vx = open(Profil_Vx_input, 'r')
                    data_vx = f_vx.readlines()

                    # enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_vx at each 0.5s)
                    Vx_profil_100steps_averaged = data_vx[24][2:-3].split('], [')# split string profil_vx to list of string
                    Vx_profil_100steps_averaged = [p.split(', ') for p in Vx_profil_100steps_averaged]# convert each string to float value
                    Vx_profil_100steps_averaged = [[np.float64(x) for x in y] for y in Vx_profil_100steps_averaged]

                    # enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_vx at each 0.5s)
                    Vx_profil_instantanenous = data_vx[36][2:-3].split('], [')# split string profil_vx to list of string
                    Vx_profil_instantanenous = [p.split(', ') for p in Vx_profil_instantanenous]# convert each string to float value
                    Vx_profil_instantanenous = [[np.float64(x) for x in y] for y in Vx_profil_instantanenous]

                    # enlever les [] aux 2 bouts et split to 1d list (each string = 1 profil_vx at each 0.5s)
                    Vy_profil_instantanenous = data_vx[53][2:-2].split('], [')# split string profil_vx to list of string
                    Vy_profil_instantanenous = [p.split(', ') for p in Vy_profil_instantanenous]# convert each string to float value
                    Vy_profil_instantanenous = [[np.float64(x) for x in y] for y in Vy_profil_instantanenous]

                    x_position_from_gate = float(data_vx[7]) 
                    x_position = tank_width/grain_size + x_position_from_gate
                    f_vx.close()

                    # creat Dead_zone.txt and write out data
                    dead_zone_ouput = "Dead_zone" + filename_input[14:-4]+"txt"
                    f_out = open(dead_zone_ouput, 'w')
                    f_out.write(filename_input)

                    # Write out defaults attributes of simulation
                    f_out.write('\n\nnum_grains\n')
                    f_out.write(str(num_of_grains))
                    f_out.write("\ntank_width\n")
                    f_out.write(str(tank_width/grain_size))
                    f_out.write('\nx_wall_position_from_gate\n')
                    f_out.write(str(x_position_from_gate))
                    f_out.write("\ntime simulation\n")
                    f_out.write(str(list(t)))

                    #grains_dead_zone_dyn_all_steps  = []
                    #grains_dead_zone_line_all_steps = []
                    
                    c_init = 0
                    c_end  = 0.1
                    d_c    = 0.025
                    coef_thresold = np.arange(c_init,c_end + d_c,d_c)
                    
                    for coef in coef_thresold:                       
                        f_out.write('\n=========== coef_thresold = {} ==========\n'.format(coef))
                        dead_zone_vx = []
                        dead_zone_vn = []
                        num_grains_dead_zone_vx = []
                        num_grains_dead_zone_vn = []

                        L_dead_zone_all_steps_vx = []
                        S_dead_zone_all_steps_vx = []

                        L_dead_zone_all_steps_vn = []
                        S_dead_zone_all_steps_vn = []

                        Vx_thresold_all_steps    = []
                        Vn_thresold_all_steps    = []
                        
                        ratio_S_all_steps_vx     = []
                        ratio_S_all_steps_vn     = []
                        for i in range(len(t)):
                            #f_out.write('\n\n====== t ={} ======='.format(t[i]))
                            Vx_t     = np.array(Vx_profil_instantanenous[i])
                            Vx_t     = Vx_t[np.nonzero(Vx_t)]
                            Vx_t     = np.mean(Vx_t)
                            
                            Vy_t     = np.array(Vy_profil_instantanenous[i])
                            Vy_t     = Vy_t[np.nonzero(Vy_t)]
                            Vy_t     = np.mean(Vy_t)
                        
                            V_norm_t = mt.sqrt(Vx_t*Vx_t+Vy_t*Vy_t)
                            #print(Vx_t,V_norm_t)
                            
                            #detected dead_zone by Vx
                            Vx_thresold = coef*Vx_t
                            grains_dead_zone_dyn_vx, grains_dead_zone_line_vx, L_dead_zone_vx, S_dead_zone_vx, ratio_S_vx = dead_zone(x_position,Vx_thresold, t[i], dx)
                            
                            dead_zone_vx.append(grains_dead_zone_dyn_vx)
                            num_grains_dead_zone_vx.append(len(grains_dead_zone_dyn_vx))
                            Vx_thresold_all_steps.append(Vx_thresold)
                            L_dead_zone_all_steps_vx.append(L_dead_zone_vx)
                            S_dead_zone_all_steps_vx.append(S_dead_zone_vx)
                            ratio_S_all_steps_vx.append(ratio_S_vx)

                            # Detected dead_zone by V_norme
                            Vn_thresold = coef*V_norm_t
                            grains_dead_zone_dyn_vn, grains_dead_zone_line_vn, L_dead_zone_vn, S_dead_zone_vn, ratio_S_vn = dead_zone(x_position,Vn_thresold, t[i], dx)

                            dead_zone_vn.append(grains_dead_zone_dyn_vn)
                            num_grains_dead_zone_vn.append(len(grains_dead_zone_dyn_vn))
                            Vn_thresold_all_steps.append(Vn_thresold)
                            L_dead_zone_all_steps_vn.append(L_dead_zone_vn)
                            S_dead_zone_all_steps_vn.append(S_dead_zone_vn)
                            ratio_S_all_steps_vn.append(ratio_S_vn)

                        f_out.write('\n=== Dead_zone with VX_thresold at 3s ===')
                        f_out.write('\n grains_dead_zone_dyn_vx\n')
                        p = [[float(x) for x in y] for y in dead_zone_vx[-1]]
                        f_out.write(str(p))
                        #f_out.write('\ngrains_dead_zone_line_vx\n')
                        #s = [[float(x) for x in y] for y in grains_dead_zone_line_vx]
                        #f_out.write(str(grains_dead_zone_line_vx))

                        f_out.write('\n Vx_thresold\n')
                        f_out.write(str(Vx_thresold_all_steps))
                        f_out.write('\n L_dead_zone_all_steps_Vx\n')
                        f_out.write(str(L_dead_zone_all_steps_vx))
                        f_out.write('\n S_dead_zone_all_steps_Vx\n')
                        f_out.write(str(S_dead_zone_all_steps_vx))   
                        f_out.write('\n num_grains_dead_zone_vx\n')
                        f_out.write(str(num_grains_dead_zone_vx))
                        f_out.write('\n ratio_S_all_steps_vx\n')
                        f_out.write(str(ratio_S_all_steps_vx))

                        f_out.write('\n\n=== Dead_zone with V_NORME_thresold at 3s ===')
                        f_out.write('\n grains_dead_zone_dyn_vn\n')
                        p1 = [[float(x) for x in y] for y in dead_zone_vn[-1]]
                        f_out.write(str(p1))
                        #f_out.write('\ngrains_dead_zone_line \n')
                        #s1 = [[float(x) for x in y] for y in grains_dead_zone_line_vn]
                        #f_out.write(str(grains_dead_zone_line_vn))
                        
                        f_out.write('\n V_thresold\n')
                        f_out.write(str(Vn_thresold_all_steps))
                        f_out.write('\n L_dead_zone_all_steps_Vn\n')
                        f_out.write(str(L_dead_zone_all_steps_vn))
                        f_out.write('\n S_dead_zone_all_steps_Vn\n')
                        f_out.write(str(S_dead_zone_all_steps_vn))
                        f_out.write('\n num_grains_dead_zone_vn\n')
                        f_out.write(str(num_grains_dead_zone_vn))
                        f_out.write('\n ratio_S_all_steps_vn\n')
                        f_out.write(str(ratio_S_all_steps_vn))

                else:
                    print('File ',Profil_Vx_input,' does not exists')    
            else:
                print('File ',filename_input,' does not exists') 