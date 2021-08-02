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
-Open one by one 2d_flow_wall.fhd5 (simulation with wall)
-Detecte contact with wall and calculate time evolution of:
+ normal_force_on_wall
+ tangential_force_on_wall
+ norme_force_on_wall
-Stock results in force_on_wall.txt 

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

# ordre : e-mu-angle
e_init  = 0.0
e_end   = 0.5
d_e     = 0.25
e = np.arange(e_init, e_end+d_e, d_e)

mu_init = 0.3
mu_end  = 0.7
d_mu    = 0.2
mu = np.arange(mu_init, mu_end+d_mu, d_mu)

ang_init = 12
ang_end  = 40
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

hstep      = 1e-4
grain_size = 1e-3
dt_force   = hstep*1

def force_wall_lin(t=0.2,d_t=1*hstep): #percussion/dt =force
    """
    calcul les forces totales at time t des billes sur les obstacles (tank,wall, ground) 
    idée:
    _ detection de contact : les 2 derniers dans cf[] sont égaux => un des deux sont statiques
    _ objet statique : objet A (premier objet)
    _ vérifier le sens de vecteur normal => quel obstacle (wall, ground, tank)
    _ somme de la force / y de contact le plus haut (force linéique)

    retrun : force (normal, tangential, norm) of (tank_front,tank_behind,ground)
    """
    contact               = [] #stock data in cf[] which 2 last value equal (contact with static objet)
    contact_wall          = [] 
    normal_force_wall_lin = tangential_force_wall_lin = norm_force_wall_lin = 0

    contact_t    = []
    raw_times_cf = []

    k=0
    # Contact detection
    bool_contact = (cf[:,-1]==cf[:,-2])
    contact = cf[bool_contact]

    #contact detection at t
    raw_times_contact = contact[:,0]
    test_indice_cf_t  = np.logical_and(raw_times_contact>t-hstep, raw_times_contact<t+hstep)
    contact_t         = contact[test_indice_cf_t]

    #Detection contact with wall
    normal_force     = contact_t[:,8]
    position_contact = contact_t[:,2]

    bool_normal_force_on_wall               = (normal_force == -1) #force turns to the left in simulation => (force on wall and force on front_tank)
    bool_position_contact_out_of_tank_width = (position_contact>tank_width*1.1) # these force but out of range of th tank => force on wall
    good_force_on_wall = bool_normal_force_on_wall & bool_position_contact_out_of_tank_width
    contact_wall       = contact_t[good_force_on_wall]

    if (len(contact_wall) == 0):   
        normal_force_wall     = 0
        tangential_force_wall = 0
        norm_force_wall       = 0
    else:
        #chercher les longueurs des contacts pour calculer la force linéique
    
        y_wall =[]

        for i in range (len(contact_wall)): 
            y_wall.append(contact_wall[i][3])
        y_wall         = np.array(y_wall)
        y_wall_max     = max(y_wall)
        y_wall_min     = min(y_wall)
        L_contact_wall = y_wall_max  - y_wall_min
        

        #L_contact minimum = grain_size
        if (L_contact_wall < grain_size):
            L_contact_wall = grain_size

        else:
            #calculer les forces totals sur les CAL    
            contact_wall_sum      = sum(contact_wall[0:len(contact_wall)])     


            #calculer les forces linéiques 
            #ATTENTION: these force is represented in global reference (X,Y,Z) of the system 
            #but not in local reference (n,p,v) at contact point !!!
            #not work for case that Wall inclined => use local reference (cf[20],cf[21]) 
            
            normal_force_wall_lin     = contact_wall_sum[11] /(d_t*L_contact_wall)

            tangential_force_wall_lin = contact_wall_sum[12] /(d_t*L_contact_wall)
        
            norm_force_wall_lin       = np.linalg.norm(contact_wall_sum[11:13])/(d_t*L_contact_wall)
 
  
    return normal_force_wall_lin,tangential_force_wall_lin,norm_force_wall_lin


# we check if exist and open all files hdf5 one by one,
# run the post-trait computation and write to a output file .txt for each simulation
for i_e in range(len(e)):
    for i_mu in range(len(mu)):
        for i_ang in range(len(ang)):
            # use for all file hdf5
            filename_input = '2d_sphere_flow_wall_{}_N-{}-e-{}-mu-{}-angle-{}.hdf5'.format(str(int(x_position_from_gate)), str(num_of_grains), str(e[i_e]), str(mu[i_mu]), str(ang[i_ang]))
            if (path.exists(filename_input) == True):
                with h5py.File(filename_input, "r") as f:

                    ls    = list(f.keys())
                    #f_out.write('list of group : \n',ls,'\n')
                    group = np.array(f.get('data'))
                    #f_out.write('list of datasets : \n', group,'\n')

                    base_items = list(f.items())
                    #f_out.write('infos of items in the base director  : \n', base_items,'\n')
                    data       = f.get('data')
                    data_items = list(data.items())
                    #f_out.write('infos of data in "data" :\n ',data_items,'\n')

                    dyn = np.array(data.get('dynamic'))
                    v   = np.array(data.get('velocities'))
                    cf  = np.array(data.get('cf'))

                    f.close()

                    # creat Dead_zone.txt and write out data
                    force_wall_filename = "force_on_wall_"+str(x_position_from_gate)+"d" + filename_input[23:-4]+"txt"
                    force_wall_ouput     = open(force_wall_filename, 'w')
                    force_wall_ouput.write(filename_input)
                    
                    # Write out defaults attributes of simulation
                    force_wall_ouput.write('\n\nnum_grains\n')
                    force_wall_ouput.write(str(num_of_grains))
                    force_wall_ouput.write("\ntank_width\n")
                    force_wall_ouput.write(str(tank_width/grain_size))
                    force_wall_ouput.write('\nx_wall_position_from_gate\n')
                    force_wall_ouput.write(str(x_position_from_gate))
                    force_wall_ouput.write("\ntime simulation\n")
                    force_wall_ouput.write(str(list(t)))

                    normal_force_wall_inst     =[]
                    tangential_force_wall_inst =[]
                    norm_force_wall_inst       =[]

                    for i in range(len(t)):
                        normal_force_wall_lin,tangential_force_wall_lin,norm_force_wall_lin = force_wall_lin(t[i],dt_force)

                        normal_force_wall_inst.append(abs(normal_force_wall_lin))
                        tangential_force_wall_inst.append(abs(tangential_force_wall_lin))
                        norm_force_wall_inst.append(abs(norm_force_wall_lin))
                    """
                    normal_force_wall_all_time = np.abs(np.asarray(normal_force_wall_all_time))
                    tangential_force_wall_all_time =np.abs(np.asarray(tangential_force_wall_all_time))
                    norm_force_wall_all_time = np.abs(np.asarray(norm_force_wall_all_time))
                    """
                    force_wall_ouput.write("\nnormal_force_wall_inst\n")
                    force_wall_ouput.write(str(normal_force_wall_inst))
                    force_wall_ouput.write("\ntangential_force_wall_int \n")
                    force_wall_ouput.write(str(tangential_force_wall_inst))
                    force_wall_ouput.write("\nnorm_force_wall_inst \n")
                    force_wall_ouput.write(str(norm_force_wall_inst))

                    force_wall_ouput.write("\n\n===== Time_average =====\n")

                    normal_force_wall_time_ave         =[] #normal force_time_averaged 100 steps before
                    normal_force_wall_time_ave_std     =[] #standard deviation
                    
                    tangential_force_wall_time_ave     =[]
                    tangential_force_wall_time_ave_std =[]
                    
                    norm_force_wall_time_ave           =[]
                    norm_force_wall_time_ave_std       =[]

                    for i in range(len(t)):

                        normal_force_100steps         = []
                        normal_force_100steps_std     = []
                        
                        tangentiel_force_100steps     = []
                        tangentiel_force_100steps_std = []
                        
                        norm_force_100steps           = []
                        norm_force_100steps_std       = []

                        t_averaged = np.arange(t[i]-100*hstep,t[i]+hstep,hstep) # calculate average in 10 next steps of simulation
                        
                        for j in range(len(t_averaged)):
                            normal_force_wall_lin_ta,tangential_force_wall_lin_ta,norm_force_wall_lin_ta = force_wall_lin(t_averaged[j],dt_force)
                            
                            normal_force_100steps.append(abs(normal_force_wall_lin_ta))
                            tangentiel_force_100steps.append(abs(tangential_force_wall_lin_ta))
                            norm_force_100steps.append(abs(norm_force_wall_lin_ta))
                            
                        normal_force_100steps     = np.array(normal_force_100steps)
                        tangentiel_force_100steps = np.array(tangentiel_force_100steps)
                        norm_force_100steps       = np.array(norm_force_100steps)
                        normal_force_100steps     = normal_force_100steps[np.nonzero(normal_force_100steps)]
                        tangentiel_force_100steps = tangentiel_force_100steps[np.nonzero(tangentiel_force_100steps)]
                        norm_force_100steps       = norm_force_100steps[np.nonzero(norm_force_100steps)]

                        normal_force_wall_time_ave.append(np.mean(normal_force_100steps)) 
                        normal_force_wall_time_ave_std.append(np.std(normal_force_100steps)) 

                        tangential_force_wall_time_ave.append(np.mean(tangentiel_force_100steps))
                        tangential_force_wall_time_ave_std.append(np.std(tangentiel_force_100steps)) 

                        norm_force_wall_time_ave.append(np.mean(norm_force_100steps)) 
                        norm_force_wall_time_ave_std.append(np.std(norm_force_100steps)) 

                    force_wall_ouput.write("\nnormal_force_wall_time_ave\n")
                    force_wall_ouput.write(str(normal_force_wall_time_ave))
                    force_wall_ouput.write("\n normal_force_wall_time_ave_std\n")
                    force_wall_ouput.write(str(normal_force_wall_time_ave_std))

                    force_wall_ouput.write("\ntangential_force_wall_time_ave \n")
                    force_wall_ouput.write(str(tangential_force_wall_time_ave))
                    force_wall_ouput.write("\n tangential_force_wall_time_ave_std \n")
                    force_wall_ouput.write(str(tangential_force_wall_time_ave_std))

                    force_wall_ouput.write("\nnorm_force_wall_time_ave \n")
                    force_wall_ouput.write(str(norm_force_wall_time_ave))
                    force_wall_ouput.write("\n norm_force_wall_time_ave_std \n")
                    force_wall_ouput.write(str(norm_force_wall_time_ave_std))

                    force_wall_ouput.close()