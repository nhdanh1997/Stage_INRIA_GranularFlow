import h5py
import numpy as np
import math
import matplotlib.pyplot as plt
filename = "2d_sphere_flow-uniform-mu-0.5-angle-0.52.hdf5"

with h5py.File(filename, "r") as f:

    ls    = list(f.keys())
    print('list of group : \n',ls,'\n')
    group =  np.array(f.get('data'))
    print('list of datasets : \n', group,'\n')
    
    base_items = list(f.items())
    print('infos of items in the base director  : \n', base_items,'\n')
    data       = f.get('data')
    data_items = list(data.items())
    #print('infos of data in "data" :\n ',data_items,'\n')
    
    dyn = np.array(data.get('dynamic'))
    v   = np.array(data.get('velocities'))
    cf  = np.array(data.get('cf'))


    ref          = data.get('/data/ref')
    ref_items    = np.array(list(ref.items()))
    #print('infos in "data/ref" :',len(ref_items),'\n')
    #print(ref_items[1][0],"====")

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

grain_size = 0.001
ground_thickness = 10 * grain_size
ground_y_shift = 0.1
y0 = - ground_thickness - ground_y_shift

hstep=1e-4


def fraction_volume_dx_dy(x1_ratio=200,dx=5,y1_ratio=0,dy=1,t=0.):
    
    #return ratio_S = (sum of surface of all grains / dx*dy) at time t 
    # fraction_volume
    
    ratio_S = 0.

    x1 = x1_ratio*grain_size
    x2 = (x1_ratio + dx)*grain_size

    y1 = y0 + y1_ratio*grain_size
    y2 = (y1_ratio + dy)*grain_size
    
    dynt = [] #sotck data of dynamics a time = t
    raw_times_dyn =[] # stock all time steps  t[] in dynt[]

    iden = [] #iden of grains satisfy the requirement in box dx_dy
    rayon_spheres_dxdy =[] # stocks rayon of grains satisfy the requirement in box dx_dy

    

    for i in range(len(dyn)):
        raw_times_dyn.append(dyn[i,0])
    raw_times_dyn = np.asarray(raw_times_dyn)

    test_indice_dyn_t = np.logical_and(raw_times_dyn>t-hstep, raw_times_dyn<t+hstep)
    #return matrix that element be true for all t satisfy
    #print(test_indice_dyn_t)

    for i in range(len(dyn)):
        if(test_indice_dyn_t[i]==True):
            dynt.append(dyn[i,:])
            #print(dynt[i,:])
    
    dynt = np.array(dynt)
    #print("==========")

    for i in range(len(dynt)):
        if (dynt[i,2] > x1 and dynt[i,2] < x2 and dynt[i,3] > y1 and dynt[i,3] < y2):
            # iden: identity of the grains between [x1,x2] at t
            iden.append(dynt[i,1])
            
    #assert(len(iden)!=0):"none of grains between [x1,x2] et this time t= ",t)
       
    #print(iden)

    for i in range(len(iden)):
        for j in range(len(id_spheres)): 
            if(iden[i] == id_spheres[j]):
                rayon_spheres_dxdy.append(float(rayon_spheres[j][0]))
    
    # Sum of the squares of radius
    for i in range(len(rayon_spheres_dxdy)):
        ratio_S += rayon_spheres_dxdy[i]*rayon_spheres_dxdy[i] 


    ratio_S = math.pi*ratio_S / (dx*grain_size*dy*grain_size)

    return ratio_S

    
# Calculate height moyenne h 
"""
x1_ratio : position of x1 scale with grain_size
dx = larger of the box

idea:
-calculer (somme des surface /dx*dy) par rectangle dx_dy de bas en haut jusqu’à quand <0.15 => trouver un hauteur ht
-calculer le moyenne des ht dans tous le temps t

return : h : flow height at this box
"""

porosity   = 0.15


x_position       = 200 #(50 from gate) attention: tank width = 150*grain_size
y1_ratio       = 1 # start at 0
dx             = 5
dy             = 1

t              = 0.8

profil_density       =[] 
height               =[] #(len of h is flow_height at this time and postion)
ratio_S              = 1

while(ratio_S >= porosity):
    ratio_S = fraction_volume_dx_dy(x_position,dx,y1_ratio,dy,t)
    profil_density.append(ratio_S)
    print('ratio_s = ',ratio_S)
    
    y1_ratio+=1 #level up box_dx_dy until reach flow_heigt (ratio_s < 0.15)
    height.append(y1_ratio)


plt.figure(1)
plt.plot(profil_density,height,'o',label = ' profil de Vx')
plt.legend()
plt.title('Profil_cpmpacity at x_pos={} band t ='.format(x_position-150,t))
plt.xlabel('Fraction_volume')
plt.ylabel('y/d')
plt.show()

