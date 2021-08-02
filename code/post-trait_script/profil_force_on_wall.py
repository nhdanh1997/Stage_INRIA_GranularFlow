import h5py
import numpy as np
from numpy import linalg
import math
import matplotlib.pyplot as plt

filename = "2d_sphere_flow-uniform-mu-0.5-angle-0.52.hdf5"

with h5py.File(filename, "r") as f:

    for key in f.keys():
        # Get the HDF5 group
        group = f[key]

    # group.visit(print) # all the objects and data

    # Checkout what keys(dataset) are inside that group.
    for key in group.keys():
        print(key)

    solv = group['solv'].value
    static = group['static'].value
    v = group['velocities'].value
    # b_cond = group['boundary_conditions'].value
    cf = group['cf'].value
    dyn = group['dynamic'].value


    f.close()


hstep=1e-4

grain_size=1e-3
ground_thickness = 10 * grain_size
ground_y_shift = 0.1
y0 = - ground_thickness - ground_y_shift

def force_wall(t=0.2,dt=1*hstep): #percussion/dt =force
    """
    calcul les forces totales at time t des billes sur les obstacles (tank,wall, ground) 
    idée:
    _ detection de contact : les 2 derniers dans cf[] sont égaux => un des deux sont statiques
    _objet statique : objet A (premier objet)
    _ vérifier le sens de vecteur normal => quel obstacle (wall, ground, tank)
    _ somme de la force / y de contact le plus haut (force linéique)

    retrun : force (normal, tangential, norm) of (tank_front,tank_behind,ground)
    """
    contact            =[] #stock data in cf[] which 2 last value equal (contact with static objet)
    contact_wall       =[] 
    normal_force_wall_lin = tangential_force_wall_lin = norm_force_wall_lin = 0

    contact_t=[]
    raw_times_cf=[]

    k=0
    # Contact detection
    for i in range(len(cf)):
        if(cf[i,-1]==cf[i,-2]): #2 derniers value égaux => un des deux est statique
            contact.append(cf[i,:])
            #print(contact[k][8],contact[k][9],contact[k][10]) #(8,9,10= x,y,zs?)
            k+=1
    contact = np.array(contact)
    #print(len.array(contact)
    #print(len(contact))
    
    for i in range(len(contact)):
        raw_times_cf.append(contact[i,0])

    raw_times_cf  =  np.asarray(raw_times_cf)
    #print(raw_times_cf)
    test_indice_cf_t     = np.logical_and(raw_times_cf>t-hstep, raw_times_cf<t+hstep)
    #return matrix that element be true for all t satisfy
    #print(test_indice_cf_t)

    for i in range(len(contact)):
        if(test_indice_cf_t[i]==True):
            contact_t.append(contact[i,:])

    contact_t=np.array(contact_t)
    #print(len(contact_t))


    #Detection contact with wall
    for i in range(len(contact_t)):
        if(contact_t[i][8]==-1 and contact_t[i][2] > grain_size*150): 
        # normal_x =-1 and position out of range of tank => contact with wall
            contact_wall.append(contact_t[i,:])
    contact_wall = np.asarray(contact_wall)

    return contact_wall



t1   =0
t2   =0.2
#dt  =1000*hstep
dt   =0.05
t    = np.arange(t1,t2+dt,dt)

contact_wall = []
contact_wall = force_wall(t2,dt)

h_ratio      = []
normal_force = []

print(contact_wall)

for i in range(len(contact_wall)):
    h_ratio.append((contact_wall[i,3]-y0)/grain_size)
    normal_force.append(contact_wall[i,12]/dt)

plt.figure(1)
plt.plot(normal_force,h_ratio,'o',label = 'normal_force_wall_all_time')
plt.legend()
plt.title('Normal_force_wall versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')
plt.show()
"""
"""