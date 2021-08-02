import h5py
import numpy as np
from numpy import linalg
import math
import matplotlib.pyplot as plt

filename = "2d_sphere_flow_wall_400_N-5000-e-0.5-mu-0.5-angle-30.hdf5"

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

def force_wall_lin(t=0.2,dt=1*hstep): #percussion/dt =force
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
            
            normal_force_wall_lin     = contact_wall_sum[11] /(dt*L_contact_wall)

            tangential_force_wall_lin = contact_wall_sum[12] /(dt*L_contact_wall)
        
            norm_force_wall_lin       = np.linalg.norm(contact_wall_sum[11:13])/(dt*L_contact_wall)
 
  
    return normal_force_wall_lin,tangential_force_wall_lin,norm_force_wall_lin


t1   =1
t2   =3
#dt  =1000*hstep
dt   =0.5
t    = np.arange(t1,t2+dt,dt)

normal_force_wall_all_time     =[]
tangential_force_wall_all_time =[]
norm_force_wall_all_time       =[]

for i in range(len(t)):
    normal_force_wall_lin,tangential_force_wall_lin,norm_force_wall_lin = force_wall_lin(t[i],dt)

    normal_force_wall_all_time.append(norm_force_wall_lin)
    tangential_force_wall_all_time.append(tangential_force_wall_lin)
    norm_force_wall_all_time.append(norm_force_wall_lin)
"""
normal_force_wall_all_time = np.abs(np.asarray(normal_force_wall_all_time))
tangential_force_wall_all_time =np.abs(np.asarray(tangential_force_wall_all_time))
norm_force_wall_all_time = np.abs(np.asarray(norm_force_wall_all_time))
"""
print(t)
print("normal_force_wall_all_time = \n",normal_force_wall_all_time)
print("tangential_force_wall_all_time = \n",tangential_force_wall_all_time)
print("norm_force_wall_all_time = \n",norm_force_wall_all_time)


plt.figure(1)
plt.plot(t,normal_force_wall_all_time,'-.',label = 'normal_force_wall_all_time')
plt.legend()
plt.title('Normal_force_wall versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')


plt.figure(2)
plt.plot(t,tangential_force_wall_all_time,'-',color='black',label = 'tangential_force_wall_all_time')
plt.legend()
plt.title('tangential_force_wall versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')

plt.figure(3)
plt.plot(t,norm_force_wall_all_time,'--',color='red',label = 'norm_force_wall_all_time')
plt.legend()
plt.title('Norm_force_wall versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')
plt.show()