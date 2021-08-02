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

def force(t=0.2,dt=1*hstep): #percussion/dt =force
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
    contact_t          =[] #stock data in contact[] at time t
    contact_ground     =[] #stock data in contact_t[] wich contact with ground 
    contact_wall       =[] 
    contact_tank_front =[]
    contact_tank_behind=[]
    contact_gate       =[]

    contact_t=[]
    raw_times_cf=[]

    normal_force_tank_front_lin     = normal_force_tank_behind_lin     = normal_force_ground_lin    =0    
    tangential_force_tank_front_lin = tangential_force_tank_behind_lin = tangential_force_ground_lin=0
    norm_force_tank_front_lin       = norm_force_tank_behind_lin       = norm_force_ground_lin            =0
    
    
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


    #Detection contact with wich obstacle
    for i in range(len(contact_t)):
        if(contact_t[i][8]==-1 and contact_t[i][2]<200*grain_size): 
        # normal_x =-1 and in range of tank => contact with tank_front
            contact_tank_front.append(contact_t[i,:])
        elif(contact_t[i][8]==1):# normal_x =1 => contact with tank_behind
            contact_tank_behind.append(contact_t[i,:])
        elif(contact_t[i][9]==1):# normal_y 1 => contact with ground
            contact_ground.append(contact_t[i,:])    
    

    if (len(contact_tank_front) == 0):   
        normal_force_tank_front_lin     = 0
        tangential_force_tank_front_lin = 0
        norm_force_tank_front_lin = 0

    elif (len(contact_tank_behind) == 0):
        normal_force_tank_behind_lin = 0
        tangential_force_tank_behind_lin = 0
        norm_force_tank_behind_lin = 0

    elif (len(contact_ground) == 0):
        normal_force_ground_lin = 0
        tangential_force_ground_lin = 0
        norm_force_ground_lin = 0

    else:
        #chercher les longueurs des contacts pour calculer la force linéique
    
        y_front =[]
        y_behind=[]
        x_ground=[]

        for i in range (len(contact_tank_front)): 
            y_front.append(contact_tank_front[i][3])
        y_front         = np.array(y_front)
        y_front_max     = max(y_front)
        y_front_min     = min(y_front)
        L_contact_front = y_front_max  - y_front_min
        
        for i in range (len(contact_tank_behind)): 
            y_behind.append(contact_tank_behind[i][3])
        y_behind         = np.array(y_behind)
        y_behind_max     = max(y_behind)
        y_behind_min     = min(y_behind)
        L_contact_behind = y_behind_max - y_behind_min

        for i in range(len(contact_ground)):
            x_ground.append(contact_ground[i][2])
        x_ground         = np.array(x_ground)
        x_ground_max     = max(x_ground)
        x_ground_min     = min(x_ground)
        L_contact_ground = x_ground_max - x_ground_min

        if (L_contact_front < grain_size):
            L_contact_front = grain_size
        elif (L_contact_behind < grain_size):
            L_contact_behind = grain_size
        elif (L_contact_ground < grain_size):
            L_contact_ground = grain_size
        else:
            #calculer les forces totals sur les CAL
            contact_tank_front_sum  = sum(contact_tank_front[0:len(contact_tank_front)])  
            contact_tank_behind_sum = sum(contact_tank_behind[0:len(contact_tank_behind)])
            contact_ground_sum      = sum(contact_ground[0:len(contact_ground)])     

            #calculer les forces linéiques 
            #ATTENTION: these force is represented in global reference (X,Y,Z) of the system 
            #but not in local reference (n,p,v) at contact point !!!
            #not work for case that Wall inclined => use local reference (cf[20],cf[21]) 
            normal_force_tank_front_lin    = contact_tank_front_sum[11] /(dt*L_contact_front)
            normal_force_tank_behind_lin   = contact_tank_behind_sum[11]/(dt*L_contact_behind)
            normal_force_ground_lin        = contact_ground_sum[12]     /(dt*L_contact_ground)
        
            tangential_force_tank_front_lin    = contact_tank_front_sum[12] /(dt*L_contact_front)
            tangentail_force_tank_behind_lin   = contact_tank_behind_sum[12]/(dt*L_contact_behind)
            tangential_force_ground_lin        = contact_ground_sum[11]     /(dt*L_contact_ground)
            
            norm_force_tank_front_lin    = np.linalg.norm(contact_tank_front_sum[11:13]) /(dt*L_contact_front)
            norm_force_tank_behind_lin   = np.linalg.norm(contact_tank_behind_sum[11:13])/(dt*L_contact_behind)
            norm_force_ground_lin        = np.linalg.norm(contact_ground_sum[11:13])     /(dt*L_contact_ground)
    """

            print("======les longeurs de contact ========")
            print("L_contact_front = ", L_contact_front)
            print("L_contact_behind = ", L_contact_behind)
            print("L_contact_ground = ", L_contact_ground)

        print("========= les forces totals sur les CAL========")
        print("num_contact_tank_front = ",len(contact_tank_front))
        print("force_tank_front_sum selon x,y ",contact_tank_front_sum[11],contact_tank_front_sum[12])
        print("num_contact_tank_behind = ",len(contact_tank_behind))
        print("force_tank_behind_sum selon x,y ",contact_tank_behind_sum[11],contact_tank_behind_sum[12])
        print("num_contact_ground = ",len(contact_ground))
        print("force_ground_sum selon x,y = ",contact_ground_sum[11],contact_ground_sum[12])


        print("========= les forces linéiques sur les CAL========")
        print("force_tank_front_lin = ",normal_force_tank_front_lin)
        print("force_tank_behind_lin = ",normal_force_tank_behind_lin)
        print("force_ground_lin = ",normal_force_ground_lin)
    """
    return normal_force_tank_front_lin, normal_force_tank_behind_lin, normal_force_ground_lin, \
           tangential_force_tank_front_lin, tangential_force_tank_behind_lin, tangential_force_ground_lin, \
           norm_force_tank_front_lin, norm_force_tank_behind_lin, norm_force_ground_lin


t1   =0
t2   =1
#dt  =1000*hstep
dt   =0.05
t    = np.arange(t1,t2+dt,dt)

normal_force_tank_front_all_time  =[]
normal_force_tank_behind_all_time =[] 
normal_force_ground_all_time      =[]

tangential_force_tank_front_all_time =[]
tangential_force_tank_behind_all_time=[] 
tangential_force_ground_all_time     =[]

norm_force_tank_front_all_time =[]
norm_force_tank_behind_all_time=[] 
norm_force_ground_all_time     =[]

normal_force_tank_front_lin     = normal_force_tank_behind_lin     = normal_force_ground_lin    =0    
tangential_force_tank_front_lin = tangential_force_tank_behind_lin = tangential_force_ground_lin=0
norm_force_tank_front_lin       = norm_force_tank_behind_lin       = norm_force_ground_lin      =0

for i in range(len(t)):
    normal_force_tank_front_lin,normal_force_tank_behind_lin,normal_force_ground_lin, \
    tangential_force_tank_front_lin, tangential_force_tank_behind_lin, tangential_force_ground_lin, \
    norm_force_tank_front_lin, norm_force_tank_behind_lin, norm_force_ground_lin  =force(t[i],dt)
    #print("============", t[i], "===========")
    #Normal force
    normal_force_tank_front_all_time.append(normal_force_tank_front_lin)
    normal_force_tank_behind_all_time.append(normal_force_tank_behind_lin)
    normal_force_ground_all_time.append(normal_force_ground_lin)
    #tangential force
    tangential_force_tank_front_all_time.append(tangential_force_tank_front_lin)
    tangential_force_tank_behind_all_time.append(tangential_force_tank_behind_lin)
    tangential_force_ground_all_time.append(tangential_force_ground_lin)
    #norm force 
    norm_force_tank_front_all_time.append(norm_force_tank_front_lin)
    norm_force_tank_behind_all_time.append(norm_force_tank_behind_lin)
    norm_force_ground_all_time.append(norm_force_ground_lin)


normal_force_tank_front_all_time  = np.abs(np.asarray(normal_force_tank_front_all_time))
normal_force_tank_behind_all_time = np.abs(np.asarray(normal_force_tank_behind_all_time))
normal_force_ground_all_time      = np.abs(np.asarray(normal_force_ground_all_time))

tangential_force_tank_front_all_time  = np.abs(np.asarray(tangential_force_tank_front_all_time))
tangential_force_tank_behind_all_time = np.abs(np.asarray(tangential_force_tank_behind_all_time))
tangential_force_ground_all_time      = np.abs(np.asarray(tangential_force_ground_all_time))

norm_force_tank_front_all_time  = np.abs(np.asarray(norm_force_tank_front_all_time))
norm_force_tank_behind_all_time = np.abs(np.asarray(norm_force_tank_behind_all_time))
norm_force_ground_all_time      = np.abs(np.asarray(norm_force_ground_all_time))

print(t)
print("normal_force_tank_front_all_time = \n",normal_force_tank_front_all_time)
print("normal_force_tank_behind_all_time = \n",normal_force_tank_behind_all_time)
print("normal_force_ground_all_time = \n",normal_force_ground_all_time)

#normal_force
plt.figure(1)
plt.plot(t,normal_force_tank_front_all_time,'-.',label = 'normal_force_tank_front_lin')
plt.plot(t,normal_force_tank_behind_all_time,'-',color='black',label = 'normal_force_tank_behind_lin')
plt.plot(t,normal_force_ground_all_time,'--',color='red',label = 'normal_force_ground_lin')
plt.legend()
plt.title('Normal_force versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')

#tangential_force
plt.figure(2)
plt.plot(t,tangential_force_tank_front_all_time,'-.',label = 'tangential_force_tank_front_lin')
plt.plot(t,tangential_force_tank_behind_all_time,'-',color='black',label = 'tangential_force_tank_behind_lin')
plt.plot(t,tangential_force_ground_all_time,'--',color='red',label = 'tangential_force_ground_lin')
plt.legend()
plt.title('Tangential_force versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')

#norm_force
plt.figure(3)
plt.plot(t,norm_force_tank_front_all_time,'-.',label = 'norm_force_tank_front_lin')
plt.plot(t,norm_force_tank_behind_all_time,'-',color='black',label = 'norm_force_tank_behind_lin')
plt.plot(t,norm_force_ground_all_time,'--',color='red',label = 'norm_force_ground_lin')
plt.legend()
plt.title('Norm_force versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')
plt.show();

#normal_force:separate graphs
plt.figure(4)
plt.plot(t,normal_force_tank_front_all_time,'-.',label = 'force_tank_front_lin')
plt.legend()
plt.title('Normal_force_tank_front versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')

plt.figure(5)
plt.plot(t,normal_force_tank_behind_all_time,'-',color='black',label = 'force_tank_behind_lin')
plt.legend()
plt.title('Normal_force_tank_behind versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')

plt.figure(6)
plt.plot(t,normal_force_ground_all_time,'--',color='red',label = 'force_ground_lin')
plt.legend()
plt.title('Normal_force_ground versus time with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('F (N.m-1)')

plt.show();