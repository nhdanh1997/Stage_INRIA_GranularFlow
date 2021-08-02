import h5py
import numpy as np
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
    #b_cond = group['boundary_conditions'].value
    cf = group['cf'].value
    dyn = group['dynamic'].value
 
    f.close()

hstep=1e-4

def vitesse_int_moy_dx(x1_ratio=200, dx=5, t=0.25):
    """
    #calcul instantaneous average velocity of grains between x1 and x2 at time t	
    x1= 200 (150(tank) + 50)
    """
    grain_size = 0.001
    x1 = x1_ratio*grain_size
    x2 = (x1_ratio + dx)*grain_size

    dynt = [] #stocks data of dynamics at time = t
    vt   = [] #stocks data of velocities at time = t
    
    raw_times_dyn =[] # stock all time steps t in dynt[]
    raw_times_v  =[] # stock all time steps t  in vt[]

    iden = [] #iden of grains satisfy the requirement in box dx_dy


    k1 = k2 = k3 = 0
    sumVx = sumVy = sumMz = 0
    print("length of dynamics and velocity")
    print(len(dyn), 'and', len(v))

    print("-----dynt[]------")
    for i in range(len(dyn)):
        raw_times_dyn.append(dyn[i,0])
    raw_times_dyn = np.asarray(raw_times_dyn)

    test_indice_dyn_t  = np.logical_and(raw_times_dyn>t-hstep, raw_times_dyn<t+hstep)
    #return matrix that element be true for all t satisfy

    for i in range(len(dyn)):
        if(test_indice_dyn_t[i]==True):
            dynt.append(dyn[i,:])
    dynt = np.array(dynt)
    print(dynt)


    print("-----vt[]------")
    for i in range(len(v)):
        raw_times_v.append(v[i,0])
    raw_times_v = np.asarray(raw_times_v)

    test_indice_vt = np.logical_and(raw_times_v>t-hstep, raw_times_v<t+hstep)
    #return matrix that element be true for all t satisfy

    for i in range(len(v)):
        if(test_indice_vt[i]==True):
            vt.append(v[i,:])
    vt = np.array(vt)
    print(vt)


    print("-------iden[] of grains at t and between [x1,x2]--------")
    for i in range(len(dynt)):
        if (dynt[i][2] > x1 and dynt[i][2] < x2):
            # iden: identity of the grains between [x1,x2] at t
            iden.append(dynt[i][1])
    #assert (len(iden) != 0), "none of grains between [x1,x2] et this time t"

    print(iden)

    if(len(iden) == 0):
        moyenne_Vx = 0
        moyenne_Vy = 0
        moyenne_Mz = 0
    else:
        for i in range(len(iden)):
            # take the grains in vt[] with iden similar to iden[] and calculate the average
            for j in range(len(vt)):
                if(vt[j][1] == iden[i]):
                    sumVx += vt[j,2]
                    sumVy += vt[j,3]
                    sumMz += vt[j,7]
        moyenne_Vx = sumVx/len(iden)
        moyenne_Vy = sumVy/len(iden)
        moyenne_Mz = sumMz/len(iden)

    return moyenne_Vx, moyenne_Vy, moyenne_Mz

x_pos = 200
dx = 5

t1   =0.1
t2   =1.
dt   =0.1
t    = np.arange(t1,t2+dt,dt)

Vx_average_all_step =[]
Vy_average_all_step =[]
Mz_average_all_step =[]


for i in range(len(t)):
    moyenne_Vx, moyenne_Vy, moyenne_Mz = vitesse_int_moy_dx(x_pos,dx,t[i])
    print("\n---------------t =",t[i] ," ------------------- \n")
    print("Vx_instant_average = ", moyenne_Vx)
    print("Vy_instant_average = ", moyenne_Vy)
    print("Mz_instant_average = ", moyenne_Mz)
    Vx_average_all_step.append(moyenne_Vx)
    Vy_average_all_step.append(moyenne_Vy)
    Mz_average_all_step.append(moyenne_Mz)

plt.figure(1)
plt.bar(t,Vx_average_all_step,align="center",width=0.03,label='Vx_average')
plt.legend()
plt.title('time evolution of Vx_average with dx ={}, dt = {}'.format(dx,dt))
plt.xlabel('time t(s)')
plt.ylabel('Vx_average')
plt.show();

plt.figure(2)
plt.bar(t,Vy_average_all_step,align="center",width=0.03,label='Vy_average')
plt.legend()
plt.title('time evolution of Vy_average with dx ={}, dt = {}'.format(dx,dt))
plt.xlabel('time t(s)')
plt.ylabel('Vy_average')
plt.show();

plt.figure(3)
plt.bar(t,Mz_average_all_step,align="center",width=0.03,label='Mz_average')
plt.legend()
plt.title('time evolution of Mz_average with dx ={}, dt = {}'.format(dx,dt))
plt.xlabel('time t(s)')
plt.ylabel('Mz_average')
plt.show();