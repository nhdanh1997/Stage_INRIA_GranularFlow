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

grain_size = 0.001
ground_thickness = 10 * grain_size
ground_y_shift = 0.1
y0 = - ground_thickness - ground_y_shift

x_position=200 #ATTENTION: TANK SIZE =150 => 200 ~ 50 FROM GATE OF TANK)
hstep=1e-4


def vitesse_int_moy_dx_dy(x1_ratio=x_position, dx=5, y1_ratio=0, dy=1, t=0.25):
    """
    #calcul instantaneous average velocity of grains bin dx_dy at time t	
    x,y_ratio : scale for grain_size meanwhile x,y :real value position
    x1_ratio= 200 (150(tank) + 50)
    """

    x1 = x1_ratio*grain_size
    x2 = (x1_ratio + dx)*grain_size

    y1 = y0 + y1_ratio*grain_size
    y2 = (y1_ratio + dy)*grain_size

    dynt = [] # stock data of dynamics at time t
    vt   = []   # stock data of velocities at time t
    iden = [] # stock iden of billes satisfy the requirement(in the box dx_dy)
    k1 = k2 = k3 = 0
    sumVx = sumVy = sumMz = 0
    #print("length of dynamics and velocity")
    #print(len(dyn), 'and', len(v))

    raw_times_dyn=[]
    for i in range(len(dyn)):
        raw_times_dyn.append(dyn[i,0])

    times_dyn,indices_dyn     = np.unique(raw_times_dyn,return_index=True)
    #print(times_dyn,'\n','===========','\n',indices_dyn)
    #print(len(times_dyn),len(indices_dyn))

    num_of_grains     = indices_dyn[1]- indices_dyn[0]
    #print(num_of_grains)

    iden_first_dyn        =  np.searchsorted(raw_times_dyn,t)
    #print(iden_first_dyn)

    # idée: par example au temps t = 0.3
    #chercher le premier index de t=0.3 dans dyn[] (par searchsorted)
    #Pour chaque t, le data contient de N billes(ici 10k)
    # => Prendre tous ces data de N billes dans le dynt[] 
    for i in range(iden_first_dyn,iden_first_dyn + num_of_grains):
        dynt.append(dyn[i,:])
        #print(dynt[k][:])
        k1=k1+1
    #print(k1)# k should be (num_of_grains to test)
    

    #stock in vt[] : velocities data at time = t of all grains
    raw_times_v=[]
    for i in range(len(v)):
        raw_times_v.append(v[i,0])

    times_v,indices_v     = np.unique(raw_times_v,return_index=True)
    #print(times_v,'\n','===========','\n',indices_v)
    #print(len(times_v),len(indices_v))

    iden_first_v        =  np.searchsorted(raw_times_v,t)
    #print(iden_first_v)

    for i in range(iden_first_v,iden_first_v + num_of_grains):
        vt.append(v[i,:])
        #print(vt[k1][:])
        k2=k2+1
    #print(k2)# k should be (num_of_grains to test)

    #print("-------iden[] of grains at t and between [x1,x2]--------")
    for i in range(len(dynt)):
        if (dynt[i][2] > x1 and dynt[i][2] < x2 and dynt[i][3] > y1 and dynt[i][3] < y2):
            # iden: identity of the grains between [x1,x2] at t
            iden.append(dynt[i][1])
    #assert (len(iden) != 0), "none of grains between [x1,x2] et this time t"

    #print(iden)

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

#mVx,mVy,mMz= vitesse_int_moy_dx_dy()
#print("mVx =", mVx)




### INT MAIN ###
#calculate Velocity average following time of simulation t[] in dx_dy 


dx        = 5
dy        = 1

t1        = 0.8
t2        = t1+10*hstep
t_10steps = np.arange(t1,t2+hstep,hstep) # calculate average in 10 next steps of simulation

flow_height =22
h1_ratio    = np.arange(0,flow_height+1,1) #flow height calculated by flow_height.py
moyenne_Vx= moyenne_Vy= moyenne_Mz=0




print("\n===== Calculate Profil_vitesse by time_averaging [",t1,",",t2,"]========\n")

#caculate Vitesse_int_moy in dx_dy at t[i] and after that, vitesse_moyenne following all t[]
Vx_profil =[]
for y1_ratio in h1_ratio:

    for i in range(len(t_10steps)):

        mVx, mVy, mMz = vitesse_int_moy_dx_dy(x_position, dx, y1_ratio, dy, t_10steps[i])
        if mVx==0 :
            print("no value at t = ",t_10steps[i])
        moyenne_Vx += mVx
        moyenne_Vy += mVy
        moyenne_Mz += mMz
    moyenne_Vx=moyenne_Vx/len(t_10steps)
    moyenne_Vy=moyenne_Vy/len(t_10steps)
    moyenne_Mz=moyenne_Mz/len(t_10steps)

    print("----------velocity_average by time at y1_ratio =",y1_ratio, " ------------------------")
    print("Vx_instant_average = ", moyenne_Vx)
    Vx_profil.append(moyenne_Vx)

print(Vx_profil)




print("\n ========== Calculate Profil_vitesse at t =",t1,"============== \n")

Vx_profil_t1=[]
moyenne_Vx_t1 = moyenne_Vy_t1 = moyenne_Mz_t1=0
for y1_ratio in h1_ratio:

    moyenne_Vx_t1,moyenne_Vy_t1,moyenne_Mz_t1 = vitesse_int_moy_dx_dy(x_position, dx, y1_ratio, dy, t1)
    if moyenne_Vx_t1==0 :
        print("no value at t = ",t1)


    print("----------velocity_average by time at y1_ratio =",y1_ratio, " ------------------------")
    print("Vx_instant_average = ", moyenne_Vx_t1)
    Vx_profil_t1.append(moyenne_Vx_t1)

print(Vx_profil_t1)




plt.figure(1)
plt.plot(Vx_profil,h1_ratio,'o',label = ' profil de Vx')
plt.legend()
plt.title('Profil_Vx at x_pos={} by time-averaging t between [{},{}]'.format(x_position-150,t1,t2))
plt.xlabel('Vx_moyennne')
plt.ylabel('y1_ratio')

plt.figure(2)
plt.plot(Vx_profil_t1,h1_ratio,'o',label = ' profil de Vx')
plt.legend()
plt.title('Profil_Vx at x_pos={} at t1 = {}'.format(x_position-150,t1))
plt.xlabel('Vx_moyennne')
plt.ylabel('y1_ratio')
plt.show();