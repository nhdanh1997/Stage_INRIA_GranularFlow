import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from statistics import mode
filename = "2d_sphere_flow_N-3000-e-0.5-mu-0.5-angle-30.hdf5"
print(filename,"\n")


with h5py.File(filename, "r") as f:

    ls    = list(f.keys())
    #print('list of group : \n',ls,'\n')
    group =  np.array(f.get('data'))
    #print('list of datasets : \n', group,'\n')
    
    base_items = list(f.items())
    #print('infos of items in the base director  : \n', base_items,'\n')
    data       = f.get('data')
    data_items = list(data.items())
    #print('infos of data in "data" :\n ',data_items,'\n')
    
    dyn = np.array(data.get('dynamic'))
    v   = np.array(data.get('velocities'))
    cf  = np.array(data.get('cf'))


    f.close()

grain_size       = 0.001
ground_thickness = 10 * grain_size
ground_y_shift   = 0.09
y0               = - ground_thickness - ground_y_shift

hstep=1e-4

def vitesse_int_moy_dx_dy(x1_ratio=50, dx=5, y1_ratio=0, dy=1, t=0.25):
    """
    #calcul instantaneous average velocity of grains between x1 and x2 at time t	
    
    x,y_ratio : scale for grain_size meanwhile x,y :real value position
    x1_ratio= 200 (150(tank) + 50)
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
    

#mVx,mVy,mMz= vitesse_int_moy_dx_dy()
#print("mVx =", mVx)

### INT MAIN ###
#calculate Velocity average following time of simulation t[] in dx_dy 


num_of_grains = 3000

n_row = num_of_grains/100
n_col = 100
distribution      = ('uniform', 0.15) 
tank_width_number = n_col+2
tank_width        = tank_width_number*grain_size*(1.0+distribution[1])

x_position_from_gate = 300

#time of simulation
t_ini= 1.
t_end= 3.
dt   = 0.5
t1   = np.arange(t_ini,t_end+dt,dt)

x_position = tank_width/grain_size + x_position_from_gate 
dx = 5
dy = 1

#flow_height time-averaged and instantaneous calculated by flow_height.py
h_time_average =[2.0, 11.425742574257425, 15.794117647058824, 14.138613861386139, 16.762376237623762]
h              =[2, 12, 16, 14, 17]
h1_ratio   = np.arange(0,max(h)+1,1) #flow height calculated by flow_height.py

"""
"""
################ 1ere SOLUTION #########################
# time evolution of profil_vitesse average in next 10 step at each moment t
print("\n=====  Profil_vitesse time-average========\n")

moyenne_Vx = moyenne_Vy= moyenne_Mz=0
Vx_profil_10steps_all_steps=[]
#caculate Vitesse_int_moy in dx_dy at t[i] and after that, vitesse_moyenne following all t[]
for t in range(len(t1)):
    Vx_profil = []
    t_10steps = np.arange(t1[t]-100*hstep,t1[t]+hstep,hstep) # calculate average in 10 next steps of simulation
    for y1_ratio in h1_ratio:

        for i in range(len(t_10steps)):

            mVx, mVy, mMz = vitesse_int_moy_dx_dy(x_position, dx, y1_ratio, dy, t_10steps[i])
            #if mVx!=0 :
                #print("no value at t = ",t_10steps[i])
            moyenne_Vx += mVx
            moyenne_Vy += mVy
            moyenne_Mz += mMz
        moyenne_Vx=moyenne_Vx/len(t_10steps)
        moyenne_Vy=moyenne_Vy/len(t_10steps)
        moyenne_Mz=moyenne_Mz/len(t_10steps)

        #print("----------velocity_average by time at y1_ratio =",y1_ratio, " ------------------------")
        #print("Vx_instant_average = ", moyenne_Vx)
        Vx_profil.append(moyenne_Vx)
    Vx_profil_10steps_all_steps.append(Vx_profil) 

print(Vx_profil_10steps_all_steps) # matrix 2d (n_row = len(t1),n_col = len(h1_ratio) 
Vx_profil_10steps_all_steps = np.array(Vx_profil_10steps_all_steps)

#calculate time evolution of fow_rate with average 10next steps each t
time_evolution_flow_rate_10steps = []
for i in range(len(Vx_profil_10steps_all_steps)):
    flow_rate = sum(Vx_profil_10steps_all_steps[i])*h_time_average[i]*grain_size/len(Vx_profil_10steps_all_steps[i])
    time_evolution_flow_rate_10steps.append(flow_rate)

print("==flow_rate time averaged == ")
print(time_evolution_flow_rate_10steps)


# Find t when flow_rate become constant    
t_all_10steps = np.arange(0,len(time_evolution_flow_rate_10steps),1)
print(len(t_all_10steps),len(time_evolution_flow_rate_10steps))

counter_flow_rate_10steps    = Counter(time_evolution_flow_rate_10steps)
const_flow_rate_10steps      = mode(time_evolution_flow_rate_10steps) # value at this flow_rate become constant
print("flow_rate_10steps_constant = ",const_flow_rate_10steps)
print("num of value flow_rate_constant = ",counter_flow_rate_10steps[const_flow_rate_10steps])
iden_flow_rate_10steps_const = np.searchsorted(time_evolution_flow_rate_10steps,const_flow_rate_10steps)
#print(iden_flow_rate_10steps_const)
T_const_flow_rate_10steps    = t_all_10steps[iden_flow_rate_10steps_const]
print("T_const =",T_const_flow_rate_10steps )



################ 2em SOLUTION #########################
#time evolution of profil_vitesse at each moment t
print("\n =======  Profil_vitesse instantaneous ======= \n")

Vx_profil_t1_all_steps=[]
moyenne_Vx_t1 = moyenne_Vy_t1 = moyenne_Mz_t1=0

for i in range(len(t1)):
    #print("\n value at time t = ", t1[i])
    Vx_profil_t1=[]
    for y1_ratio in h1_ratio:

        moyenne_Vx_t1,moyenne_Vy_t1,moyenne_Mz_t1 = vitesse_int_moy_dx_dy(x_position, dx, y1_ratio, dy, t1[i])
        #if moyenne_Vx_t1!=0 :
            #print("value at t = ",t1[i],", h = ",y1_ratio)

        #print("----------velocity_average by time at y1_ratio =",y1_ratio, " ------------------------")
        #print("Vx_instant_average = ", moyenne_Vx_t1)
        Vx_profil_t1.append(moyenne_Vx_t1)
    #print(Vx_profil_t1)
    Vx_profil_t1_all_steps.append(Vx_profil_t1)    

print(Vx_profil_t1_all_steps) # matrix 2d (n_row = len(t1),n_col = len(h1_ratio) 
Vx_profil_t1_all_steps=np.array(Vx_profil_t1_all_steps)

#caclutate time evolution of flow_rate each t
time_evolution_flow_rate = []
for i in range(len(Vx_profil_t1_all_steps)):
    flow_rate = sum(Vx_profil_t1_all_steps[i])*h[i]*grain_size/len(Vx_profil_t1_all_steps[i])
    time_evolution_flow_rate.append(flow_rate)

print("== flow_rate== ")
print(time_evolution_flow_rate)


# Find t when flow_rate become constant    
t_all = np.arange(0,len(time_evolution_flow_rate),1)
print(len(t_all),len(time_evolution_flow_rate))

counter_flow_rate    = Counter(time_evolution_flow_rate)
const_flow_rate      = mode(time_evolution_flow_rate) # value at this flow_rate become constant
print("flow_rate_constant = ",time_evolution_flow_rate)
print("num of value flow_rate_constant = ",counter_flow_rate[const_flow_rate])
iden_flow_rate_const = np.searchsorted(time_evolution_flow_rate,const_flow_rate)
#print(iden_flow_rate_const)
T_const_flow_rate    = t_all[iden_flow_rate_const]
print("T_const =",T_const_flow_rate )


#3D time evolution of profil_vitesse average in next 10 step at each moment t
"""
"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

z1    = h1_ratio
y1    = t1
Z1,Y1 = np.meshgrid(z1,y1)
X1    = np.zeros((len(y1),len(z1)))

for i in range(len(y1)):
    X1[i]=Vx_profil_10steps_all_steps[i,:]

ax.scatter3D(X1,Y1,Z1, color = "green")
#ax.view_init(23,-100)
plt.suptitle('Evolution of Profil_Vx time-averaged at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
plt.title(filename, fontsize=10)
plt.xlabel('Vx_time-averaged (m/s)')
plt.ylabel('t(s)')
ax.set_zlabel("y/d")



#3D time evolution of profil_vitesse at each moment t
fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')

z2    = h1_ratio
y2    = t1
Z2,Y2 = np.meshgrid(z2,y2)
X2    =np.zeros((len(y2),len(z2)))
 
#print(len(y2))
for i in range(len(y2)):
    X2[i]=Vx_profil_t1_all_steps[i,:]

ax.scatter3D(X2,Y2,Z2, color = "green")
#ax.set_xlim(0, 23)
#ax.set_zlim(0, 2)
#ax.view_init(1,-100)
plt.suptitle('Evolution of Profil_Vx at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
plt.title(filename, fontsize=10)
plt.xlabel('Vx (m/s)')
plt.ylabel('t(s)')
ax.set_zlabel("y/d")

"""
"""
#graph results of time evolution of flow_rate
width = dt/4
plt.figure()
plt.bar(t1-width/2,time_evolution_flow_rate_10steps,width=width,label='flow_rate_time-average')
plt.bar(t1+width/2,time_evolution_flow_rate,width=width,label='flow_rate')
plt.legend()
plt.suptitle('Time evolution of Flow_rate at x ={x_pos}, dt = {dt}'.format(x_pos= x_position_from_gate,dt=dt),fontsize=15)
plt.title(filename, fontsize=10)
plt.xlabel('time t(s)')
plt.ylabel('m3.s-1')
plt.show();

