import h5py
import numpy as np
from numpy import linalg
import math
import matplotlib.pyplot as plt

filename = "wall_test.hdf5"

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

grain_size = 0.001
ground_thickness = 10 * grain_size
ground_y_shift = 0.1
y0 = - ground_thickness - ground_y_shift

wall_x_position = 250  # ATTENTION: TANK SIZE =150 => 200 ~ 50 FROM GATE OF TANK)
wall_height = 13
wall_thickness = 10
hstep = 1e-4


def dead_zone(V_thresold, t=0.5, dy=1,dx =1):
    """
    This function detect the grains in dead_zone
    t : time
    dy: to test with cas different V_thresold follow height
    dx: step for detection dead_zone_line 
    
    ideal:
    follow the wall_height with step = dy
    if v < vt => grain quasi-statique situate in dead_zone
       v > vt => grain out of range of dead zone
       vt for difference height position y alongside the wall
    vt(thresold velocity) = 0.05*u with u : the depth-averaged velocity (bagnol profil)     

    return  grains_dead_zone_dyn : dynamics of grains in dead zone
            grains_dead_zone_line : coordination of grains in dead zone line
            L_dead_zone: length of dead zone
            S_dead_zone: area of dead zone

    """
    height = np.arange(0, len(V_thresold)+1, 1)

    dynt = []  # sotck data of dynamics a time = t
    raw_times_dyn = []  # stock all time steps  t[] in dynt[]

    vt = []  # sotck data of velocities a time = t
    raw_times_v = []  # stock all time steps  t[] in dynt[]

    grains_dead_zone_dyn = []  # stocks dynamics grains detected in dead_zone

    L_dead_zone = 0
    S_dead_zone = 0

    # Build dynt[]
    for i in range(len(dyn)):
        raw_times_dyn.append(dyn[i, 0])
    raw_times_dyn = np.asarray(raw_times_dyn)

    test_indice_dyn_t = np.logical_and(raw_times_dyn > t-hstep, raw_times_dyn < t+hstep)
    # return matrix that element be true for all t satisfy
    # print(test_indice_dyn_t)

    for i in range(len(dyn)):
        if(test_indice_dyn_t[i] == True):
            dynt.append(dyn[i, :])
            # print(dynt[i,:])
    dynt = np.array(dynt)

    # Build vt[]
    for i in range(len(v)):
        raw_times_v.append(v[i, 0])
    raw_times_v = np.asarray(raw_times_v)

    test_indice_v_t = np.logical_and(
        raw_times_v > t-hstep, raw_times_v < t+hstep)
    # return matrix that element be true for all t satisfy
    # print(test_indice_dyn_t)

    # Build grains_dead_zone_dyn
    for i in range(len(v)):
        if(test_indice_v_t[i] == True):
            vt.append(v[i, :])
    vt = np.array(vt)

    # Build grains_dead_zone_dyn
    for i in range(len(vt)):
        if(vt[i, 2] < (sum(V_thresold)/len(V_thresold)) and dynt[i, 2] < 200*grain_size):
            #(i,'====', vt[i,2])
            grains_dead_zone_dyn.append(dynt[i, :])

    grains_dead_zone_dyn = np.asarray(grains_dead_zone_dyn)
    print(len(grains_dead_zone_dyn))
    """
    """
    # Calculate L_dead_zone
    x_grains_dead_zone = []
    for i in range(len(grains_dead_zone_dyn)):
        x_grains_dead_zone.append(grains_dead_zone_dyn[i][2])
    x_grains_dead_zone = np.array(x_grains_dead_zone)
    L_dead_zone = max(x_grains_dead_zone) - min(x_grains_dead_zone)

    L_dead_zone_grain_size_scaled = L_dead_zone/grain_size

    x_ratio_dead_zone = np.arange(0, L_dead_zone_grain_size_scaled+1, 1)
    print(x_ratio_dead_zone)

     #Caculate S_dead_zone by building grains_dead_zone_line
    #Build grains_dead_zone_line
    
    grains_dead_zone_line = []  # stocks all grains at tallest position in each pas dx

    for j in range(len(x_ratio_dead_zone)-dx):
        x1 = (wall_x_position -
              x_ratio_dead_zone[j]-wall_thickness*0.5)*grain_size
        x2 = (wall_x_position -
              x_ratio_dead_zone[j+dx]-wall_thickness*0.5)*grain_size
        grains_in_dx = []  # stocks grains in each pas dx to find the grain at tallest position
        k = 0
        for i in range(len(grains_dead_zone_dyn)):
            if(grains_dead_zone_dyn[i, 2] < x1 and grains_dead_zone_dyn[i, 2] > x2):
                grains_in_dx.append(grains_dead_zone_dyn[i, :])
                k += 1
                print(j, "====", k)
        grains_in_dx = np.array(grains_in_dx)
        print(len(grains_in_dx))

        if(len(grains_in_dx) != 0):
            y_grains_in_dx = grains_in_dx[:, 3]
            tallest_grain = max(y_grains_in_dx)
            grains_dead_zone_line.append([x1, tallest_grain])
    grains_dead_zone_line = np.array(grains_dead_zone_line)
    print(grains_dead_zone_line)

    #Caculate S_dead_zone
    for i in range(len(grains_dead_zone_line)-1):
        h = grains_dead_zone_line[i+1,0]- grains_dead_zone_line[i,0]
        # print(h)
        S_trapeze = (grains_dead_zone_line[i+1,1]+grains_dead_zone_line[i,1])*h
        # print(S_parallepipede)
        S_dead_zone += S_trapeze
    S_dead_zone = 0.5*S_dead_zone
    #print('S_dead_zone =',S_dead_zone )

    return grains_dead_zone_dyn, grains_dead_zone_line, L_dead_zone_grain_size_scaled, S_dead_zone


# Thresold velocity
# t=0.2
V_thresold_t01 = 0.05*np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.8957518759895774, 0.9329946156471006, 0.9586240958284449, 0.9831273763076119, 1.0013674905425625, 1.0069515426953635, 1.0206746784123508, 1.0391508681433541, 1.0709348917007446, 0, 0, 0], [0.629127667138451, 0.6572620298181261, 0.7067438066005707, 0.7342914268374443, 0.7445247230075654, 0.7628892239402322, 0.783287878219898, 0.7859816133975983, 0.7893734931945801, 0.8517604470252991, 0, 0], [0.5489887081478771, 0.5847040292094735, 0.6141084638135187, 0.6424858694275221, 0.6604870438575745, 0.6764473207294941, 0.6940112809340159, 0.723680283342089, 0.7787586152553558, 0, 0, 0], [0.5190331492573023, 0.5655379210199628, 0.5819473601877689, 0.6079402852941442, 0.6340156594912211, 0.6492217928171158, 0.6607640060511503, 0.6746038283620562, 0.7189599474271139, 0.7412778735160828, 0, 0], [0.4241938606968948, 0.4793618632399518, 0.4951381772756577, 0.5239948984235525, 0.5471581816673279, 0.5808665454387665, 0.6143896381060282, 0.6428987383842468, 0, 0, 0, 0], [0.3632000440900976, 0.3980704645315806, 0.42801018861623913, 0.4306018021371629, 0.3754631181557973, 0, 0, 0, 0, 0, 0, 0], [0.37807474136352537, 0.43901134350083093, 0.48347866038481396, 0.5031113624572754, 0.518687903881073, 0, 0, 0, 0, 0, 0, 0]])
V_thresold = V_thresold_t01
print(V_thresold)



# print(len(V_thresold07))
dy = 1
dx = 1


dt=0.1
t=np.arange(0.1,(len(V_thresold)+1)*dt,dt)

L_dead_zone_all_steps = []
S_dead_zone_all_steps = []

for i in range(len(t)):
    grains_dead_zone_dyn, grains_dead_zone_line, L_dead_zone, S_dead_zone = dead_zone(V_thresold[i,:], t[i], dy,dx)
    L_dead_zone_all_steps.append(L_dead_zone)
    S_dead_zone_all_steps.append(S_dead_zone)

"""
print(grains_dead_zone_line)
print(grains_dead_zone_dyn)


print("L_dead_zone (grain_size scaled) = ", L_dead_zone_grain_size_scaled)
print("S_dead_zone = ", S_dead_zone)

x_dead_zone = grains_dead_zone_dyn[:, 2] / grain_size
y_dead_zone = (grains_dead_zone_dyn[:, 3]-y0)/grain_size


plt.figure(1)
plt.plot(x_dead_zone, y_dead_zone, 'o', color='green')
plt.legend()
plt.title('dead zone at t = {}'.format(t))
plt.xlabel('x (grain_size scale)')
plt.ylabel('y (grain_size scale)')

x_grains_dead_zone_line = grains_dead_zone_line[:, 0]/grain_size
y_grains_dead_zone_line = (grains_dead_zone_line[:, 1]-y0)/grain_size

plt.figure(2)
plt.plot(x_grains_dead_zone_line, y_grains_dead_zone_line, 'o', color='green')
plt.legend()
plt.title('dead_zone_line at t = {}'.format(t))
plt.xlabel('x (grain_size scale)')
plt.ylabel('y (grain_size scale)')
"""

plt.figure(1)
plt.plot(t, L_dead_zone_all_steps, 'o', color='green')
plt.legend()
plt.title('Time evolution of L_dead_zone with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('L_dead_zone (scaled zith grain_size)')

plt.figure(2)
plt.plot(t, S_dead_zone_all_steps, 'o', color='green')
plt.legend()
plt.title('Time evolution of S_dead_zone with dt = {}'.format(dt))
plt.xlabel('t (s)')
plt.ylabel('S_dead_zone ')


plt.show()
