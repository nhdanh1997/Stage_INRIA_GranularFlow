import h5py
import numpy as np
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

    """

    # test for flow_height

    ref          = data.get('/data/ref')
    ref_items    = np.array(list(ref.items()))
    print('infos in "data/ref" :',len(ref_items),'\n')
    print(ref_items[1][0],"====")

    sphere         = np.array(ref.get(ref_items[1][0]))
    print("sphere : \n", sphere ,"\n")
    print(f['/data/ref/SphereCS_0_0'].attrs['id'])

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
        #print(iden)

    print(len(id_spheres))   

    print(rayon_spheres)
    print("===============")
    print(len(name_spheres))


    #Other methode to get dataset value
    for key in f.keys():
        # Get the HDF5 group
        group = f[key]

    #group.visit(print) # all the group and subgroup

    # Checkout what keys(dataset) are inside that group.
    for key in group.keys():
        print(key)
    
    solv = group['solv'].value
    static = group['static'].value
    v = group['velocities'].value
    #b_cond = group['boundary_conditions'].value
    cf = group['cf'].value
    dyn = group['dynamic'].value
    #ext_func = group['external_functions'].value

    print("----group-------")
    # (dyn: 1: pas de temps, 2: identifient de sphere, 3 et 4,5:pos x et y, z ; 6:, quaternion([cos(theta/2.), 0. , 0., sin(theta/2.)]
    print(dyn[900])
    
    
    """
    """
    #test for profil_vitesse
    t=0.2
    
    raw_times_dyn=[]
    for i in range(len(dyn)):
        raw_times_dyn.append(dyn[i,0])

    times_dyn,indices_dyn     = np.unique(raw_times_dyn,return_index=True)
    print(times_dyn,'\n','===========','\n',indices_dyn)
    print(len(times_dyn),len(indices_dyn))

    num_of_grains     = indices_dyn[1]- indices_dyn[0]
    print(num_of_grains)

    iden_first_dyn        =  np.searchsorted(raw_times_dyn,t)
    print(iden_first_dyn)

    dynt = []
    vt=[]
    k=0
    for i in range(iden_first_dyn,iden_first_dyn + num_of_grains):
        dynt.append(dyn[i,:])
        #print(dynt[k][:])
        k=k+1
    print(k)# k should be (num_of_grains to test)
    

    k1=0
    raw_times_v=[]
    for i in range(len(v)):
        raw_times_v.append(v[i,0])

    times_v,indices_v     = np.unique(raw_times_v,return_index=True)
    #print(times_v,'\n','===========','\n',indices_v)
    #print(len(times_v),len(indices_v))

    iden_first_v        =  np.searchsorted(raw_times_v,t)
    #print(iden_first_v)

    #k1=0
    for i in range(iden_first_v,iden_first_v + num_of_grains):
        vt.append(v[i,:])
        #print(vt[k1][:])
        k1=k1+1
    #print(k1)# k should be (num_of_grains to test)
    """
    """
    #test for force_on_obstacle
    t=2e-2
    hstep=1e-4
    indices_cf_t=[]
    contact_t=[]
    raw_times_cf=[]

    for i in range(len(cf)):
        raw_times_cf.append(cf[i,0])

    raw_times_cf=np.asarray(raw_times_cf)
    print(raw_times_cf)

    test_indice_cf_t     = np.logical_and(raw_times_cf>t-hstep, raw_times_cf<t+hstep)

    for i in range(len(cf)):
        if(test_indice_cf_t[i]==True):
            contact_t.append(cf[i])

    contact=np.asarray(contact)
    #print(contact)
    #print(len(contact))
    """

    # TEST FOR density
    x_position =150
    t1   =0.3
    t2   =0.4
    dt   =0.1
    t    = np.arange(t1,t2+dt,dt)
    profil_density_all_steps = np.array([[13.929906200308565, 13.176186967062106, 12.397500617206209, 11.68542094661808, 11.01157245610988, 10.39089918661346, 9.914295144285711, 8.920733265226776, 7.922558591334422, 7.324313937013362, 6.665456354695031, 5.963794549691564, 5.144832034384943, 4.346545675152083, 3.5487479461432816, 2.7611928370889802, 2.282302904543015, 1.5023460131287039, 0.7957327923019902, 0.0],[14.667581688919123, 13.922825976654334, 13.333542231638386, 12.70869965190736, 12.058162449871949, 11.318565465375205, 10.467097438174843, 9.721017043437621, 8.774652237916868, 7.941576919563031, 7.093012620737844, 6.346243120745671, 5.683078520904106, 4.993014148407943, 4.152184861430369, 3.3601667582812196, 2.4853100109370603, 1.8709418567267517, 1.1549485197974416, 0.15314629331927734, 0.0]])
    print(len(profil_density_all_steps[0]),len(profil_density_all_steps[1]))
    h=len(profil_density_all_steps[1])

    len_max = 0
    for i in range(len(profil_density_all_steps)):
        if (len(profil_density_all_steps[i])> len_max):
            len_max = len(profil_density_all_steps[i])

    for i in range(len(profil_density_all_steps)):
        while(len(profil_density_all_steps[i])<len_max):
            profil_density_all_steps[i].append(0.)

    profil_density_all_steps = np.array(profil_density_all_steps)
    print(len(profil_density_all_steps[0]),len(profil_density_all_steps[1]))
    print(len(profil_density_all_steps))
    #3D Time evolution of density (fraction volume)
    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d')

    z1 = np.arange(0,h,1)
    y1 = t
    Z1,Y1 = np.meshgrid(z1,y1)
    X1=np.zeros((len(y1),len(z1)))

    print(len(y1),len(z1))

    for i in range(len(y1)):
        X1[i]=profil_density_all_steps[i]

    ax.scatter3D(X1,Y1,Z1, color = "green")
    #ax.view_init(23,-100)
    plt.title("Time evolution of Profil_density  at x_position = {x_pos}".format(x_pos = x_position-150))
    plt.xlabel('Density (m/s)')
    plt.ylabel('t(s)')
    ax.set_zlabel("h (scaled with grain_size)")
    plt.show()

    f.close()