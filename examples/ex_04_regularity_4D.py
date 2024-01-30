import matplotlib.pyplot as plt
import numpy as np
import nafflib


# Let's look at the regularity of the frequencies with proper identification of the spectral lines
#------------
num_turns = int(1e4)
x_points = np.linspace(1e-6,0.5,50)
px_point = 0.35*x_points

# Henon parameters:
#---
coupling = 0.1
Qx = (3-np.sqrt(5))/2
Qy = np.sqrt(2)-1
#------------


# Choosing number of harmonics
num_harmonics = 20

# Initial tunes (on closed orbit)
Q_vec = [Qx,Qy]


# CORRECT APPROACH:
#======================
# Frequencies to look at:
look_at = { 'x':[   ( 1, 0, 0),
                    ( 3, 0,-1)],
            'y':[   ( 0, 1, 0),
                    ( 1, -1, 0)]}
Q_ref = { 'x':[0,0],
          'y':[0,0]}

fig, axes = plt.subplots(1, 2, figsize=(18, 6))
plt.suptitle(f'Henon map 4D, Regularity. (Qx,Qy) = ({Qx:.3f},{Qy:.3f})')
for x0,px0 in zip(x_points,px_point):
    
    # Tracking
    x,px,y,py = nafflib.henon_map_4D(x0,px0,x0,px0, Qx = Qx,
                                                    Qy = Qy,
                                                    coupling = coupling,
                                                    num_turns= num_turns )
    # Saving in a dict
    dct = {'x':x,'px':px,'y':y,'py':py}

    # Looking at the regularity of the frequencies    
    for plane,ax in zip(['x','y'],axes):
        plt.sca(ax)

        # Extracting harmonics
        A,Q = nafflib.harmonics(dct[plane],dct[f'p{plane}'],num_harmonics,window_order=4)

        # Finding linear combinations based on previous tune (assuming regularity)
        Q_vec = [nafflib.tune(dct[f'{plane}'],dct[f'p{plane}']) for plane in ['x','y']]
        combin,err,freq = nafflib.find_linear_combinations(Q,fundamental_tunes= Q_vec)
        
        J_linear = np.sqrt(x0**2 + px0**2)


        # Plotting
        #----------
        Q_idx = combin.index(look_at[plane][0])
        Q_found = Q[Q_idx]
        if x0 == x_points[0]:
            Q_ref[plane][0] = Q_found
        plt.plot(J_linear,Q_found-Q_ref[plane][0],'.',color='k',alpha=0.9)
        #----------
        Q_idx = combin.index(look_at[plane][1])
        Q_found = Q[Q_idx]
        if x0 == x_points[0]:
            Q_ref[plane][1] = Q_found
        plt.plot(J_linear,Q_found-Q_ref[plane][1],'.',color='C0',alpha=0.9)
        #----------


# Adding labels
for plane,ax in zip(['x','y'],axes):
    plt.sca(ax)
    plt.plot(np.nan,np.nan,'.',color='k' ,alpha=0.9,label=f'Line {look_at[plane][0]}')
    plt.plot(np.nan,np.nan,'.',color='C0',alpha=0.9,label=f'Line {look_at[plane][1]}')
    
    plt.legend()
    
    plt.xlabel(rf'Linear action $J_{plane}$')
    plt.ylabel(rf'$Q - Q_0$')


# Initial tunes (on closed orbit)



# INCORRECT APPROACH:
#======================
Q_vec = [Qx,Qy]
# Frequencies to look at:
look_at = { 'x':[ 0,
                  6],
            'y':[ 0,
                  6]}
Q_ref = { 'x':[0,0],
          'y':[0,0]}

fig, axes = plt.subplots(1, 2, figsize=(18, 6))
plt.suptitle(f'Henon map 4D, Regularity. (Qx,Qy) = ({Qx:.3f},{Qy:3f})')
for x0,px0 in zip(x_points,px_point):
    
    # Tracking
    x,px,y,py = nafflib.henon_map_4D(x0,px0,x0,px0, Qx = Qx,
                                                    Qy = Qy,
                                                    coupling = coupling,
                                                    num_turns= num_turns )
    # Saving in a dict
    dct = {'x':x,'px':px,'y':y,'py':py}

    # Looking at the regularity of the frequencies    
    for plane,ax in zip(['x','y'],axes):
        plt.sca(ax)

        # Extracting harmonics
        A,Q = nafflib.harmonics(dct[plane],dct[f'p{plane}'],num_harmonics,window_order=4)

        # Finding linear combinations based on previous tune (assuming regularity)
        Q_vec = [nafflib.tune(dct[f'{plane}'],dct[f'p{plane}']) for plane in ['x','y']]
        combin,err,freq = nafflib.find_linear_combinations(Q,fundamental_tunes= Q_vec)
        
        J_linear = np.sqrt(x0**2 + px0**2)


        # Plotting
        #----------
        Q_idx = look_at[plane][0]
        Q_found = Q[Q_idx]
        if x0 == x_points[0]:
            Q_ref[plane][0] = Q_found
        plt.plot(J_linear,Q_found-Q_ref[plane][0],'.',color='k',alpha=0.9)
        #----------
        Q_idx = look_at[plane][1]
        Q_found = Q[Q_idx]
        if x0 == x_points[0]:
            Q_ref[plane][1] = Q_found
        plt.plot(J_linear,Q_found-Q_ref[plane][1],'.',color='C0',alpha=0.9)
        #----------


# Adding labels
for plane,ax in zip(['x','y'],axes):
    plt.sca(ax)
    plt.plot(np.nan,np.nan,'.',color='k' ,alpha=0.9,label=f'{look_at[plane][0]}th Line')
    plt.plot(np.nan,np.nan,'.',color='C0',alpha=0.9,label=f'{look_at[plane][1]}th Line')
    
    plt.legend()
    
    plt.xlabel(rf'Linear action $J_{plane}$')
    plt.ylabel(rf'$Q - Q_0$')


plt.show()
        
