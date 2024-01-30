import matplotlib.pyplot as plt
import numpy as np
import nafflib



#========================
#-----
# Henon parameters
num_turns = int(1e5)
x_points  = np.array([0.1,0.3,0.45,0.5])
px_point  = 0.35*x_points

# Sample of Q values for different topology
Q_list = [0.2064898024701758,
          0.3761365735491556,
          0.1261960823639152]
#========================
        

# Choosing number of harmonics
num_harmonics = 50

# Choosing map tune
Q0 = Q_list[2]

# Choosing turns to plot
N_start = 0
N_stop  = 100



# Iterating and  PLOTTING
#========================
plt.figure(figsize=(10,10))
plt.suptitle(f'Henon map, Q = {Q0}')
for x0,px0 in zip(x_points,px_point):
    x,px = nafflib.henon_map(x0,px0,Q0,num_turns)

    # Extracting harmonics
    A,Q = nafflib.harmonics(x,px,num_harmonics,window_order=4)

    # Reconstructing signal
    x_r,px_r = nafflib.generate_signal(A,Q,np.arange(num_turns))

    plt.plot(x[N_start:N_stop],px[N_start:N_stop],'o',mfc='none',color='C4' ,alpha=0.9)
    plt.plot(x_r[N_start:N_stop],px_r[N_start:N_stop],'.'   ,color='C0' ,alpha=0.9)

plt.plot(np.nan,np.nan,'o',mfc='none',color='C4' ,alpha=0.9,label='Henon')
plt.plot(np.nan,np.nan,'.',color='C0' ,alpha=0.9,label='nafflib')

plt.axhline(0,color='k',alpha=0.4)
plt.axvline(0,color='k',alpha=0.4)
plt.axis('square');
plt.xlim([-0.75,0.75])
plt.ylim([-0.75,0.75])
plt.xlabel('x')
plt.ylabel('px')
plt.legend()
plt.show()