import matplotlib.pyplot as plt
import numpy as np
import nafflib



#========================
#-----
# Henon parameters
x_points  = np.array([0.05,0.1,0.15])
px_point  = 0.35*x_points

# Sample of Q values for different topology
Q_list = [0.2064898024701758,
          0.3761365735491556,
          0.1261960823639152]
#========================
        

# Choosing number of harmonics
num_harmonics = 50

# Choosing map tune
Q0 = Q_list[0]






# Iterating and  PLOTTING
#========================
fig, axes = plt.subplots(1, 2, figsize=(18, 6))
plt.suptitle(f'Henon map, Q = {Q0:.3f}, n_harm = {num_harmonics}')
for x0,px0,color in zip(x_points,px_point,['C0','C3','C4']):

    # Let's generate the map for 1e6 turns
    num_turns = int(1e6)
    x,px = nafflib.henon_map(x0,px0,Q0,num_turns)

    # But exctract the harmonics only for 2e4 turns
    n_naff = int(2e4)
    A,Q = nafflib.harmonics(x[:n_naff],px[:n_naff],num_harmonics,window_order=4)

    
    # And now we compare at the start and at the end of the tracking
    show_turns = 100
    x_s,px_s   = nafflib.generate_signal(A,Q,np.arange(0,show_turns))
    x_e,px_e   = nafflib.generate_signal(A,Q,np.arange(num_turns-show_turns,num_turns))

    err_s = np.abs((x_s-1j*px_s) - (x[:show_turns] - 1j*px[:show_turns]))
    err_e = np.abs((x_e-1j*px_e) - (x[-show_turns:] - 1j*px[-show_turns:]))


    # Plotting
    plt.sca(axes[0])
    plt.plot(x[:show_turns],px[:show_turns],'.',color=color,alpha=0.2)
    plt.plot(x[-show_turns:],px[-show_turns:],'.',color=color,alpha=0.2)

    plt.sca(axes[1])
    plt.plot(np.concatenate((err_s,err_e)),'.',color=color)


plt.sca(axes[1])
plt.axvline(show_turns,color='k',alpha=0.4,label=f'First {show_turns} | Last {show_turns}')
plt.xticks(np.arange(0,2*show_turns+30,30),list(np.arange(0,show_turns,30)) + list(-np.arange(0,show_turns,30))[::-1])
plt.xlabel(f'Turn out of {num_turns:.2e}')
plt.ylabel(rf'Absolute error, $|(x_r - ip_{{x_r}}) - (x - ip_x)|$')
plt.yscale('log')
plt.ylim(1e-16,1e-1)
plt.grid()
plt.legend()


plt.sca(axes[0])
plt.axhline(0,color='k',alpha=0.4)
plt.axvline(0,color='k',alpha=0.4)
plt.axis('square');
plt.xlim([-0.75,0.75])
plt.ylim([-0.75,0.75])
plt.xlabel('x')
plt.ylabel('px')

plt.show()