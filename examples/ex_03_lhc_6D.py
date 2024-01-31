import matplotlib.pyplot as plt
import numpy as np
import nafflib



#=================================
num_turns = int(3e4)
files = [nafflib.__path__[0] + f'/../tests/data/LHC_particle_{s}_momentum_{i}sigma.csv'  for s in ['on','off'] for i in [1,3,5]]
#-------------------
def read_csv(filename):
    filecontent = np.genfromtxt(filename,
                            delimiter   =',',
                            skip_header = 1,
                            converters  = {col: lambda s: complex(s.decode()) for col in [1,3,5]},
                            unpack=True)
    data = {}
    data['Ax']      = filecontent[1]
    data['Qx']      = filecontent[2]
    data['Ay']      = filecontent[3]
    data['Qy']      = filecontent[4]
    data['Azeta']   = filecontent[5]
    data['Qzeta']   = filecontent[6]
    return data
#-------------------
example_signals = []
for file in files:
    data = read_csv(file)
    data['set'] = ['on_momentum' if 'on_momentum' in file else 'off_momentum'][0]
    for plane in ['x','y','zeta']:
        z,pz = nafflib.generate_signal(data[f'A{plane}'],data[f'Q{plane}'],np.arange(num_turns))
        data[f'{plane}']  = z
        data[f'p{plane}'] = pz
    example_signals.append(data)
#=================================
    



# Choosing number of harmonics and particls set
particle_set = ['on_momentum','off_momentum'][1]
num_harmonics = 50

# Choosing turns to plot
N_start = 0
N_stop  = 100

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
plt.suptitle(f'LHC, Particles {particle_set}')
N = np.arange(num_turns)
for data in example_signals:
    # Choosing particles set
    if data['set'] != particle_set:
        continue

    # Extracting fundamental frequencies
    Q_vec = [nafflib.tune(data[f'{plane}'],data[f'p{plane}'],window_order=4) for plane in ['x','y','zeta']]
    
    # Plotting plane-by-plane
    for plane,ax in zip(['x','y','zeta'],axes):
        plt.sca(ax)
        A,Q      = data[f'A{plane}'],data[f'Q{plane}']
        z_r,pz_r = data[f'{plane}'] ,data[f'p{plane}']

        # Finding harmonics
        A_r,Q_r = nafflib.harmonics(z_r,pz_r,num_harmonics, window_order=4)

        # Reconstructing signal (original is already a reconstruction...)
        z_rr,pz_rr = nafflib.generate_signal(A_r,Q_r,N)

        # Plotting
        plt.plot(z_r[N_start:N_stop] ,pz_r[N_start:N_stop],'o',mfc='none',color='C4' ,alpha=0.9)
        plt.plot(z_rr[N_start:N_stop],pz_rr[N_start:N_stop],'.'   ,color='C0' ,alpha=0.9)

        # If one wants to see the harmonics:
        combin,err,freq = nafflib.find_linear_combinations(Q_r,fundamental_tunes= Q_vec)
        # print(combin)

        
# Adding labels
for plane,ax in zip(['x','y','zeta'],axes):
    plt.sca(ax)

    plt.plot(np.nan,np.nan,'o',mfc='none',color='C4' ,alpha=0.9,label='LHC trk (100 harmonics)')
    plt.plot(np.nan,np.nan,'.',color='C0' ,alpha=0.9,label=f'nafflib ({num_harmonics} harmonics)')
    plt.legend()

    plt.axis('equal')
    if plane == 'zeta':
        plane = '\zeta'
    plt.xlabel(rf'$\tilde {plane}/\sqrt{{\varepsilon_{plane}}}$')
    plt.ylabel(rf'$\tilde p_{plane}/\sqrt{{\varepsilon_{plane}}}$')
    
    if plane == '\zeta':
        plt.xlim(-0.5,0.5)
        plt.ylim(-0.5,0.5)
    else:
        plt.xlim(-7,7)
        plt.ylim(-7,7)

plt.show()