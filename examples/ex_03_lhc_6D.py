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
    data['name'] = file
    for plane in ['x','y','zeta']:
        z,pz = nafflib.generate_signal(data[f'A{plane}'],data[f'Q{plane}'],np.arange(num_turns))
        data[f'{plane}']  = z
        data[f'p{plane}'] = pz
    example_signals.append(data)
#=================================
    


# TODO


# fig, axes = plt.subplots(1, 3, figsize=(18, 6))
# plt.suptitle('Particles on momentum')
# for part in particles_on_p:
#     for plane,ax in zip(['x','y','zeta'],axes):
#         plt.sca(ax)
#         plt.plot(part[plane],part[f'p{plane}'],'.',color='k',alpha=0.1)

# # Adding labels
# for plane,ax in zip(['x','y','zeta'],axes):
#     plt.sca(ax)
#     plt.axis('equal')
#     if plane == 'zeta':
#         plane = '\zeta'
#     plt.xlabel(rf'$\tilde {plane}/\sqrt{{\varepsilon_{plane}}}$')
#     plt.ylabel(rf'$\tilde p_{plane}/\sqrt{{\varepsilon_{plane}}}$')
    
#     if plane == '\zeta':
#         plt.xlim(-0.5,0.5)
#         plt.ylim(-0.5,0.5)
#     else:
#         plt.xlim(-6,6)
#         plt.ylim(-6,6)

# fig, axes = plt.subplots(1, 3, figsize=(18, 6))
# plt.suptitle('Particles off momentum')
# for part in particles_off_p:
#     for plane,ax in zip(['x','y','zeta'],axes):
#         plt.sca(ax)
#         plt.plot(part[plane],part[f'p{plane}'],'.',color='k',alpha=0.1)

# # Adding labels
# for plane,ax in zip(['x','y','zeta'],axes):
#     plt.sca(ax)
#     plt.axis('equal')
#     if plane == 'zeta':
#         plane = '\zeta'
#     plt.xlabel(rf'$\tilde {plane}/\sqrt{{\varepsilon_{plane}}}$')
#     plt.ylabel(rf'$\tilde p_{plane}/\sqrt{{\varepsilon_{plane}}}$')
    
#     if plane == '\zeta':
#         plt.xlim(-0.5,0.5)
#         plt.ylim(-0.5,0.5)
#     else:
#         plt.xlim(-6,6)
#         plt.ylim(-6,6)
