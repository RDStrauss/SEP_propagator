import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------
data = np.loadtxt("output.txt")

time = data[:,2]
omni_directional_intensity = data[:,3]
anisotropy = data[:,4]
z = data[0,0]
r = data[0,1]

#------------------------------------------------------------------------
#Load data

data = np.loadtxt("Ani.csv", skiprows = 1, delimiter = ',')
time_ani = data[:,0]
ani = data[:,1]


data = np.loadtxt("DI.csv", skiprows = 1, delimiter = ',')
time_di = data[:,0]
di = data[:,1]

#load parameters
param = np.loadtxt("model_setup.txt")

lambdas = param[0]
energy = param[1]
species = param[2]
injection_swtich = param[3]
V_sw = param[4]
acceleration_time = param[5]
escape_time = param[6]
N = param[7]
M = param[8]
Z_max = param[9]
mulimiter = param[10]
zlimiter = param[11]

#normalize model di to observations
omni_directional_intensity = omni_directional_intensity/np.amax(omni_directional_intensity)*np.amax(di)

#Plot the figures
fig = plt.figure(figsize=(5,6))

subplot2 = fig.add_subplot(211)
#subplot2.set_xlabel('Time (hr)')
#subplot2.set_title('aap')
subplot2.set_ylabel('Intensity (arb. units)')
subplot2.set_yscale('log')
subplot2.set_xscale('linear')
subplot2.set_ylim([0.005,3.])
#subplot2.ticklabel_format(useOffset=False)
subplot2.set_xlim([-0.5,12])


subplot2.plot(time,omni_directional_intensity, color = 'red')
subplot2.scatter(time_di, di, color = 'blue', facecolor = 'none')
subplot2.set_xticklabels([])

# add labels


subplot3 = fig.add_subplot(413)  
subplot3.set_xlabel('Time')
subplot3.set_ylabel('Anisotropy')
subplot3.set_yscale('linear')
subplot3.set_xscale('linear')
subplot3.set_ylim([-0.2,3.0])
#subplot3.ticklabel_format(useOffset=False)
subplot3.set_xlim([-0.5,12])
  
subplot3.plot(time, anisotropy, color = 'red')
subplot3.scatter(time_ani, ani, color = 'blue', facecolor = 'none')
subplot3.axhline(y = 0, color = 'black', linestyle = '--')

plt.subplots_adjust(hspace=0.15)
fig.align_ylabels([subplot2,subplot3])
fig.tight_layout

fig.savefig('Output.png', dpi=300, bbox_inches='tight')

#--------------------------------------------------------------------

plt.show()
