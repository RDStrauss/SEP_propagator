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

#add observational background
omni_directional_intensity = omni_directional_intensity + np.amin(di)

#Plot the figures
fig = plt.figure(figsize=(5,6))

#f, (subplot1, subplot2) = plt.subplots(2, 1, sharex=True,figsize=(5,10))

subplot1 = plt.subplot2grid((5, 1), (0, 0), rowspan=1, colspan=1)
subplot2 = plt.subplot2grid((5, 1), (1, 0), rowspan=3, colspan=1)
subplot3 = plt.subplot2grid((5, 1), (4, 0), rowspan=1, colspan=1)

#time = np.linspace(start = 0, stop = 10, num = 200)

#injection profile
ra_profile = 0.25/time*np.exp(-acceleration_time/time - time/escape_time)
ra_profile = ra_profile/np.amax(ra_profile)

subplot1.set_title('2010 Feb 07 STEREO-B/SEPT 65- 105 keV Electrons')
subplot1.set_ylabel('Injection Profile')
subplot1.set_yscale('linear')
subplot1.set_xscale('linear')
#subplot1.set_ylim([50.,10000.])
#subplot1.ticklabel_format(useOffset=False)
subplot1.set_xlim([-0.5,10])


subplot1.plot(time,ra_profile, color = 'red')
subplot1.set_xticklabels([])

#subplot2 = fig.add_subplot(211)
#subplot2.set_xlabel('Time (hr)')
#subplot2.set_title('aap')
subplot2.set_ylabel('Intensity (cm$^2$ sr s MeV)$^{-1}$')
subplot2.set_yscale('log')
subplot2.set_xscale('linear')
subplot2.set_ylim([50.,10000.])
#subplot2.ticklabel_format(useOffset=False)
subplot2.set_xlim([-0.5,10])


subplot2.plot(time,omni_directional_intensity, color = 'red')
subplot2.scatter(time_di - 2.5, di, color = 'blue', facecolor = 'none')
subplot2.set_xticklabels([])

# add labels
V_sw = V_sw*1.5e8/60./60.
s = '$\lambda$ = ' + str(round(lambdas,3)) + ' AU;' + ' E = ' + str(round(energy,3)) + ' MeV'
plt.text(2.0, 5.5, s)
s = '$V_{sw}$ = ' + str(round(V_sw,3)) + ' km/s'
plt.text(2.0, 5.0, s)
s = 'particle = ' + str(int(species)) 
plt.text(2.0, 4.5, s)
s = '(1 = electrons; 2 = protons)'
plt.text(2.0, 4.0, s)
s = '$t_a$, $t_e$ = ' + str(round(acceleration_time,1)) + ',' + str(round(escape_time,2)) + ' hrs'  
plt.text(2.0, 3.5, s)

#subplot3 = fig.add_subplot(313,sharex=True)  
subplot3.set_xlabel('Event time - 2.5 (hr)')
subplot3.set_ylabel('Anisotropy')
subplot3.set_yscale('linear')
subplot3.set_xscale('linear')
subplot3.set_ylim([-0.2,2.2])
#subplot3.ticklabel_format(useOffset=False)
subplot3.set_xlim([-0.5,10])
  
subplot3.plot(time, anisotropy, color = 'red')
subplot3.scatter(time_ani - 2.5, ani, color = 'blue', facecolor = 'none')
subplot3.axhline(y = 0, color = 'black', linestyle = '--')

plt.subplots_adjust(hspace=0.15)
fig.align_ylabels([subplot1,subplot2,subplot3])
fig.tight_layout

fig.savefig('Output.pdf', dpi=300, bbox_inches='tight')

#--------------------------------------------------------------------

plt.show()
