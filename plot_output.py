import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------
data = np.loadtxt("output.txt")

time = data[:,2]
omni_directional_intensity = data[:,3]
anisotropy = data[:,4]
z = data[0,0]
r = data[0,1]

times = np.zeros(5)
pitch_angle_distribution = np.zeros((5,99))
data = np.loadtxt("pitch_angle_distribution.txt")
mu = data[:,0]
times[0] = data[0,1]
times[1] = data[0,3]
times[2] = data[0,5]
times[3] = data[0,7]
times[4] = data[0,9]
pitch_angle_distribution[0,:] = data[:,2]
pitch_angle_distribution[1,:] = data[:,4]
pitch_angle_distribution[2,:] = data[:,6]
pitch_angle_distribution[3,:] = data[:,8]
pitch_angle_distribution[4,:] = data[:,10]

fig = plt.figure()

subplot1 = fig.add_subplot(221)
#subplot1.set_xlabel('Time (hr)')
#subplot1.set_title('aap')
subplot1.set_ylabel('Omni-directional intensity')
subplot1.set_yscale('log')
subplot1.set_xscale('linear')
#subplot1.set_ylim([0.001,10.])
#subplot1.ticklabel_format(useOffset=False)
#subplot1.set_xlim([-0.5,10])

subplot1.plot(time,omni_directional_intensity, color = 'red')

subplot1 = fig.add_subplot(223)
  
subplot1.set_xlabel('Time (hr)')
subplot1.set_ylabel('Anisotropy')
subplot1.set_yscale('linear')
subplot1.set_xscale('linear')
#subplot1.set_ylim([0.001,10.])
#subplot1.ticklabel_format(useOffset=False)
#subplot1.set_xlim([-0.5,10])
  
subplot1.plot(time, anisotropy, color = 'blue')
subplot1.axhline(y = 0, color = 'black', linestyle = '--')

subplot1 = fig.add_subplot(222)

#subplot1.set_title('aap')
subplot1.set_xlabel('$\mu$')
subplot1.set_ylabel('Distribution function')
subplot1.yaxis.tick_right()
subplot1.yaxis.set_ticks_position('both')
subplot1.yaxis.set_label_position("right")
subplot1.set_yscale('linear')
subplot1.set_xscale('linear')
#subplot1.set_ylim([0.001,10.])
#subplot1.ticklabel_format(useOffset=False)
subplot1.set_xlim([-1,1])

s = 't = ' + str(round(times[0],1)) + ' hr'
subplot1.plot(mu, pitch_angle_distribution[0,:]/np.max(pitch_angle_distribution[0,:]), color = 'blue', linestyle = '-', label = s)

s = 't = ' + str(round(times[1],1)) + ' hr'
subplot1.plot(mu, pitch_angle_distribution[1,:]/np.max(pitch_angle_distribution[1,:]), color = 'blue', linestyle = '--', label = s)

s = 't = ' + str(round(times[2],1)) + ' hr'
subplot1.plot(mu, pitch_angle_distribution[2,:]/np.max(pitch_angle_distribution[2,:]), color = 'red', linestyle = '-', label = s)

s = 't = ' + str(round(times[3],1)) + ' hr'
subplot1.plot(mu, pitch_angle_distribution[3,:]/np.max(pitch_angle_distribution[3,:]), color = 'red', linestyle = '--', label = s)

s = 't = ' + str(round(times[4],1)) + ' hr'
subplot1.plot(mu, pitch_angle_distribution[4,:]/np.max(pitch_angle_distribution[4,:]), color = 'green', linestyle = '-', label = s)
plt.legend(fontsize = 8)

subplot1 = fig.add_subplot(224)
subplot1.set_xlim([-1,1])
subplot1.set_ylim([-1,1])
subplot1.axis('off')

s = 'z = ' + str(round(z,2)) + ' AU;' + ' r = ' + str(round(r,2)) + ' AU'
plt.text(-0.9, 0.5, s)


data = np.loadtxt("model_setup.txt")

lambdas = data[0]
energy = data[1]
species = data[2]
injection_swtich = data[3]
V_sw = data[4]
acceleration_time = data[5]
escape_time = data[6]
N = data[7]
M = data[8]
Z_max = data[9]
mulimiter = data[10]
zlimiter = data[11]

s = '$\lambda$ = ' + str(round(lambdas,3)) + ' AU;' + ' E = ' + str(round(energy,3)) + ' MeV'
plt.text(-0.9, 0.25, s)

s = '$V_{sw}$ = ' + str(round(V_sw,3)) + ' km/s'
plt.text(-0.9, 0.0, s)

s = 'particle = ' + str(round(species,0)) 
plt.text(-0.9, -0.25, s)

s = '(1 = electrons; 2 = protons)'
plt.text(-0.9, -0.5, s)

s = 'injection = ' + str(round(injection_swtich,0))
plt.text(-0.9, -0.75, s)

s = '(1 = delta-like; 2 = Reid-Axford)'
plt.text(-0.9, -1., s)

s = '$t_a$, $t_e$ = ' + str(round(acceleration_time,1)) + ',' + str(round(escape_time,2)) + ' hrs'  
plt.text(-0.9, -1.25, s)

#s = '(N,M,Zmax,Zl,Ml) = (' + str(round(N,0)) + ',' + str(round(M,0)) + ',' + str(round(Z_max,1)) + ',' + str(round(zlimiter,0)) + ',' + str(round(mulimiter,0)) + ')'
#plt.text(-0.9, -1.25, s)

fig.savefig('Output.pdf', dpi=300, bbox_inches='tight')

#--------------------------------------------------------------------

plt.show()
