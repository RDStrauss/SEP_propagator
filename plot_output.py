import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------
data = np.loadtxt("output_at_1AU.txt")

time = data[:,1]
omni_directional_intensity = data[:,2]
anisotropy = data[:,3]
z = data[0,0]

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
s = 'Output at z = ' + str(round(z,2)) + ' AU'
subplot1.set_title(s)
subplot1.set_ylabel('Omni-directional intensity')
subplot1.set_yscale('log')
subplot1.set_xscale('linear')
#subplot1.set_ylim([0.001,10.])
#subplot1.ticklabel_format(useOffset=False)
#subplot1.set_xlim([-0.5,10])

subplot1.plot(time,omni_directional_intensity, color = 'red')

subplot1 = fig.add_subplot(223)
  
subplot1.set_xlabel('Time (hr)')
subplot1.set_ylabel('anisotropy')
subplot1.set_yscale('linear')
subplot1.set_xscale('linear')
#subplot1.set_ylim([0.001,10.])
#subplot1.ticklabel_format(useOffset=False)
#subplot1.set_xlim([-0.5,10])
  
subplot1.plot(time, anisotropy, color = 'blue')
subplot1.axhline(y = 0, color = 'black', linestyle = '--')

subplot1 = fig.add_subplot(122)

s = 'Output at z = ' + str(round(z,2)) + ' AU'
subplot1.set_title(s)
subplot1.set_xlabel('$\mu$')
#subplot1.set_ylabel('anisotropy')
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
plt.legend()

fig.savefig('Output.pdf', dpi=300)

#--------------------------------------------------------------------

plt.show()
