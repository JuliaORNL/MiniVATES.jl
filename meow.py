import numpy as np

import matplotlib.cm as cm
import matplotlib.pyplot as plt

Z = np.loadtxt('meow.txt')
Z = np.reshape(Z,(200,200))

fig, ax = plt.subplots()
im = ax.imshow(Z, interpolation='none', cmap=cm.viridis,
               origin='lower', extent=[-10, 10, -10, 10],
#               vmax=np.nanmax(Z), vmin=np.nanmin(Z))
               vmax=1e-5, vmin=0)

print(np.nanmin(Z),np.nanmax(Z))
plt.xlabel('[H,H,0] (r.l.u.)')
plt.ylabel('[H,-H,0] (r.l.u.)')
plt.title('MDNorm Proxy Output')
#plt.legend(loc='upper right')
#plt.xlim(0,5000)
#plt.ylim(0,0.2)
fig = plt.gcf()
fig.set_size_inches(6, 6)
plt.savefig('mdnorm_proxy.png')
plt.show()

