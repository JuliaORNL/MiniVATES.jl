import numpy as np
import sys
import matplotlib.cm as cm
import matplotlib.pyplot as plt

N = 1
if len(sys.argv) > 1:
    N = int(sys.argv[1])

for i in range(N):
    with open(f"purr_{i}.txt") as f:
        shape = [int(line) for line in f.readlines()]
    X = np.loadtxt(f"meow_{i}.txt")
    X = np.reshape(X,shape)
    if i == 0:
        Z = X
    else:
        Z = np.append(Z, X, axis=0)

Z = np.transpose(Z)
fig, ax = plt.subplots()
im = ax.imshow(Z, interpolation='none', cmap=cm.viridis,
           origin='lower',
           # extent=[-10, 10, -10, 10],
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

