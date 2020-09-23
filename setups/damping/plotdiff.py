from pylab import *
import numpy as np

### Data
path    = "../../outputs/damping/"
path_org    = "../../../tmp/public/outputs/damping/"



nx      = 16
ny      = 16
nz      = 16
dt      = 0.1

rhod1 = []
rhod2 = []
rhog  = []

rhod1_org = []
rhod2_org = []
rhog_org  = []

time     = []
t = 0
for n in range(20000):
    rhod1.append(  fromfile( path+"dust1vy{:}.dat".format(n)).reshape(nz,ny,nx) [0,0,0] )
    rhod2.append(  fromfile( path+"dust2vy{:}.dat".format(n)).reshape(nz, ny,nx) [0,0,0] )
    rhog.append(   fromfile( path+"gasvy{:d}.dat".format(n)).reshape(nz, ny,nx)  [0,0,0] )

    rhod1_org.append(  fromfile( path_org+"dust1vy{:}.dat".format(n)).reshape(nz,ny,nx) [0,0,0] )
    rhod2_org.append(  fromfile( path_org+"dust2vy{:}.dat".format(n)).reshape(nz, ny,nx) [0,0,0] )
    rhog_org.append(   fromfile( path_org+"gasvy{:d}.dat".format(n)).reshape(nz, ny,nx)  [0,0,0] )

    time.append(t)

    t += dt
    
rhod1 = array(rhod1)
rhod2 = array(rhod2)
rhog = array(rhog)

rhod1_org = array(rhod1_org)
rhod2_org = array(rhod2_org)
rhog_org= array(rhog_org)

## Plot  
fig = figure(figsize=(5.5,5.5))

ax1 = fig.add_axes([0.05,0.9, 0.9,0.9])

ax1.plot(time,  (abs(rhod1 - rhod1_org))/abs(rhod1) , '-')
ax1.plot(time,  (abs(rhod2 - rhod2_org))/abs(rhod2) , '-',)
ax1.plot(time,  (abs(rhog  - rhog_org ))/abs(rhog) , '-',lw=2, label='Error obtained with FARGO3D')

### python data from gauss.py
np_load_old = np.load
np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, **k)
pydata_vydiff = np.load("python_data_vy.npy")
pydata_time   = np.load("python_data_time.npy")

ax1.plot(pydata_time,  pydata_vydiff , '--', label='Error obtained from gauss.py')

ax1.legend(frameon=False)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel("Output")
ax1.set_ylabel("relative difference")

savefig("FARGO3D_diff.png",  bbox_inches='tight', dpi=250)
