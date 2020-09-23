from pylab import *

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
time      = [0,]

for n in range(20000):
    rhod1.append(  fromfile( path+"dust1vy{:}.dat".format(n)).reshape(nz,ny,nx) [0,0,0] )
    rhod2.append(  fromfile( path+"dust2vy{:}.dat".format(n)).reshape(nz, ny,nx) [0,0,0] )
    rhog.append(   fromfile( path+"gasvy{:d}.dat".format(n)).reshape(nz, ny,nx)  [0,0,0] )

    rhod1_org.append(  fromfile( path_org+"dust1vy{:}.dat".format(n)).reshape(nz,ny,nx) [0,0,0] )
    rhod2_org.append(  fromfile( path_org+"dust2vy{:}.dat".format(n)).reshape(nz, ny,nx) [0,0,0] )
    rhog_org.append(   fromfile( path_org+"gasvy{:d}.dat".format(n)).reshape(nz, ny,nx)  [0,0,0] )

    t+=dt
    time.append(t)
    
rhod1 = array(rhod1)
rhod2 = array(rhod2)
rhog = array(rhog)

rhod1_org = array(rhod1_org)
rhod2_org = array(rhod2_org)
rhog_org= array(rhog_org)

## Plot  
fig = figure(figsize=(5.5,5.5))

ax1 = fig.add_axes([0.05,0.9, 0.9,0.9])

ax1.plot(time, rhod1, '-')
ax1.plot(time, rhod2, '-')
ax1.plot(time,  rhog , '-')

ax1.plot(time, rhod1_org, '--')
ax1.plot(time, rhod2_org, '--')
ax1.plot(time, rhog_org , '--')

savefig("damping_FARGO3D.png",  bbox_inches='tight', dpi=250)
