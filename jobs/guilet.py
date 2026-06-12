#!/usr/bin/python

import sys
import re

p = re.compile ("\[(\d+)\]")
q= p.findall(sys.argv[1])
index= int(q[0])
index_nu  = index/10
index_eta = index%10
nu = 6.0e-4*10.0**((float(index_nu)-5.0)/5.0)
eta= 6.0e-4*10.0**((float(index_eta)-5.0)/5.0)
prandtl = nu/eta
f=open("../in/guilet"+str(index)+".par","w")
f.write(
"""
##############
Setup guilet
##############

### Disk parameters

AspectRatio     	0.1            Thickness over Radius in the disc
Sigma0			6e-4		Surface Density at r=1
Nu			"""+str(nu)+"""		Uniform kinematic viscosity
SigmaSlope		0.5		Slope of the surface density
FlaringIndex		0.0
VerticalDamping		0.00
MagneticSlope		-1.0
Beta			50.0
Prandtl			"""+str(prandtl)+"""
Noise			0.0

### Planet parameters

PlanetConfig		planets/guilet2.cfg
ThicknessSmoothing 	0.2
RocheSmoothing 		0.0
Eccentricity		0.0
Inclination		0.0
ReleaseDate		0.0
ReleaseRadius		0.0
ForcedCircular		no
ExcludeHill		no
IndirectTerm		no

MassTaper		6.28

### Numerical method parameters

Gamma			1.66666667

### Mesh parameters

Nx			320		Azimuthal number of zones
Ny               	256		Radial number of zones
Xmin			-1.57079632679489661922	
Xmax			1.57079632679489661922
Ymin			0.5		Inner boundary radius
Ymax			2.0		Outer boundary radius
PeriodicY		no		Do not change it
PeriodicZ		yes		Do not change it
OmegaFrame     		1.00014998875168718366
Frame			F

### Output control parameters

Ntot			600		Total number of time steps
Ninterm	 		5		Time steps between outputs
DT			0.314159265359	Time step length. 2PI = 1 orbit
OutputDir		/state/partition1/frederic/guilet1_"""+str(index)+"""

#Plotting parameters

PlotLine		field[0,:,:]
field			gasdens
cmap			cubehelix
funcarchfile		std/func_arch.cfg
""")
f.close()

