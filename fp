### Run only this parameter file is FARGO3D was build using the setup 'fargo'
Setup			p3diso

### Disk parameters

AspectRatio     	0.05            Thickness over Radius in the disc
Sigma0			6.3661977237e-4	Surface Density at r=1
Nu			1.0e-6		Uniform kinematic viscosity
SigmaSlope		1.5		Slope of surface
FlaringIndex		0.0
Noise			0.0
### Planet parameters

PlanetConfig		planets/SuperEarth.cfg
#PlanetConfig		planets/Jupiter.cfg
ThicknessSmoothing 	0.1

### Numerical method parameters

Disk			YES
OmegaFrame     		1.0
Frame			C
IndirectTerm		No

### Mesh parameters

Nx			120		Azimuthal number of zones
Ny               	80		Radial number of zones
Nz			15
Xmin			-3.14159265359	
Xmax			3.14159265359
Ymin			0.4		Inner boundary radius
Ymax			2.5		Outer boundary radius
Zmin			1.42079632679
Zmax			1.57079632679489661922

### Output control parameters

Ntot			200		Total number of time steps
Ninterm	 		1		Time steps between outputs
DT			0.31416
OutputDir		testfargo

#GRIDINFO
-0.61229056 0.8 1.42079632679 0.61229056 1.2 1.57079632679 1 
-0.3874048 0.8425952 1.42079632679 0.3874048 1.1274048 1.57079632679 2 
#-0.185184 0.914816 1.48561232679 0.185184 1.085184 1.57079632679 3 
#-0.03872 0.96128 1.53207632679 0.03872 1.03872 1.57079632679 4 
#-0.0176 0.9824 1.55319632679 0.0176 1.0176 1.57079632679 5 
#-0.008 0.992 1.56279632679 0.008 1.008 1.57079632679 6 
#-2.5 0.7 1.47 1.5  1.3 1.57079632679489661922 1
#-0.5 1. 1.47 2.5  1.6 1.57079632679489661922 1
#-0.5 0.85 1.50 0.5 1.15 1.57079632679489661922 2
#-0.25 0.905 1.53 0.1 1.075 1.57079632679489661922 3
#-0.1 0.905 1.53 0.25 1.075 1.57079632679489661922 3
