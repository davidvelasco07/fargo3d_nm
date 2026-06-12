### Run only this parameter file is FARGO3D was build using the setup 'fargo'
Setup			fargo

### Disk parameters

AspectRatio     	0.05            Thickness over Radius in the disc
Sigma0			6.3661977237e-4	Surface Density at r=1
Nu			0.0		Uniform kinematic viscosity
SigmaSlope		0.0		Slope of the surface density
FlaringIndex		0.0

### Planet parameters

PlanetConfig		planets/jupiter.cfg
ThicknessSmoothing 	0.6
RocheSmoothing 		0.0
Eccentricity		0.0
ExcludeHill		no
IndirectTerm		Yes

#Noise			0.8
### Mesh parameters

Nx			64		Azimuthal number of zones
Ny               	32		Radial number of zones
Xmin			-3.14159265358979323844	
Xmax			3.14159265358979323844
Ymin			0.4		Inner boundary radius
Ymax			2.5		Outer boundary radius
OmegaFrame     		1.0005
Frame			G

### Output control parameters

Ntot			100		Total number of time steps
Ninterm	 		10		Time steps between outputs
DT			0.05
OutputDir		@outputs/fargo

#Plotting parameters

Log			yes

#GRIDINFO
-1 0.7 1 1.4 1
-0.5 0.85 0.5 1.2 2
