import numpy as np
from optparse import OptionParser

class grid(object):
        """ grids are defined by gridid, level, cpus, Nx, Ny, Nz"""
        def __init__(self, l):
                line=lines[l].split()
                self.index = int(line[0])
		self.level = int(line[1])
		self.cpus  = int(line[2])
		self.N = np.array(lines[l+1].split(),int)
		self.i = np.array(lines[l+2].split(),int)
		self.X = np.array(lines[l+3].split(),float)
		self.Y = np.array(lines[l+4].split(),float)
		self.Z = np.array(lines[l+5].split(),float)
		self.corners = np.array(lines[l+6].split(),int)
		self.ghcorners = np.array(lines[l+7].split(),int)
		self.cpuN =[]
		self.pcorners =[]
		self.iface =[]		
		for cpu in range(self.cpus):
			self.cpuN.append(np.array(lines[l+8+cpu*3].split(),int))
			self.pcorners.append(np.array(lines[l+9+cpu*3].split(),int))
			self.iface.append(np.array(lines[l+10+cpu*3].split(),int))

parser = OptionParser()
parser.add_option("-d", "--dir",      action="store", type="string", dest="d", default="")
parser.add_option("-i", "--dir2",     action="store", type="string", dest="i", default="")
parser.add_option("-o", "--output",     action="store", type="string", dest="o", default="1")
options, args = parser.parse_args()
directory = options.d
directory2 = options.i
sfd = options.o.zfill(5)

grids=[]
data = open (directory+"output"+sfd+"/Descriptor"+str(int(sfd))+".dat","r")
data2 = open (directory2+"output"+sfd+"/Descriptor"+str(int(sfd))+".dat","r")
global lines 
lines = data.readlines()
line=lines[1].split();ndim=int(line[0]);coordtype=int(line[1])
line=lines[2].split();Ngrid=int(line[0]);Ncpu=int(line[1]);Nghost=int(line[2]);Nvar=int(line[3]);
l=3+Nvar
for i in range(Ngrid):
	grids.append(grid(l))
	l+=8+3*grids[i].cpus
data.close()
sorted_grids = sorted(grids, key=lambda x: x.level)

maxlevel=sorted_grids[-1].level
minlevel=sorted_grids[0].level
cpu=0
failed = 0
for grid in sorted_grids:
	f  = open(directory +"output"+sfd+"/gasdensity"+str(int(sfd))+"_"+str(grid.cpuN[cpu][4])+"_"+str(grid.level)+".dat","r")
	f2 = open(directory2+"output"+sfd+"/gasdensity"+str(int(sfd))+"_"+str(grid.cpuN[cpu][4])+"_"+str(grid.level)+".dat","r")			
	n  = np.fromfile(f).reshape(grid.cpuN[cpu][2],grid.cpuN[cpu][1],grid.cpuN[cpu][0])
	m = np.fromfile(f2).reshape(grid.cpuN[cpu][2],grid.cpuN[cpu][1],grid.cpuN[cpu][0])
	diff = n-m
	ratio = n/m
	relative = np.fabs(diff)/n
	print "Desnity for grid ", grid.index, " of level ", grid.level, ": Difference = ", diff.max(), "Ratio = ",ratio.max(), "Relative difference = ", relative.max()
	if relative.max() > 10e-3:
		failed |= 1
	f  = open(directory +"output"+sfd+"/gasvelocity"+str(int(sfd))+"_"+str(grid.cpuN[cpu][4])+"_"+str(grid.level)+".dat","r")
	f2 = open(directory2+"output"+sfd+"/gasvelocity"+str(int(sfd))+"_"+str(grid.cpuN[cpu][4])+"_"+str(grid.level)+".dat","r")			
	n  = np.fromfile(f).reshape(ndim,grid.cpuN[cpu][2],grid.cpuN[cpu][1],grid.cpuN[cpu][0])
	m = np.fromfile(f2).reshape(ndim,grid.cpuN[cpu][2],grid.cpuN[cpu][1],grid.cpuN[cpu][0])
	diff = n-m
	ratio = n/m
	relative = np.fabs(diff)/n
	print "Vx for grid ", grid.index, " of level ", grid.level, ": Difference = ", diff[0].max(), "Ratio = ",ratio[0].max(), "Relative difference = ", relative[0].max()
	if relative[0].max() > 10e-13:
		failed |= 2	
	print "Vy for grid ", grid.index, " of level ", grid.level, ": Difference = ", diff[1].max(), "Ratio = ",ratio[1].max(), "Relative difference = ", relative[1].max()
	if relative[1].max() > 10e-13:
		failed |= 4	
	print "Vz for grid ", grid.index, " of level ", grid.level, ": Difference = ", diff[2].max(), "Ratio = ",ratio[2].max(), "Relative difference = ", relative[2].max()
	if relative[2].max() > 10e-13:
		failed |= 8		

					
