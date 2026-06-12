#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Rectangle, Wedge
from matplotlib.widgets import Slider, RadioButtons, CheckButtons
import scipy.ndimage as ndimage
import math
from optparse import OptionParser
import sys

if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

root = Tk.Tk()
root.wm_title("Embedding in TK")


fig = plt.figure(figsize=(15,10))
ax1 = plt.subplot()

# a tk.DrawingArea
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.show()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
toolbar = NavigationToolbar2TkAgg(canvas, root)
toolbar.update()
canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

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

def slicing(Plane,m3d,index):
	if Plane == "XY":
		return m3d[index,:,:]
	elif Plane == "YZ":
		return m3d[:,:,index]
	elif Plane == "ZX":
		return m3d[:,index,:].transpose()

def initiate(grids,minlevel,Plane):
	if Plane == "XY":
		return grids[minlevel].Z[Nghost]
	elif Plane == "YZ":
		return grids[minlevel].X[Nghost]
	elif Plane == "ZX":
		return grids[minlevel].Y[Nghost]

def get_index(Plane):
	if Plane == "XY":
		return 2
	elif Plane == "YZ":
		return 0
	elif Plane == "ZX":
		return 1

#Command line parser/////////////////////////////////////////
parser = OptionParser()
parser.add_option("-o", "--output", action="store", type="string", dest="o", default ="0")
parser.add_option("-p", "--plane",  action="store", type="string", dest="p", default ="XY")
parser.add_option("-d", "--dir",    action="store", type="string", dest="d", default="outputs/fargo/")
parser.add_option("-i", "--dir2",   action="store", type="string", dest="i", default="")
parser.add_option("-c", "--cut",    action="store", type="string", dest="c", default="")
parser.add_option("-r", "--rest",    action="store", type="string", dest="r", default="")
options, args = parser.parse_args()
sfd = options.o.zfill(5)
Plane = options.p
directory = options.d
directory2 = options.i
cut = options.c
rest = options.r
if rest != "":
	rest=rest.zfill(5)
Q = "density"
vel = 0
vector=False
draw_levels = 0
print_text = 0
draw_grids = 0
function = ""
polar = 0
aspect = 0
colormap = cm.gnuplot
colormap.set_bad('w',1.)
#////////////////////////////////////////////////////////////


#Reading Descriptor and filling grid obejects////////////////////
grids=[]
data = open (directory+"output"+sfd+"/Descriptor"+str(int(sfd))+".dat","r")
if directory2 != "":
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

if coordtype == 0:
	label=['x','y','z']
elif coordtype == 1:
	label=['azimuth','radius','vertical']
elif coordtype == 2:
	label=['azimuth','radius','colatitude']

sorted_grids = sorted(grids, key=lambda x: x.level)

maxlevel=sorted_grids[-1].level
minlevel=sorted_grids[0].level
#/////////////////////////////////////////////////////////////
def matrix_builder(level, cpu, output):
	q=[]
	matrix=[]
	if Q != 'temp':
		q.append(Q)
	else:
		q.append("energy")
		q.append("density")
	for quantity in q:
		f = open(directory+"output"+output+"/gas"+quantity+str(int(output))+"_"+str(level.cpuN[cpu][4])+"_"+str(level.level)+".dat","r")
		if directory2 != "":	
			f2 = open(directory2+"output"+output+"/gas"+quantity+str(int(output))+"_"+str(level.cpuN[cpu][4])+"_"+str(level.level)+".dat","r")	
		if vector == False:		
			m3d = np.fromfile(f).reshape(level.cpuN[cpu][2],level.cpuN[cpu][1],level.cpuN[cpu][0])
			if directory2 != "":
				n3d = np.fromfile(f2).reshape(level.cpuN[cpu][2],level.cpuN[cpu][1],level.cpuN[cpu][0])
				m3d = (m3d-n3d)#/m3d
		else:
			m4d = ((np.fromfile(f).reshape(ndim,level.cpuN[cpu][2],level.cpuN[cpu][1],level.cpuN[cpu][0])))
			if directory2 != "":
				n4d = np.fromfile(f2).reshape(ndim,level.cpuN[cpu][2],level.cpuN[cpu][1],level.cpuN[cpu][0])
				m4d = np.fabs(m4d-n4d)#/m4d
			m3d = m4d[vel,:,:,:]
		matrix.append(m3d)
	if Q == 'temp':
		m3d = matrix[0]/matrix[1]				
	return m3d


#Principal plotting function//////////////////////////////////
def plotting(cut):
	m=[]
	maximums=[]
	minimums=[]
	maxvalue=0
        minvalue=0
	minimum=0
	maximum=0
	for level in sorted_grids:
		if level.level <= maxlevel and level.level >= minlevel:
			lattice=[]
			lattice.append((level.X[-(Nghost+1)]-level.X[Nghost])/level.N[0])
			lattice.append((level.Y[-(Nghost+1)]-level.Y[Nghost])/level.N[1])
			lattice.append((level.Z[-(Nghost+1)]-level.Z[Nghost])/level.N[2])
			if Plane == "XY":
				rcorners = [level.Z,level.Y,level.X]
				dim = [2,1,0]
			elif Plane == "YZ":
				rcorners = [level.X,level.Z,level.Y]
				dim = [0,2,1]
			elif Plane == "ZX":
				rcorners = [level.Y,level.X,level.Z]
				dim = [1,0,2]
			else:
				print "Incorrect set of coordinates"
				break
			if level.level == minlevel:
				if polar == 0 or Plane != "XY" :
					ax1.axis([rcorners[2][Nghost],rcorners[2][-(Nghost+1)],rcorners[1][Nghost],rcorners[1][-(Nghost+1)]])
					ratio = (rcorners[2][Nghost]-rcorners[2][-(Nghost+1)])/(rcorners[1][Nghost]-rcorners[1][-(Nghost+1)])
				else:
					ax1.axis([-rcorners[1][-(Nghost+1)],rcorners[1][-(Nghost+1)],-rcorners[1][-(Nghost+1)],rcorners[1][-(Nghost+1)]])
					ratio=1
                                
				maxvalue = rcorners[0][-(Nghost+1)]
				minvalue = rcorners[0][Nghost]
					
			if cut>=rcorners[0][Nghost] and cut<=rcorners[0][-(Nghost+1)]:
				for cpu in range(level.cpus):
					if ( cut>=rcorners[0][Nghost+level.pcorners[cpu][dim[0]]] and cut<=rcorners[0][Nghost+level.pcorners[cpu][3+dim[0]]]):
						if level.N[dim[0]] > 1:
							diff, index = np.modf((cut-rcorners[0][Nghost+level.pcorners[cpu][dim[0]]]-0.5*lattice[dim[0]])/float(lattice[dim[0]]))
						else:
							diff = 0
							index = 0
                                                #diff = 0        
						corner1 = Nghost+level.pcorners[cpu][dim[2]]
						corner2 = Nghost+level.pcorners[cpu][dim[2]+3]
						corner3 = Nghost+level.pcorners[cpu][dim[1]]
						corner4 = Nghost+level.pcorners[cpu][dim[1]+3]
						m3d = matrix_builder(level, cpu, sfd)
						if rest != "":
							m3d -= matrix_builder(level, cpu, rest)
						if diff < 0:
							m2d = diff*slicing(Plane,m3d,index+1)+(1-diff)*slicing(Plane,m3d,index)
						elif index == level.cpuN[cpu][dim[0]]-1:	
							m2d = (1-diff)*slicing(Plane,m3d,index-1)+diff*slicing(Plane,m3d,index)
						else:
							m2d = (1-diff)*slicing(Plane,m3d,index)+diff*slicing(Plane,m3d,index+1)
						if function == "log":
							m2d = np.log10(m2d+0.0001)
						m2d = np.ma.array (m2d, mask=np.isnan(m2d))
						y, x = np.meshgrid(rcorners[1][corner3:corner4+1],rcorners[2][corner1:corner2+1])
						if polar == 1 and Plane == "XY":
							y, x = y * np.sin(x), y * np.cos(x)				
						if draw_grids == 0:
                                                        m.append(ax1.pcolormesh(x,y,m2d.T,cmap=colormap,zorder = level.index))
                                                        a=1
						else:
							m.append(ax1.pcolormesh(x,y,m2d.T,cmap=colormap,edgecolor="w",lw=0.001,zorder = level.index,alpha=1.0))	
                                                
                                                if print_text == 1:
                                                        if  polar == 0 or Plane != "XY":
                                                                ax1.text(rcorners[2][corner1],rcorners[1][corner3],"CPU " + str(level.cpuN[cpu][3]),color="b",fontsize=10-(level.level-minlevel),zorder=level.level,weight="heavy")
                                                                ax1.add_patch(Rectangle((rcorners[2][corner1], rcorners[1][corner3]), rcorners[2][corner2]-rcorners[2][corner1], rcorners[1][corner4]-rcorners[1][corner3],ec="b",lw=1.5,fc="none",zorder=level.index))
                                                        else:
                                                                 ax1.add_patch(Wedge(0,rcorners[1][corner4],rcorners[2][corner1]*(180/np.pi),rcorners[2][corner2]*(180/np.pi),rcorners[1][corner4]-rcorners[1][corner3],fill=False,ec="b",zorder=level.level))
                                                maximums.append(m2d.max())
						minimums.append(m2d.min())
                                if draw_levels == 1:
                                        if polar == 0 or Plane != "XY":
                                                ax1.add_patch(Rectangle([rcorners[2][Nghost],rcorners[1][Nghost]],rcorners[2][-Nghost-1]-rcorners[2][Nghost],rcorners[1][-Nghost-1]-rcorners[1][Nghost],zorder=level.index,ec="r",fill=False))
                                        else:
                                                ax1.add_patch(Wedge(0,rcorners[1][-Nghost-1],rcorners[2][Nghost]*(180/np.pi),rcorners[2][-Nghost-1]*(180/np.pi),rcorners[1][-Nghost-1]-rcorners[1][Nghost],fill=False,ec="r",zorder=level.level))
					
	if len(m) != 0:	
		maximum = max(maximums)
		minimum = min(minimums)
		if aspect == 1:
			ax1.set_aspect(float(ratio))
		else:
			ax1.set_aspect('auto')
		ax1.set_title("Output:"+str(int(sfd))+"\n minlevel = "+str(minlevel)+" maxlevel = "+str(maxlevel)+" minvalue = "+str(minimum)+" maxvalue = "+str(maximum))
		if polar == 1 and Plane == "XY":
			ax1.set_xlabel("X")
			ax1.set_ylabel("Y")
		else:
			ax1.set_xlabel(label[dim[2]])
			ax1.set_ylabel(label[dim[1]])
		for image in m:
			image.set_clim(vmax=maximum,vmin=minimum)
	return maxvalue, minvalue, m, [minimum,maximum]
#////////////////////////////////////////////////////////////////////////////////////////

def pick_value(x,y,cut):
        global Plane
	for level in sorted_grids:
		if level.level <= maxlevel and level.level >= minlevel:
			lattice=[]
			lattice.append((level.X[-(Nghost+1)]-level.X[Nghost])/level.N[0])
			lattice.append((level.Y[-(Nghost+1)]-level.Y[Nghost])/level.N[1])
			lattice.append((level.Z[-(Nghost+1)]-level.Z[Nghost])/level.N[2])
			if Plane == "XY":
				rcorners = [level.Z,level.Y,level.X]
				dim = [2,1,0]
			elif Plane == "YZ":
				rcorners = [level.X,level.Z,level.Y]
				dim = [0,2,1]
			elif Plane == "ZX":
				rcorners = [level.Y,level.X,level.Z]
				dim = [1,0,2]
			else:
				print "Incorrect set of coordinates"
				break
				
			if cut>=rcorners[0][Nghost] and cut<=rcorners[0][-(Nghost+1)]:
				for cpu in range(level.cpus):
					if (  cut>=rcorners[0][Nghost+level.pcorners[cpu][dim[0]]] and cut<=rcorners[0][Nghost+level.pcorners[cpu][3+dim[0]]] and\
						y>=rcorners[1][Nghost+level.pcorners[cpu][dim[1]]] and   y<=rcorners[1][Nghost+level.pcorners[cpu][3+dim[1]]] and\
					       	x>=rcorners[2][Nghost+level.pcorners[cpu][dim[2]]] and   x<=rcorners[2][Nghost+level.pcorners[cpu][3+dim[2]]]):
						if level.N[dim[0]] > 1:
							diff, index = np.modf((cut-rcorners[0][Nghost+level.pcorners[cpu][dim[0]]]-0.5*lattice[dim[0]])/float(lattice[dim[0]]))
							index=int(index)
						else:
							diff = 0
							index = 0
						corner1 = Nghost+level.pcorners[cpu][dim[2]]
						corner2 = Nghost+level.pcorners[cpu][dim[2]+3]
						corner3 = Nghost+level.pcorners[cpu][dim[1]]
						corner4 = Nghost+level.pcorners[cpu][dim[1]+3]
						m3d = matrix_builder(level, cpu, sfd)
						if rest != "":
							m3d -= matrix_builder(level, cpu, rest)
						i = int((x-rcorners[2][corner1])/lattice[dim[2]])
						j = int((y-rcorners[1][corner3])/lattice[dim[1]])
						if Plane == "XY":
							print "Level ",level.level," CPU ", level.cpuN[cpu][3], "CPU GRID", level.cpuN[cpu][4]," position ", [i,j,index], m3d[index,j,i]
						elif Plane == "YZ":
							print "Level ",level.level," CPU ", level.cpuN[cpu][3], "CPU GRID", level.cpuN[cpu][4]," position ", [index,i,j], m3d[j,i,index]
						elif Plane == "ZX":
							print "Level ",level.level," CPU ", level.cpuN[cpu][3], "CPU GRID", level.cpuN[cpu][4]," position ", [j,index,i], m3d[i,index,j]
							
if(cut == ""):
	initial_cut = initiate(sorted_grids,0,Plane)
else:
	initial_cut = float(cut)
maxvalue, minvalue, images, clim = plotting(initial_cut)
axcbar = plt.axes([ax1.get_position().x1, 0.1, 0.03, 0.8])
cbar=fig.colorbar(images[0],cax=axcbar)

axcolor = 'khaki'
axDIM = plt.axes([0.12, 0.025, 0.25, 0.025], axisbg='slateblue')
axRB = plt.axes([ax1.get_position().x0-0.12, 0.025, 0.09, 0.35], axisbg=axcolor)
axCB = plt.axes([ax1.get_position().x0-0.12, 0.55, 0.09, 0.12], axisbg="tan")

if ndim == 3:
	sDIM = Slider(axDIM, label[get_index(Plane)], minvalue, maxvalue, valinit=initial_cut)		
	sRB = RadioButtons(axRB, ('Density','log(Density)','T','Energy','Vx','Vy','Vz','ERad','Potential','Mpx','Mmx','Mpy','Mmy','Mpz','Mmz'),activecolor='darkorchid')
        
else:
	sDIM = Slider(axDIM, label[get_index(Plane)], initial_cut,initial_cut, valinit=initial_cut)
	sRB = RadioButtons(axRB, ('Density','log(Density)','Energy','Vx','Vy'),activecolor='blue')

sCB = CheckButtons(axCB, ('Levels', 'CPUs', 'Grid', 'Polar', 'Aspect'), (draw_levels, print_text, draw_grids, polar, aspect))
for circle in sRB.circles:
        circle.set_radius(0.025)
def update(val):
	global maxvalue, minvalue, images, clim
	ax1.cla()
	maxvalue, minvalue, images, clim = plotting(sDIM.val)
	cbar.set_clim(clim)
	cbar.update_normal(images)
	fig.canvas.draw_idle()


def press(event):
	sys.stdout.flush()
	global maxlevel,minlevel,sfd
	if event.key=='right' and maxlevel<Ngrid-1:
		maxlevel+=1
	elif event.key=='left' and maxlevel>0 and maxlevel>minlevel:
		maxlevel-=1	
        elif event.key=='up' and minlevel<Ngrid-1 and minlevel<maxlevel:
	    	minlevel+=1
	elif event.key=='down' and minlevel>0:
	    	minlevel-=1
        elif event.key=="d":
                n=int(sfd)+1
                s=str(n).zfill(5)
                try: 
                        data = open (directory+"output"+s+"/Descriptor"+str(int(s))+".dat","r")
                except:
                        exit
		sfd = s
        elif event.key=="a":
                n=int(sfd)-1
                if n>=0:
                        sfd=str(n).zfill(5)
                else:
                        exit
        else:
                exit
	update(sDIM.val)


def click(event):
    	global Plane, ax1, axDIM, sDIM, maxvalue, minvalue, images, clim, canvas
	canvas.get_tk_widget().focus_force()
	if event.inaxes == axDIM and ndim == 3:
		sDIM.on_changed(update)

	

def pick(event):
	global Plane, ax1, axDIM, sDIM, maxvalue, minvalue, images, clim, polar, cbar, var, limits
	val = event.mouseevent.ydata
	if event.mouseevent.button == 1 and event.mouseevent.inaxes == ax1 and ndim == 3:
		limits.focus()
		if  polar == 1 and coordtype != 0 and Plane == "XY":
			x = np.arctan2(event.mouseevent.ydata,event.mouseevent.xdata)
		else:
			x = event.mouseevent.xdata
		if Plane == "XY":
			Plane = "YZ"
		elif Plane == "YZ":
			Plane = "ZX"
		elif Plane == "ZX":
			Plane = "XY"
		ax1.cla()
		axDIM.cla()
		maxvalue,minvalue,images,clim = plotting(x)
		var.set("Range: ("+str(minvalue)+","+str(maxvalue)+")")
		limits.pack()
		sDIM = Slider(axDIM, label[get_index(Plane)], minvalue, maxvalue, valinit=x)
		cbar.set_clim(clim)
		cbar.update_normal(images)
		fig.canvas.draw_idle()

	elif event.mouseevent.button == 2 and event.mouseevent.inaxes == ax1:
		if polar == 1 and coordtype != 0 and Plane == "XY":
			x = np.arctan2(event.mouseevent.ydata,event.mouseevent.xdata)
			y = np.sqrt(event.mouseevent.ydata**2+event.mouseevent.xdata**2)
			print x, y
		else:
			x = event.mouseevent.xdata
			y = event.mouseevent.ydata
		pick_value(x,y,sDIM.val)
		
	elif event.mouseevent.button == 3 and event.mouseevent.inaxes == ax1 and ndim == 3:
		if  polar == 1 and coordtype != 0 and Plane == "XY":
			y = np.sqrt(event.mouseevent.ydata**2+event.mouseevent.xdata**2)
		else:
			y = event.mouseevent.ydata
		if Plane == "XY":
			Plane = "ZX"
		elif Plane == "YZ":
			Plane = "XY"
		elif Plane == "ZX":
			Plane = "YZ"
		ax1.cla()
		axDIM.cla()
		maxvalue,minvalue,images,clim = plotting(y)
		sDIM = Slider(axDIM, label[get_index(Plane)], minvalue, maxvalue, valinit=y)
		cbar.set_clim(clim)
		cbar.update_normal(images)
		fig.canvas.draw()
		
	if event.mouseevent.button == 1 and event.mouseevent.inaxes == cbar.ax:
		current_clim = cbar.get_clim()
		val = (current_clim[1]-current_clim[0])*val+current_clim[0]
		cbar.set_clim(vmax=val)
   		cbar.update_normal(images)
		for image in images:
			image.set_clim(vmax=val)
    		fig.canvas.draw()
	if event.mouseevent.button == 3 and event.mouseevent.inaxes == cbar.ax:
		current_clim = cbar.get_clim()
		val = (current_clim[1]-current_clim[0])*val+current_clim[0]
		cbar.set_clim(vmin=val)
   		cbar.update_normal(images)
		for image in images:
			image.set_clim(vmin=val)
    		fig.canvas.draw()
	if event.mouseevent.button == 2 and event.mouseevent.inaxes == cbar.ax:
		cbar.set_clim(clim)
   		cbar.update_normal(images)
		for image in images:
			image.set_clim(clim)
    		fig.canvas.draw()
	

def change_field(label):
	global Q, vel, ax1, maxvalue, minvalue, images, clim, cut, function, vector
	if label == "Density":
		Q="density"
                vector=False
		function=""
	elif label == "log(Density)":
		Q="density"
                vector=False
		function="log"
	elif label == "T":
		Q="temp"
                vector=False
		function=""
	elif label == "Energy":
		Q="energy"
                vector=False
		function=""
	elif label == "Vx":
		Q="velocity"
		vel = 0
                vector=True
		function=""
	elif label == "Vy":
		Q="velocity"
		vel = 1
                vector=True
		function=""
	elif label == "Vz":
		Q="velocity"
		vel = 2
                vector=True
		function=""
	elif label == "ERad":
                Q="energyrad"
                vector=False
		function=""
	elif label == "Potential":
		Q="potential"
                vector=False
		function=""
        elif label == "Mpx":
		Q="momentap"
		vel = 0
                vector=True
		function=""
	elif label == "Mpy":
		Q="momentap"
		vel = 1
                vector=True
		function=""
	elif label == "Mpz":
		Q="momentap"
		vel = 2
                vector=True
		function=""  
        elif label == "Mmx":
		Q="momentam"
		vel = 0
                vector=True
		function=""
	elif label == "Mmy":
		Q="momentam"
		vel = 1
                vector=True
		function=""
	elif label == "Mmz":
		Q="momentam"
		vel = 2
                vector=True
		function=""  
	
        ax1.cla()
	maxvalue, minvalue, images, clim = plotting(sDIM.val)
	cbar.set_clim(clim)
	cbar.update_normal(images)
	fig.canvas.draw_idle()

def change_drawing(label):
	global draw_levels, print_text, draw_grids, polar, aspect, ax1, maxvalue, minvalue, images, clim, cut, function
	if label == "Levels":
		draw_levels = not draw_levels
	if label == "CPUs":
		print_text = not print_text
	if label == "Grid":
		draw_grids = not draw_grids
	if label == "Polar":
		polar = not polar
	if label == "Aspect":
		aspect = not aspect
	ax1.cla()
	maxvalue, minvalue, images, clim = plotting(sDIM.val)
	cbar.set_clim(clim)
	cbar.update_normal(images)
	fig.canvas.draw_idle()

def read_entry():
	global ax1, maxvalue, minvalue, images, clim, entry, sDIM
	try:
		s = float(entry.get())
	except:
		canvas.get_tk_widget().focus_force()
		return
	if s <= maxvalue and s >= minvalue:
		axDIM.cla()	
		sDIM = Slider(axDIM, label[get_index(Plane)], minvalue, maxvalue, valinit=s)
		ax1.cla()
		maxvalue, minvalue, images, clim = plotting(sDIM.val)
		cbar.set_clim(clim)
		cbar.update_normal(images)
		fig.canvas.draw_idle()
		canvas.get_tk_widget().focus_force()
		

cbar.ax.set_picker(5)
ax1.set_picker(5)
fig.canvas.mpl_connect('pick_event', pick)
fig.canvas.mpl_connect('key_press_event', press)
fig.canvas.mpl_connect('button_press_event', click)
sDIM.on_changed(update)
sRB.on_clicked(change_field)
sCB.on_clicked(change_drawing)
entry = Tk.Entry(root, bd =5)
entry.pack(side=Tk.LEFT)
button = Tk.Button(root, text='Slice', command=read_entry)
button.pack(side=Tk.LEFT)
var = Tk.StringVar()
limits = Tk.Label(root, textvariable=var)
var.set("Range: ("+str(minvalue)+","+str(maxvalue)+")")
limits.pack(side=Tk.LEFT)
Tk.mainloop()
