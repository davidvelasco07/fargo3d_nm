import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata

class field:
    """This class allows "allocate and plot" a 2D FARGO3D scalar field."""
    def __init__(self, input_file):
        """input_file = String, FARGO3D output field
           option = XY, XZ, YZ
        """
        self.file = input_file

        #Open files sections----------------

        try: 
            file = open(self.file,'r')
        except:
            print 'Error!!!',self.file, 'cannot be opened.'
        try:
            domain_x = open("domain_x.dat",'r')
        except:
            print 'Error!!!',domain_x.dat, 'cannot be opened.'
        try:
            domain_y = open("domain_y.dat",'r')
        except:
            print 'Error!!!',domain_y.dat, 'cannot be opened.'
        try:
            domain_z = open("domain_z.dat",'r')
        except:
            print 'Error!!!',domain_y.dat, 'cannot be opened.'
        try:
            dims = open("dims.dat",'r')
        except:
            print 'Error!!!',dims.dat, 'cannot be opened.'

        #end open file section--------------

        DIMS = dims.readlines() #Read all dims file
        DIMS = DIMS[1].split()

        self.data = np.fromfile(self.file, "d", )

        file.close()
        domain_y.close()
        domain_z.close()
        dims.close()

        if(int(DIMS[6]) == 1):
            self.dim1 = int(DIMS[7])
            self.dim2 = int(DIMS[8])
            self.x_axis = 'Y'
            self.y_axis = 'Z'
            self.xmin = float(DIMS[2])
            self.xmax = float(DIMS[3])
            self.ymin = float(DIMS[4])
            self.ymax = float(DIMS[5])
        if(int(DIMS[7]) == 1):
            self.dim1 = DIMS[6]
            self.dim2 = DIMS[8]
            self.x_axis = 'X'
            self.y_axis = 'Z'
            self.xmin = float(DIMS[0])
            self.xmax = float(DIMS[1])
            self.ymin = float(DIMS[4])
            self.ymax = float(DIMS[5])
        if(int(DIMS[8]) == 1):
            self.dim1 = int(DIMS[6])
            self.dim2 = int(DIMS[7])
            self.x_axis = 'X'
            self.y_axis = 'Y'
            self.xmin = float(DIMS[0])
            self.xmax = float(DIMS[1])
            self.ymin = float(DIMS[2])
            self.ymax = float(DIMS[3])

        self.data = np.reshape(self.data, [self.dim2,self.dim1])

    def plot(self, imout = None, zmin = None, zmax = None, colorbox = False, dpi = 72):
        """Ploting field function"""
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_xlabel(self.x_axis)
        ax1.set_ylabel(self.y_axis)
        im = ax1.imshow(self.data,origin="lower", \
                            extent=[self.xmin,self.xmax,self.ymin,self.ymax],\
                            vmin = zmin, vmax = zmax)
        if(colorbox == True):
            cbar = fig.colorbar(im)
        if(imout == None):
            plt.show()
        else:
            plt.savefig(imout, dpi =dpi, format = 'png')
            print imout, 'created'
        plt.close(fig) #crucial, prevents memory leaks
    
    def meridian(self, n = 500, imout = None, zmin = None, zmax = None, colorbox = False, dpi = 72):
        """cartesian projection of spherical data"""
        r = np.arange(self.xmin,self.xmax,\
                          (self.xmax-self.xmin)/self.dim1)
        t = np.arange(self.ymin,self.ymax,\
                          (self.ymax-self.ymin)/self.dim2)
        x = np.ndarray([self.dim1*self.dim2], dtype = float)
        y = np.ndarray([self.dim1*self.dim2], dtype = float)
        z = np.ndarray([self.dim1*self.dim2], dtype = float)
        k = 0
        for i in range(self.dim1):
            for j in range(self.dim2):
                x[k] = r[i]*np.sin(t[j])
                y[k] = r[i]*np.cos(t[j])
                z[k] = self.data[j,i]
                k += 1
        xx = np.arange(self.xmin,self.xmax,(self.xmax-self.xmin)/n)
        yy = np.arange(self.xmax*np.cos(self.ymax),self.xmax*np.cos(self.ymin),\
                           (self.xmax*(np.cos(self.ymin)-np.cos(self.ymax))/n))
        zz = griddata(x,y,z,xx,yy)
        fig_pol = plt.figure()
        ax1 = fig_pol.add_subplot(111)
        ax1.set_xlabel(self.x_axis)
        ax1.set_ylabel(self.y_axis)
        im = ax1.imshow(zz, origin="lower",\
                       extent=[self.xmin,self.xmax,self.xmax*np.cos(self.ymax), \
                                   self.xmax*np.cos(self.ymin)], vmin = zmin, vmax = zmax)
        if(colorbox == True):
            cbar = fig_pol.colorbar(im,shrink=0.6, aspect = 30)
        if(imout == None):
            plt.show()
        else:
            plt.savefig(imout, dpi = dpi, format = 'png')
            print imout, 'created'
        plt.close(fig_pol) #crucial, prevents memory leaks

    def polar(self, n = 500, imout = None, zmin = None, zmax = None, colorbox = False, dpi = 72):
        """cartesian projection of spherical data
        zmin = float, minval of colorbox
        zmax = float, maxval of colorbox
        colorbox = boolean, (True or False). Show or hide colorbox
        """
        newdim = (self.dim2+int(self.ymin*self.dim2/(self.ymax-self.ymin))) #Correction for black centre
        r = np.arange(0.0,self.ymax,self.ymax/newdim)
        t = np.arange(self.xmin,self.xmax,\
                          (self.xmax-self.xmin)/self.dim1)
        x = np.ndarray([newdim*self.dim1], dtype = float)
        y = np.ndarray([newdim*self.dim1], dtype = float)
        z = np.ndarray([newdim*self.dim1], dtype = float)
        k = 0
        for i in range(newdim):
            for j in range(self.dim1):
                x[k] = r[i]*np.cos(t[j])
                y[k] = r[i]*np.sin(t[j])
                if(r[i]>self.ymin):
                    z[k] = self.data[i-(newdim-self.dim2),j]
                else:
                    z[k] = 0.0
                k += 1
        xx = np.arange(-self.ymax,self.ymax,2.0*self.ymax/n)
        yy = np.arange(-self.ymax,self.ymax,2.0*self.ymax/n)
        zz = griddata(x,y,z,xx,yy)
        fig_pol = plt.figure()
        ax1 = fig_pol.add_subplot(111)
        ax1.set_xlabel(self.x_axis)
        ax1.set_ylabel(self.y_axis)
        im_pol = ax1.imshow(zz, origin="lower",\
            extent=[-self.ymax,self.ymax,-self.ymax,self.ymax], vmin = zmin, vmax = zmax)
        if(colorbox == True):
            cbar = fig_pol.colorbar(im_pol)
        if(imout == None):
            plt.show()
        else:
            plt.savefig(imout, dpi = dpi, format = 'png')
            print imout, 'created'
        plt.close(fig_pol) #crucial, prevents memory leaks


    def animate(self, init_frame = 0, end_frame = 10, projection = "plot", \
                    n = 500, zmin = None, zmax = None, colorbox = False, dpi = 72):
        for i in range(init_frame, end_frame+1):
            filename = ((self.file).split('0'))[0] + "{:06d}".format(i) + '.dat'
            try: 
                self.data = np.fromfile(filename, "d")
                self.data = np.reshape(self.data, [self.dim2,self.dim1])
            except:
                print 'Error!!!',filename, 'cannot be opened.'            
                continue
            if(projection == "plot"):
                self.plot(imout =((self.file).split('0'))[0] + "{:06d}".format(i) + '.png', \
                              zmin=zmin, zmax=zmax, colorbox=colorbox, dpi = dpi)
            if(projection == "polar"):
                self.polar(n=n,imout =((self.file).split('0'))[0] + "{:06d}".format(i) + '.png', \
                              zmin=zmin, zmax=zmax, colorbox=colorbox, dpi = dpi)
            if(projection == "meridian"):
                self.meridian(n=n,imout =((self.file).split('0'))[0] + "{:06d}".format(i) + '.png', \
                              zmin=zmin, zmax=zmax, colorbox=colorbox, dpi = dpi)


if (__name__ == "__main__"):
    f = field("Density000000.dat")
    f.animate(0,100,0.0006,0.0007,False)
