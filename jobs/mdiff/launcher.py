import re
import grid
import numpy as np
import os

Parameters = {}        
Parameters['A'] = [-1.0]
Parameters['B'] = [-0.3183098861837907]
Parameters['ETA'] = np.linspace(5e-6,5e-3,100)

g = grid.Grid(Parameters)

print "Grid created..."

pbs = open("template.pbs",'r')
PBS = pbs.readlines()
pbs.close()

for par in g.parnames:
    ifile = par[:-4]
    pbs = open(ifile,'w')
    for line in PBS:
       # print line.lower(),
        if re.search("%NAME",line):
            line = line.replace("%NAME",par.split("/")[-1][:-4])
            pbs.write(line)
            continue
        if re.search("%PARFILE",line):
            line = line.replace("%PARFILE","jobs/mdiff/"+par)
            pbs.write(line)
            continue
        pbs.write(line)
    pbs.close()
    print "qsub " + ifile
    os.system("qsub "+ifile)
