import re
import numpy as np
import itertools
import os

class Grid():
    
    def __init__(self,parameters,template_par='template.par', template_pbs='template.pbs'):

        """
        Class for generating a set of cube of .par & .pbs
        files. parameters is a dictionary with all the parameters you
        want to modify.  template is a .par standard file taken as a
        model. Only will be modified the OutputDir param, and all the
        parameters that matches with the keys of the dictionary
        parameters.
        ===================================================
        Important Notes: 1) OutputDir is automatically changed.
                         2) The script assumes you are inside projects/setupname
        ===================================================
        """

        self.parameters = parameters
        self.template_par = template_par
        self.template_pbs = template_pbs

    def execute(self):
       
        parameters   = self.parameters
        template_par = self.template_par
        template_pbs = self.template_pbs

        template = open(template_par)
        TEMPLATE_PAR = template.readlines()
        template.close()

        template = open(template_pbs)
        TEMPLATE_PBS = template.readlines()
        template.close()

        setupname = self.__get_name(TEMPLATE_PAR)

        self.parnames = []
        outdir = setupname+"_grid/"
        os.system("mkdir "+outdir)
        self.outdir = outdir
        
        simulations = open(outdir+setupname+".dat","w")
        combinations = self.__build_combinations(parameters)

        i = 0
        for combination in combinations:
            newpar = open(outdir+setupname+str(i)+".par","w")
            simulations.write("\n"+setupname+str(i)+"\n")
            self.parnames.append(outdir+setupname+str(i)+".par")
            for line in TEMPLATE_PAR:
                keyfound = False
                for key in combination:
                    if re.match(key.lower()+"\s",line.lower()):
                        simulations.write(key +"\t"+str(combination[key])+"\n")
                        newpar.write(key +"\t"+str(combination[key])+"\n")
                        keyfound=True
                        continue
                if re.match("outputdir\s",line.lower()):
                    newpar.write("outputdir\t@outputs/"+setupname+str(i)+"\n")
                    continue
                if keyfound:
                    continue
                newpar.write(line)
            newpar.close()
            i += 1
        simulations.close()

        for par in self.parnames:
            ifile = par[:-4]+".pbs"
            pbs = open(ifile,'w')
            for line in TEMPLATE_PBS:
                if re.search("%NAME",line):
                    line = line.replace("%NAME",par.split("/")[-1][:-4])
                    pbs.write(line)
                    continue
                if re.search("%PARFILE",line):
                    line = line.replace("%PARFILE","projects/"+setupname+"/"+par)
                    pbs.write(line)
                    continue
                pbs.write(line)
            pbs.close()
            
    def qsub(self):
        for f in os.listdir(self.outdir):
            if f.endswith(".pbs"):
                line = 'qsub '+ self.outdir + f
                print line
                os.system(line)
        
        
    def __build_combinations(self,Parameters):
        """
        Returns a list of dictionaries with the cartesian product of
        the keys of the dictionary Parameters.

        For example:

        Grid.__build_combinations({'a':[1,2],'b':[3,4]})
        [{'a':1,'b':3},{'a':1,'b':4},{'a':2,'b':3},{'a':2,'b':4}]

        and the same is valid for n dimensions.
        """

        combinations = []
        nparams = len(Parameters.keys()) #Number of parameters modified
        product_arguments = ''
        keys = []
        for key in Parameters:
            product_arguments += "Parameters['" + key + "'],"
            keys.append(key) #Taking care about randomly access...
        product_arguments = product_arguments[:-1]
        exec("iterator = itertools.product("+product_arguments+")")
        combinations = []
        for combination in iterator:
            dic = {}
            for i in range(len(keys)):
                dic[keys[i]] = combination[i]
            combinations.append(dic)
        return combinations
       
    def __get_name(self,lines):
        """
        Getting the setupname from the parfile. Useful for the name of
        the .par files.
        """

        for line in lines:
            match = re.match("setup\s*(\w*)",line.lower())
            if match:
                return match.group(1)

if __name__ == '__main__':

    Parameters = {}        
    Parameters['PRANDTL'] = np.linspace(0.5,3.0,16)
    g = Grid(Parameters)
    g.execute()
    g.qsub()
