import re
import numpy as np
import itertools
import os

class Grid():
    def __init__(self,parameters,template='template.par'):
        """
        Class for generating a set of cube of .par files. parameters
        is a dictionary with all the parameters you want to modify.
        template is a .par standard file taken as a model. Only will
        be modified the outputdir param, and all the parameters that
        matches with the keys of the dictionary parameters.
        """
        
        template = open(template)
        TEMPLATE = template.readlines()
        template.close()

        setupname = self.__get_name(TEMPLATE)

        self.parnames = []

        outdir = setupname+"_grid/"
        os.system("mkdir "+outdir)
        
        simulations = open(outdir+setupname+".dat","w")
        combinations = self.__build_combinations(parameters)

        i = 0
        for combination in combinations:
            newpar = open(outdir+setupname+str(i)+".par","w")
            simulations.write("\n"+setupname+str(i)+"\n")
            self.parnames.append(outdir+setupname+str(i)+".par")
            for line in TEMPLATE:
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
    Parameters['A'] = [-1.0,1.0]
    Parameters['B'] = [-0.3183098861837907]
    Parameters['ETA'] = np.linspace(5e-6,5e-3,50)
    g = Grid(Parameters)
