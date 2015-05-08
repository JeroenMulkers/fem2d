
###############################################################################
#
# This module creates a material library based on the material files in ./MATE
#
# The material library is a dictionary of materials which are dictionaries
# of properties
#
# This should be used as follows:
#
#    from material import materials
#    materials['materialname']['propertyname'].value(solution)
#    materials['materialname']['propertyname'].isIsotropic
#    materials['materialname']['propertyname'].isLinear
#
###############################################################################

from scipy import array

###############################################################################

class MatProp:

    """ Material property class """

    def __init__(self,isLinear,isIsotropic,values):
        self.isLinear = isLinear
        self.isIsotropic = isIsotropic
        self.values = values

    def value(self,sol=None):

        """ returns the value which may depend on the solution """

        if self.isLinear and self.isIsotropic:
            return array([[self.values],[self.values]])
        elif self.isLinear and not self.isIsotropic:
            return array([[self.values[0]],[self.values[1]]])
        elif not self.isLinear and self.isIsotropic:
            for x,y in self.values:
                if x>sol: break
                xprev, yprev = x,y
            return yprev+(x-sol)*(y-yprev)/(x-xprev) # linear interpolation

###############################################################################

def readMaterialFile(filename):

    """ reads a .MAT file and returns a material """

    material = {}
    fid = open(filename,'r')
    while True:
        line = fid.readline()
        if not line:
            break
        if len(line.split()) == 1:
            propname = line.split()[0]
            if propname in ["PERMEABILITY","CONDUCTIVITY","DIFFUSIVITY"]:
                line1 = fid.readline().split()[0]
                line2 = fid.readline().split()[0]
                isLinear = True if (line1=="linear") else False
                isIsotropic = True if (line2=="isotropic") else False
                if (isLinear and isIsotropic):
                    values = float(fid.readline().split()[0])
                elif (isLinear and not isIsotropic):
                    value_x = float(fid.readline().split()[0])
                    value_y = float(fid.readline().split()[0])
                    values = [value_x,value_y]
                else:
                    values = []
                    nDataPoints = int(fid.readline().split()[0])
                    for i in range(nDataPoints):
                        values.append([ float(v) for v in fid.readline().split() ])
                material[propname] = MatProp(isLinear,isIsotropic,values)
    fid.close()
    return material

###############################################################################
# building of the material library

names = ["ANIS","AIR","CU","FE","FLUI"]
materials = {}
for name in names:
    materials[name] = readMaterialFile("MATE/"+name+".MAT")
