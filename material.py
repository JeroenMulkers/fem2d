from numpy import array, cos, sin

class MatProp:

    def __init__(self,isLinear,isIsotropic,values):
        self.isLinear = isLinear
        self.isIsotropic = isIsotropic
        self.values = values

    def value(self,sol=None):
        if self.isLinear and self.isIsotropic:
            return array([[self.values],[self.values]])
        elif self.isLinear and not self.isIsotropic:
            return array([[self.values[0]],[self.values[1]]])

###############################################################################

def readMaterialFile(filename):
    material = {}
    fid = open(filename,'r')
    while True:
        line = fid.readline()
        if not line:
            break
        if len(line.split()) == 1:
            propname = line.split()[0]
            if propname in ["PERMEABILITY","CONDUCTIVITY"]:
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
                        values.append( fid.readline() )
                material[propname] = MatProp(isLinear,isIsotropic,values)
    fid.close()
    return material

names = ["ANIS","AIR","CU","FE","FLUI"]
materials = {}
for name in names:
    materials[name] = readMaterialFile("MATE/"+name+".MAT")
