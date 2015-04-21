from numpy import array, cos, sin

class material:
    def __init__(self,name,filename):
        self.name = name
        fid = open(filename,'r')
        self.readMaterialFile(fid)
        fid.close()

    def readMaterialFile(self,fid):
        while True:
            line = fid.readline()
            if not line:
                break
            line = line.split()
            if len(line)==0:
                continue
            elif line[0] == "PERMEABILITY":
                self.readPermeability(fid)
            elif line[0] == "CONDUCTIVITY":
                self.readConductivity(fid)


    def readPermeability(self,fid):
        line1 = fid.readline().split()[0]
        line2 = fid.readline().split()[0]

        self.permeabilityLinear = True if (line1=="linear") else False
        self.permeabilityIsotropic = True if (line2=="isotropic") else False

        if (self.permeabilityLinear and self.permeabilityIsotropic):
            pxx = float(fid.readline().split()[0])
            pyy = pxx
            self.permeability = array([[pxx],[pyy]])
        elif (self.permeabilityLinear and not self.permeabilityIsotropic):
            pxx = float(fid.readline().split()[0])
            pyy = float(fid.readline().split()[0])
            self.permeability = array([[pxx],[pyy]])
        else:
            pass

    def readConductivity(self,fid):
        line1 = fid.readline().split()[0]
        line2 = fid.readline().split()[0]

        self.conductivityLinear = True if (line1=="linear") else False
        self.conductivityIsotropic = True if (line2=="isotropic") else False

        if (self.conductivityLinear and self.conductivityIsotropic):
            self.conductivity = float(fid.readline().split()[0])
        elif (self.conductivityLinear and not self.conductivityIsotropic):
            self.conductivity = [None,None]
            self.conductivity[0] = float(fid.readline().split()[0])
            self.conductivity[1] = float(fid.readline().split()[0])
        else:
            pass

names = ["ANIS","AIR","CU","FE","FLUI"]
materials = {}
for name in names:
    materials[name] = material(name,"MATE/"+name+".MAT")
