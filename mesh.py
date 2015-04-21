import numpy as np
from numpy import array, cos, sin
from material import materials
from scipy.sparse import coo_matrix
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


###############################################################################

class Node:

    def __init__(n,id,x,y,boundary=False,value=None):
        n.id = id
        n.x = x
        n.y = y
        n.boundary = boundary
        n.value = value

###############################################################################

class Region:

    def __init__(self,id,material,orientation,source=0.0):
        self.id = id
        self.material = material
        self.orientation = orientation
        self.source = source

    def calcTensor(self,solution=None):
        a = self.orientation
        v = materials[self.material].permeability
        T = array([[cos(a),-sin(a)],[sin(a),cos(a)]])
        return T*v*T.transpose()

###############################################################################

class Element:

    def __init__(self,nodes,region=None):
        self.node = nodes
        self.region = region
        n1,n2,n3 = self.node
        self.area = (n1.x*(n2.y-n3.y)+n2.x*(n3.y-n1.y)+n3.x*(n1.y-n2.y))/2
        self.b = [n2.y-n3.y,n3.y-n1.y,n1.y-n2.y];
        self.c = [n3.x-n2.x,n1.x-n3.x,n2.x-n1.x];

    def calcTensor(self):
        return self.region.calcTensor()

    def calcSourceVec(self):
        return 3*[self.region.source*self.area/3]

    def calcDiffMat(self):
        k = self.calcTensor()
        self.Ke = np.zeros((3,3))
        for j in range(3):
            coefj = array([[self.b[j]],[self.c[j]]])
            for i in range(3):
                coefi = array([self.b[i],self.c[i]])
                self.Ke[j,i] = -coefi.dot(k).dot(coefj)/(4*self.area)
        return self.Ke

###############################################################################

class Mesh:

    def __init__(self,dirname):
        self.dirname = dirname
        self.readFiles()

    def solve(self):

        NN = len(self.node)
        rhs = np.zeros(NN)
        K = lil_matrix((NN,NN))

        # loop through all the elements
        for e in self.element:

            # initial Ke and fe
            Ke = e.calcDiffMat()
            Se = e.calcSourceVec()

            for i,ni in enumerate(e.node):
                rhs[ni.id] -= Se[i]
                if ni.boundary:
                    rhs[ni.id] = ni.value
                    K[ni.id,ni.id] = 1
                else:
                    for j,nj in enumerate(e.node):
                        if nj.boundary:
                            rhs[ni.id] -= Ke[i,j]*nj.value
                            K[ni.id,nj.id] = 0
                        else:
                            K[ni.id,nj.id] += Ke[i,j]

        # solve the system of linear equations and put solution on nodes
        solution = spsolve(K,rhs)
        for i in range(len(solution)):
            self.node[i].value = solution[i]

    def readFiles(self):

        import os

        self.node = []
        self.region = []
        self.element = []

        # read nodes
        filename = os.path.join(self.dirname,"nodes.txt")
        dtype = [('x',np.float64),('y',np.float64)]
        coors = np.loadtxt(filename,dtype=dtype)
        for id, c in enumerate(coors):
            self.node.append(Node(id,c['x'],c['y']))

        # read boundary values
        fvals = os.path.join(self.dirname,"uconsvals.txt")
        vals = np.loadtxt(fvals,unpack=True)[2]
        fboundelems = open(os.path.join(self.dirname,"ucons.txt"))
        for line in fboundelems:
            line = line.split()
            self.node[int(line[0])-1].boundary = True
            self.node[int(line[0])-1].value = vals[int(line[2])-1]
        fboundelems.close()

        # read regions (materials and orientations)
        forient = open(os.path.join(self.dirname,"orientations.txt"))
        fmaterial = open(os.path.join(self.dirname,"materials.txt"))
        fsources = open(os.path.join(self.dirname,"sources.txt"))
        for id, (orient, materialname, j) in enumerate(zip(forient,fmaterial,fsources)):
            materialname = materialname.split()[0]
            j = j.split()[0]
            self.region.append( Region(id,materialname,float(orient),float(j)) )
        forient.close()
        fmaterial.close()
        fsources.close()

        # read elements
        felements = open(os.path.join(self.dirname,"elems.txt"))
        for line in felements:
            line = line.split()
            nodes = [ self.node[int(line[i])-1] for i in range(3) ]
            region = self.region[int(line[3])-1]
            self.element.append( Element(nodes,region) )
        felements.close()

    #----------------------------------------------------------------------

    def readSolution(self):
        import os
        fsolution = open(os.path.join(self.dirname,"solu2.txt"))
        for node in self.node:
            node.value = float(fsolution.readline())
        fsolution.close()

    def getUnstructuredGrid(self):
        import vtk
        grid = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(len(self.node))
        for n in self.node:
            points.InsertPoint(n.id,n.x,n.y,0.0)
        grid.SetPoints(points)
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(grid)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        return grid

    def writeVTKfile(self):

        NP = len(self.node)
        NE = len(self.element)

        filename = "output.vtk"
        fid = open(filename,"w")

        fid.write("# vtk DataFile Version 1.0\n")
        fid.write("title\n")
        fid.write("ASCII\n")
        fid.write("DATASET UNSTRUCTURED_GRID\n")

        fid.write("POINTS %d float\n"%NP)
        for n in self.node:
            fid.write("%f %f %f\n"%(n.x,n.y,0.0))

        fid.write("CELLS %d %d\n"%(NE,4*NE))
        for e in self.element:
            fid.write("%d %d %d %d\n"%(3,e.node[0].id,e.node[1].id,e.node[2].id))

        fid.write("CELL_TYPES %d\n"%NE)
        for i in range(NE):
            fid.write("5\n")

        fid.write("POINT_DATA %d\n"%NP)

        fid.write("SCALARS solution float\n")
        fid.write("LOOKUP_TABLE default\n")
        for n in self.node:
            fid.write("%f\n"%n.value)
        fid.close()

###############################################################################

if __name__ == "__main__":
    mesh = Mesh("./Prob1_thermal_convdiff/")
    mesh.solve()
    mesh.writeVTKfile()
