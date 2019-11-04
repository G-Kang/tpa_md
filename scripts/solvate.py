from bbb import lattice, block, packmol, randomize_orientation
import copy, sys, os
import math
import numpy
import scipy.special
import xml.dom.minidom
import operator
import hoomd2mol2

def SoakPolymer(initdir,polyfile,solvfile,outfile,filetype='xml',num_poly=1,num_solv=1,lf=50):
    solvfile += '.mol2'
    if num_poly==0:
        solvent=MOL2fileBlock(initdir,solvfile,mass=1.0)
        seed = numpy.random.randint(1,10000)
        box = packmol.ConstraintSimulationBoxHOOMD(1,1,1,tolerance=2.0)
        Lf = lf
        vol_box = Lf**3
        box.scaleToVolume(vol_box)
        #shell1 = packmol.ConstraintSphere("outside",0,0,0, NP_rad+2)
        #shell2 = packmol.ConstraintSphere("inside",0,0,0, 40)
        #shell3 = packmol.ConstraintSphere("outside",0,0,0, 40)
        #shell4 = packmol.ConstraintSphere("outside",0,0,0, Lf/2.0-4.0)
        #shell = packmol.ConstraintSphere("inside",0,0,0, Lf/2.0-5.0)
        origin = packmol.ConstraintFixed(0,0,0)
        if filetype=='xml':
            g = packmol.GeneratorXML(simulation_box=box, seed=seed, cleanup=False, tolerance=2.0)
            g.addBuildingBlock(solvent,num_solv)#,constraints=[shell])
            outfile += '.xml'
            g.writeOutput(open(initdir + '/' + outfile, 'w'))
        else:
            print box.getAABB()
            g = packmol.GeneratorTinker(simulation_box=box, seed=seed, filetype='tinker', cleanup=False, tolerance=2.0)
            g.addBuildingBlock(solvent,num_solv,tinkerxyz=initdir+'/'+solvfile.split('.')[0]+'.xyz')#,constraints=[shell])
            outfilexml =  outfile+'.xml'
            outfile += '.xyz'
            g.writeOutput(open(initdir + '/' + outfile, 'w'),open(initdir + '/' + outfilexml, 'w'))
        out_filemol2 = outfile.split('.')[0]+'.mol2'
        print outfile,initdir
        randomize_orientation.xml_to_pos(initdir+'/'+outfilexml,"atom")
        hoomd2mol2.writeout(initdir,outfilexml,out_filemol2)
    else:
        polyfile += '.mol2'
        DApoly=MOL2fileBlock(initdir,polyfile,mass=1.0)
        solvent=MOL2fileBlock(initdir,solvfile,mass=1.0)
        seed = numpy.random.randint(1,10000)
        box = packmol.ConstraintSimulationBoxHOOMD(1,1,1,tolerance=2.0)
        Lf = lf
        vol_box = Lf**3
        box.scaleToVolume(vol_box)
        #shell1 = packmol.ConstraintSphere("outside",0,0,0, NP_rad+2)
        #shell2 = packmol.ConstraintSphere("inside",0,0,0, 40)
        #shell3 = packmol.ConstraintSphere("outside",0,0,0, 40)
        #shell4 = packmol.ConstraintSphere("outside",0,0,0, Lf/2.0-4.0)
        shell = packmol.ConstraintSphere("inside",0,0,0, Lf/2.0-5.0)
        origin = packmol.ConstraintFixed(0,0,0)
        if filetype=='xml':
            if latflag == 1:
                g = lattice.GeneratorXML();
            else:
                g = packmol.GeneratorXML(simulation_box=box, seed=seed, cleanup=False, tolerance=2.0)
            g.addBuildingBlock(DApoly,num_poly,constraints=[origin])
            g.addBuildingBlock(solvent,num_solv,constraints=[shell])
            outfile += '.xml'
            g.writeOutput(open(initdir + '/' + outfile, 'w'))
        else:
            g = packmol.GeneratorTinker(simulation_box=box, seed=seed, filetype='tinker', cleanup=False, tolerance=2.0)
            g.addBuildingBlock(DApoly,num_poly,tinkerxyz=initdir+'/'+polyfile.split('.')[0]+'.xyz',constraints=[origin])
            g.addBuildingBlock(solvent,num_solv,tinkerxyz=initdir+'/'+solvfile.split('.')[0]+'.xyz',constraints=[shell])
            outfilexml =  outfile+'.xml'
            outfile += '.xyz'
            g.writeOutput(open(initdir + '/' + outfile, 'w'),open(initdir + '/' + outfilexml, 'w'))
        out_filemol2 = outfile.split('.')[0]+'.mol2'
        print outfile,initdir
        randomize_orientation.xml_to_pos(initdir+'/'+outfilexml,"atom")
        hoomd2mol2.writeout(initdir,outfilexml,out_filemol2)
#def addmonomer():

#def writepolymer():

class MOL2fileBlock(block.BuildingBlock):

    def __init__(self,initdir,file,d_atom=1.0,mass=None):
        # initialize the base class
        block.BuildingBlock.__init__(self);
        
        fileData = open(initdir + "/" + file, 'r').read()
        linesInFile = fileData.splitlines()
        vertices = []
        types = []
        bonds = []
        bondfl=0
        atomfl=0
        bondid='at'
        for line in linesInFile:
            if "ATOM" in line:
                atomfl=1
                continue
            elif len(line.split()) < 1:
                atomfl=0
                bondfl=0
                continue
            elif atomfl==1 and 'BOND' not in line:
                valx=float(line.split()[2])
                valy=float(line.split()[3])
                valz=float(line.split()[4])
                type=line.split()[1]
                coords = numpy.array((valx,valy,valz))
                types.append(type)
                vertices.append(coords)
            elif "BOND" in line:
                bondfl=1
                atomfl=0
                continue
            elif bondfl==1: 
                #bondty=str(line.split()[3])
                atA=int(line.split()[1])
                atB=int(line.split()[2])
                bond= [bondid,atA,atB]
                bonds.append(bond)
            else:
                continue
        xyzfile=initdir + "/" + file.split('.')[0]+'.xyz'
        print xyzfile
        xyz=os.path.isfile(xyzfile)
        if xyz==True:
            fileData = open(xyzfile, 'r').read()
            linesInFile = fileData.splitlines()
            fftyps = []
            for line in linesInFile[1:]:
                fftyp=int(line.split()[5])
                fftyps.append(fftyp)
            if len(fftyps) != len(types):
                print 'xyz file is not consistent with mol2 file'
                sys.exit(1)
        else:
            fftyps = []
            for typ in types:
                fftyps.append(0)
        # start iterating the particles
        for t,v,f in zip(types,vertices,fftyps):
            p = block.Particle(position=v,diameter=d_atom,type=t,body=-1,charge=f)
            self.addParticle(p);

        for b in bonds:
            self.addBond(block.Bond(b[0],b[1]-1,b[2]-1))
        # Set mass if specified 
        #print "done"
        #if not mass is None:
        #   if not mass > 0:
        #       raise ValueError('building block mass must be greater than 0');
        #   self.setMass(mass)
        

        centroid = self.getCentroid()
    #    d = block.Particle(position=centroid,diameter=0.5,type='D',body=0)
    #    self.addParticle(d)
    def getNumAtoms(self):
        max_beads = -1
        for b in self.iterParticles():
            max_beads += 1
        return max_beads+1
