from bbb import lattice, block, packmol, randomize_orientation
import copy, sys, os
import math
import numpy
import scipy.special
import xml.dom.minidom
import operator
import hoomd2mol2

def CreatePolymer(initdir,file,outfile,filetype='xml',num_mono=1,num_poly=1,lf=50,tw=0,latflag=0):
    file += '.mol2'
    DApoly=MOL2fileBlock(initdir,file,mass=1.0)
    DAmono=copy.deepcopy(DApoly)
    m=num_mono-1
    if m==0:
        seed = numpy.random.randint(1,10000)
        box = packmol.ConstraintSimulationBoxHOOMD(1,1,1,tolerance=2.0)
        Lf = lf
        vol_box = Lf**3
        box.scaleToVolume(vol_box)
        #shell1 = packmol.ConstraintSphere("outside",0,0,0, NP_rad+2)
        #shell2 = packmol.ConstraintSphere("inside",0,0,0, 40)
        #shell3 = packmol.ConstraintSphere("outside",0,0,0, 40)
        #shell4 = packmol.ConstraintSphere("outside",0,0,0, Lf/2.0-4.0)
        #shell5 = packmol.ConstraintSphere("inside",0,0,0, Lf/2.0-2.0)
        #origin = packmol.ConstraintFixed(0,0,0)
        if filetype=='xml':
            gx = packmol.GeneratorXML(simulation_box=box, seed=seed, cleanup=False, tolerance=2.0)
            gx.addBuildingBlock(DApoly,num_poly)
            outfile += '.xml'
            gx.copyOutput(open(initdir + '/' + outfile, 'w'))
            randomize_orientation.xml_to_pos(initdir+'/'+outfile,"atom")
            out_filemol2 = outfile.split('.')[0]+'.mol2'
            hoomd2mol2.writeout(initdir,outfile,out_filemol2)
            print( initdir+'/'+outfile)
        elif filetype=='xyz':
            gt = packmol.GeneratorTinker(simulation_box=box, seed=seed, cleanup=False, tolerance=2.0)
            gt.addBuildingBlock(DApoly,num_poly)
            outfile += '.xyz'
            gt.writeOutput(open(initdir + '/' + outfile, 'w'))
            gx = packmol.GeneratorXML(simulation_box=box, seed=seed, cleanup=False, tolerance=2.0)
            gx.addBuildingBlock(DApoly,num_poly)
            print( initdir+'/'+outfile)
            outfile[:-4] += '.xml'
            gx.copyOutput(open(initdir + '/' + outfile, 'w'))
            randomize_orientation.xml_to_pos(initdir+'/'+outfile,"atom")
            print( initdir+'/'+outfile)
    else:
        Ai=[]
        Di=[]
        for p in DAmono.iterParticles():
            if p.type == 'Ca':
                acAt= p.position
                break
        nextCD=0
        for p in DAmono.iterParticles():
            if p.type == 'Cd' and nextCD==0:
                nextCD=1
                continue
            if p.type == 'Cd' and nextCD==1:
                dnAt= p.position
                nextCD=0
                break
        off=acAt-dnAt
        odis = numpy.sqrt(off[0]**2 + off[1]**2 + off[2]**2)
        px,py,pz= off[0]/odis,off[1]/odis,off[2]/odis
        j=0
        while m > 0:
            for i,p in enumerate(DApoly.iterParticles()):
                if p.type == 'Ca':
                    Ai.append(i)
            for i,p in enumerate(DAmono.iterParticles()):
                if p.type == 'Cd':
                    Di.append(i)

            if j==0 and tw==0:
                up=(0,1,0)
            elif j==0 and tw==1:
                up=(0,0,1)
                DApoly.reorient((0,.2,-1))
                j+=1
            elif j==1 and tw==1:
                up=(0,0,-1)
                DApoly.reorient((0,1,0.2))
                j+=1
            elif j== 2 and tw==1:
                up=(0,1,0)
                DApoly.reorient((0,.2,-1))
                j=0
            else:
                print( "Epic Fail")
                sys.exit(1)
            DApoly.attach(DAmono,Ai[-1],Di[0],offset=(1.4*px,1.4*py,1.4*pz),new_up=up, bond_type='DA')
            Ai=[]
            Di=[]
            m-=1
        #add terminal H's
        else:
            DH = block.BuildingBlock()
            AH = block.BuildingBlock()
            DHp = block.Particle(position=(0,0,0),diameter=0.5,type='H',body=-1,charge=88)
            AHp = block.Particle(position=(0,0,0),diameter=0.5,type='H',body=-1,charge=88)
            DH.addParticle(DHp);
            AH.addParticle(DHp);
            for i,p in enumerate(DApoly.iterParticles()):
                if p.type == 'Ca':
                    Ai.append(i)
                elif p.type == 'Cd':
                    Di.append(i)
            DApoly.attach(DH,Di[0],0,offset=(-1.08*px,-1.08*py,-1.08*pz),new_up=(0,1,0), bond_type='DA')
            DApoly.attach(AH,Ai[-1],0,offset=(1.08*px,1.08*py,1.08*pz),new_up=(0,1,0), bond_type='DA')


        if filetype=='xml':
            if latflag == 1:
                g = lattice.GeneratorXML();
            else:
                seed = numpy.random.randint(1,10000)
                box = packmol.ConstraintSimulationBoxHOOMD(1,1,1,tolerance=2.0)
                Lf = lf
                vol_box = Lf**3
                box.scaleToVolume(vol_box)
                #shell1 = packmol.ConstraintSphere("outside",0,0,0, NP_rad+2)
                #shell2 = packmol.ConstraintSphere("inside",0,0,0, 40)
                #shell3 = packmol.ConstraintSphere("outside",0,0,0, 40)
                #shell4 = packmol.ConstraintSphere("outside",0,0,0, Lf/2.0-4.0)
                #shell5 = packmol.ConstraintSphere("inside",0,0,0, Lf/2.0-2.0)
                #origin = packmol.ConstraintFixed(0,0,0)
                g = packmol.GeneratorXML(simulation_box=box, seed=seed, cleanup=False, tolerance=2.0)
            g.addBuildingBlock(DApoly,num_poly)
            #if num_Ca > 0:
            #    g.addBuildingBlock(Caions,num_Ca,constraints=[shell1,shell2])
            #g.addBuildingBlock(Clions,num_Cl,constraints=[shell1,shell5])
            #if num_FDNA>0:
            #    g.addBuildingBlock(FDNA,num_FDNA,constraints=[shell4,shell5])
            #g.addBuildingBlock(TNP,1,constraints=[origin]);
            outfile += '.xml'
            out_filemol2 = outfile.split('.')[0]+'.mol2'
            print( outfile,initdir)
            g.writeOutput(open(initdir + '/' + outfile, 'w'))
            randomize_orientation.xml_to_pos(initdir+'/'+outfile,"atom")
            hoomd2mol2.writeout(initdir,outfile,out_filemol2)
        else:
            seed = numpy.random.randint(1,10000)
            box = packmol.ConstraintSimulationBoxHOOMD(1,1,1,tolerance=2.0)
            Lf = lf
            vol_box = Lf**3
            box.scaleToVolume(vol_box)
            #shell1 = packmol.ConstraintSphere("outside",0,0,0, NP_rad+2)
            #shell2 = packmol.ConstraintSphere("inside",0,0,0, 40)
            #shell3 = packmol.ConstraintSphere("outside",0,0,0, 40)
            #shell4 = packmol.ConstraintSphere("outside",0,0,0, Lf/2.0-4.0)
            #shell5 = packmol.ConstraintSphere("inside",0,0,0, Lf/2.0-2.0)
            #origin = packmol.ConstraintFixed(0,0,0)
            g = packmol.GeneratorTinker(simulation_box=box, seed=seed, cleanup=False, tolerance=2.0)
            g.addBuildingBlock(DApoly,num_poly)
            #if num_Ca > 0:
            #    g.addBuildingBlock(Caions,num_Ca,constraints=[shell1,shell2])
            #g.addBuildingBlock(Clions,num_Cl,constraints=[shell1,shell5])
            #if num_FDNA>0:
            #    g.addBuildingBlock(FDNA,num_FDNA,constraints=[shell4,shell5])
            #g.addBuildingBlock(TNP,1,constraints=[origin]);
            outfile += '.xyz'
            out_filemol2 = outfile.split('.')[0]+'.mol2'
            print( outfile,initdir)
            g.writeOutput(open(initdir + '/' + outfile, 'w'))
            randomize_orientation.xml_to_pos(initdir+'/'+outfile,"atom")
            hoomd2mol2.writeout(initdir,outfile,out_filemol2)
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
        print( xyzfile)
        xyz=os.path.isfile(xyzfile)
        if xyz==True:
            fileData = open(xyzfile, 'r').read()
            linesInFile = fileData.splitlines()
            fftyps = []
            for line in linesInFile[1:]:
                fftyp=int(line.split()[5])
                fftyps.append(fftyp)
            if len(fftyps) != len(types):
                print( 'xyz file is not consistent with mol2 file')
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
        #print( "done")
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
