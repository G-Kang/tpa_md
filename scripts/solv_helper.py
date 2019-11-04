import numpy as np
import numpy.linalg as linalg
from numpy.matlib import *
from scipy import spatial
import itertools
from collections import OrderedDict
import operator,sys

def read_xyz(path,orig_file,polyfile,solv_size):
    initpath='init/solvent'
    cent_size = int(open(initpath+'/'+polyfile).readlines()[0].split()[0])

    with open(path+'/'+orig_file, 'r') as xyzf:
        Masterlines = xyzf.readlines()
        solvlines = Masterlines[-int(solv_size):]
        num_atoms = int(Masterlines[0].split()[0])
        Masterlines = Masterlines[1:]
    xyzf.close()
    connectivity = []
    atom_array = np.zeros( (num_atoms,3), dtype = float)
    atom_type = []
    tink_type = []
    for i,line in enumerate(Masterlines):
        temp = line.split()
        atom_array[i,0] = float(temp[2])
        atom_array[i,1] = float(temp[3])
        atom_array[i,2] = float(temp[4])
        tink_type.append(int(temp[5]))
        connectivity.append(map(int,temp[6:]))
        atom_type.append(temp[1])
    return atom_array,tink_type,atom_type,connectivity,solvlines,num_atoms,cent_size

def solv_shell(atom_array,connectivity,num_atoms,cent_size,solv_size,cutoff):
    #cutoff = 2.5
    solv_size = int(solv_size)
    solvent_shell = []

    atom_indx=np.arange(num_atoms+1)[1:].reshape([1,num_atoms])
    num_catoms=cent_size
    cent_molecule=np.arange(num_catoms+1)[1:]
    all_sol_atoms = atom_indx[atom_indx > len(cent_molecule)]
    for item in cent_molecule:
        XA_array = np.array([[atom_array[item-1,0],atom_array[item-1,1],atom_array[item-1,2]]])
        cdistarray = spatial.distance.cdist(XA_array,atom_array)
        shellarray = atom_indx[cdistarray < cutoff]
        sharrwocent = shellarray[shellarray > len(cent_molecule)]
        solvent_shell.append(sharrwocent)
    #remove duplicates from solvent_shell
    if len(list(itertools.chain.from_iterable(solvent_shell)))<1:
        print "no molecules in shell exiting"
        sys.exit(0)
    solvent_shell=sorted(list(OrderedDict.fromkeys(itertools.chain.from_iterable(solvent_shell))))
    all_solvent=sorted(list(OrderedDict.fromkeys(all_sol_atoms)))
    #print cent_molecule
    #print solvent_shell[0:40]
    for item in solvent_shell:
        for item2 in connectivity[item-1]:
            if item2 not in solvent_shell:
                solvent_shell.append(item2)
    return sorted(solvent_shell),all_solvent,cent_molecule

def deloop_solv(mol_size,shell_atominds,all_atominds,atom_array,cent_molecule,solv_cut):
    solv_cent = []
    solv_map = []
    for j,i in enumerate(shell_atominds[::mol_size]):
        fatm = i-1
        latm = fatm+mol_size
        xm,ym,zm = atom_array[fatm:latm,0],atom_array[fatm:latm,1],atom_array[fatm:latm,2]
        #if j < 2:
        #    print fatm,latm
        #    print xm,ym,zm
        #else:
        #    sys.exit(0)
        xc,yc,zc = reduce(operator.add, xm, 0)/float(mol_size),reduce(operator.add, ym, 0)/float(mol_size), reduce(operator.add, zm, 0)/float(mol_size)
        solv_cent.append([xc,yc,zc])
        solv_map.append([fatm,latm])
    #cf
    ccutoff=float(solv_cut)
    #ccutoff=3.0
    #cb
    #ccutoff=3.6
    bad_solv=[]
    for item in cent_molecule:
            XA_array = np.array([[atom_array[item-1,0],atom_array[item-1,1],atom_array[item-1,2]]])
            cdistarray = spatial.distance.cdist(XA_array,np.array(solv_cent))
            for i,dist in enumerate(cdistarray[0]):
                if dist < ccutoff:
                    bad_solv.append(i)
    bad_solv=sorted(list(set(bad_solv)))
    print "deleting ",len(bad_solv)," solvent molecules"
    for i in bad_solv:
        fatm = shell_atominds.index(solv_map[i][0]+1)
        fatm1 = all_atominds.index(solv_map[i][0]+1)
        latm = shell_atominds.index(solv_map[i][1])
        latm1 = all_atominds.index(solv_map[i][1])
        #print fatm,fatm1
        #print latm,latm1
        del all_atominds[fatm1:latm1+1]
        del shell_atominds[fatm:latm+1]
    return sorted(shell_atominds),sorted(all_atominds)

def write_xyz(path,orig_file,atom_array,atom_type,tink_type,connectivity,solvlines,cent_molecule,solv_atoms,out_type):
    solv_connect = []
    solv_size = len(solvlines)
    for j,line in enumerate(solvlines):
        temp = line.split()
        solv_connect.append([i-int(temp[0]) for i in map(int,temp[6:])])
    #print len(solv_atoms),len(cent_molecule)
    cent_str=np.zeros([len(cent_molecule)],dtype=object)
    for indx,item in enumerate(cent_molecule):
        cstr=str(indx+1)+"\t"+str(atom_type[int(item)-1])+"\t"+str(atom_array[int(item)-1,0])+"\t"+str(atom_array[int(item)-1,1])+"\t"+str(atom_array[int(item)-1,2])+"\t"
        cstr=cstr+str(tink_type[int(item)-1])+"\t"+' '.join([str(i) for i in connectivity[int(item)-1]])
        cent_str[indx]=cstr
    cent=indx+1
    sindx=0
    solv_str=np.zeros([len(solv_atoms)],dtype=object)
    for indx,item in enumerate(solv_atoms):
        final_str=str(cent+indx+1)+"\t"+str(atom_type[int(item)-1])+"\t"+str(atom_array[int(item)-1,0])+"\t"+str(atom_array[int(item)-1,1])+"\t"+str(atom_array[int(item)-1,2])+"\t"
        final_str=final_str+str(tink_type[int(item)-1])+"\t"+' '.join([str(cent+indx+1+i) for i in solv_connect[sindx]])
        solv_str[indx]=final_str
        if indx == 0 or sindx==0:
            sindx+=1
        elif (sindx % (solv_size-1))==0:
            sindx=0
        else:
            sindx+=1
    final_file = orig_file.split('.')[0]+'_'+out_type+'.xyz'
    final = open(path+'/'+final_file,'w')
    final.write(str(len(solv_atoms)+len(cent_molecule))+"\n")
    cent_str.tofile(final,"\n")
    final.write("\n")
    solv_str.tofile(final,"\n")

    final.close()
