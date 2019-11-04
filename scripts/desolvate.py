from sys import argv
import numpy as np
import numpy.linalg as linalg
from numpy.matlib import *
from scipy import spatial
import itertools
from collections import OrderedDict
import operator
from solv_helper import *

script,path,file,sol_size,shell_type,solv_cut=argv
cutoff=2.5
sol_size=int(sol_size)
polyfile=file.split('min')[0]+'min.xyz'
#print connectivity
atom_array,tink_type,atom_type,connect,slines,natoms,cent_size=read_xyz(path,file,polyfile,sol_size)
shell_atoms,all_atoms,cent_molecule=solv_shell(atom_array,connect,natoms,cent_size,sol_size,cutoff)

#Centroid of each solvent mol
clean_solv_shell,clean_solv_all=deloop_solv(sol_size,shell_atoms,all_atoms,atom_array,cent_molecule,solv_cut)
if shell_type=='all':
    write_xyz(path,file,atom_array,atom_type,tink_type,connect,slines,cent_molecule,clean_solv_all,shell_type)
else:
    write_xyz(path,file,atom_array,atom_type,tink_type,connect,slines,cent_molecule,clean_solv_shell,shell_type)
####
#print new_solv

