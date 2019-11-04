from solvate import SoakPolymer
import os,shutil
#number of (identical) total polymer chains
num_poly=0
#box size in Angstroms
lf=30
#number of solvent molecules
num_solv=204  #30 A box 204 cf and 160 cb
path='init/solvent'
solvfile='cfsol'
if num_poly>0:
    j=1
    i=0
    num_mono=4
    polyfile='D'+str(j)+'_'+str(i)+'_poly'+str(num_mono)
    outfile=polyfile+'_'+solvfile
else:
    polyfile=''
    outfile=solvfile+'_box'
    
SoakPolymer(path,polyfile,solvfile,outfile,filetype='xyz',num_poly=num_poly,lf=lf,num_solv=num_solv)
