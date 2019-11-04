from polymerize import CreatePolymer
from shutil import copyfile
from os.path import isdir,isfile,join
from os import mkdir,listdir,system
import sys
#number of total monomers; always 1 for molecule
num_mono=1
#number of (identical) total molecules
num_poly=1
#box size in Angstroms
lf=70
#twist monomers around dihedral; Keep at 0 for molecules
twist=0
#place molecules in a lattice configuration
lattice=0
#directory of mol2 files
path='molecules'
#mol2list=[f for f in listdir(path) if (isfile(join(path,f)) and 'mol2' in f)]
mol2list=['diphenylethylene.mol2','diphenylstilbenen2.mol2','stilbenen2.mol2']
#independent trajectories
num_traj=1
#loop over all molecules of interest
for molfile in mol2list:
    prmfile= molfile.split('.')[0]+'o.prm'
    monofile= molfile.split('.')[0]+'_mono'
    if isfile(join(path,monofile)+'.xyz') == False:
        system('mol2xyz '+path+'/'+molfile) ##sybylxyz
        copyfile(path+'/'+molfile.split('.')[0]+'.xyz_2',path+'/'+monofile+'.xyz')
        system('rm '+path+'/'+molfile.split('.')[0]+'*.xyz_*')
        system('xyzmol2 '+path+'/'+monofile+'.xyz')
        continue
    #number of initial files
    else:
        for i in range(0,num_traj):
            polyfile=molfile.split('.')[0]+'_'+str(i)+'_mols'+str(num_poly)
            system('rm '+path+'/'+polyfile+'*.mol2*')
            system('rm '+path+'/'+polyfile+'*.xyz*')
            CreatePolymer(path,monofile,polyfile,filetype='xml',num_mono=num_mono,num_poly=num_poly,lf=lf,latflag=lattice)
            system('mol2xyz '+path+'/'+polyfile+'.mol2')
            system('xyzedit '+path+'/'+polyfile+'.xyz '+path+'/'+prmfile+' 12')
            copyfile(path+'/'+polyfile+'.xyz_2',path+'/'+polyfile+'_cent.xyz')
            system('rm '+path+'/'+polyfile+'.xyz_*')
            system('analyze '+path+'/'+polyfile+'_cent.xyz '+path+'/'+prmfile+' M')
            system('minimize '+path+'/'+polyfile+'_cent.xyz '+path+'/'+prmfile+' 0.2 > tmp.out')
            copyfile(path+'/'+polyfile+'_cent.xyz_2',path+'/'+polyfile+'_cent_min.xyz')
            system('rm '+path+'/'+polyfile+'_cent.xyz_*')
            system('xyzmol2 '+path+'/'+polyfile+'_cent_min.xyz')
            CreatePolymer(path,polyfile+'_cent_min',polyfile+'_cent_min',filetype='xml',num_mono=num_mono,num_poly=num_poly,lf=lf,latflag=lattice)
            print(i)
