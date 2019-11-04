from polymerize import CreatePolymer
import os,shutil,sys
import xml.dom.minidom
from solv_helper import *
#place polymer in a lattice configuration
path=sys.argv[1]
num=int(sys.argv[2])
size=int(sys.argv[3])
runid=int(sys.argv[4])
archive=sys.argv[5]
infile=sys.argv[6]
shell=sys.argv[7]
prmfile='D'+str(num)+'o'
arc=archive

def arctopos(path,file,shell,lx=65):
    posfile=file.split('.')[0]+'arc.pos'
    polyfile=file.split('min')[0]+'min.xyz'
    os.system('rm '+path+'/'+posfile)
    if shell=="True":
        if 'cbsol' in file:
            sol_size=12
        elif 'cfsol' in file:
            sol_size=5
        shell_type='shell'
        posfile=file.split('.')[0]+'shellarc.pos'
        solv_cut=1.0
        cutoff=8.0
    with open(path+'/'+file+'.arc') as infile:
        i=0
        j=0
        infile.seek(-10,os.SEEK_END)
        size=infile.tell()
        print size
        infile.seek(0,0)
        for line in infile:
            if j==0:
                tmpf='temp_'+str(num)+'_'+str(runid)
                os.system('rm '+path+'/'+tmpf+'*')
                f=open(path+'/'+tmpf+'.xyz','w')
                numatom=int(line.split()[0])
                xyzlines=line
                j+=1
            elif j==numatom:
                xyzlines+=line
                print infile.tell()
                f.write(xyzlines)
                f.close()
                #add
                if shell=="True":
                    atom_array,tink_type,atom_type,connect,slines,natoms,cent_size=read_xyz(path,tmpf+'.xyz',polyfile,sol_size)
                    shell_atoms,all_atoms,cent_molecule=solv_shell(atom_array,connect,natoms,cent_size,sol_size,cutoff)
                    clean_solv_shell,clean_solv_all=deloop_solv(sol_size,shell_atoms,all_atoms,atom_array,cent_molecule,solv_cut)
                    write_xyz(path,tmpf+'.xyz',atom_array,atom_type,tink_type,connect,slines,cent_molecule,clean_solv_shell,shell_type)
                    shfile=tmpf.split('.')[0]+'_shell'
                    s=open(path+'/'+shfile+'.xyz','r')
                    numatoms=int(s.readline().split()[0])
                    s.close()
                    os.system('xyzsybyl '+path+'/'+shfile+'.xyz > garb.out')
                    CreatePolymer(path,shfile,shfile,num_mono=1,num_poly=1,lf=lx)
                    if not os.path.isfile(path+'/'+posfile):
                        shutil.copyfile(path+'/'+shfile+'.pos',path+'/'+posfile)
                        os.system('rm '+path+'/'+shfile+'.*')
                        os.system('rm '+path+'/'+tmpf+'.*')
                    else:
                        os.system('cat '+path+'/'+shfile+'.pos >> '+path+'/'+posfile)
                        os.system('rm '+path+'/'+shfile+'.*')
                        os.system('rm '+path+'/'+tmpf+'.*')
                    if infile.tell()>size:
                        os.system('rm '+path+'/'+file+'_lastshell.*')
                        final=open(path+'/'+file+'_lastshell.xyz','w')
                        final.write(xyzlines)
                        final.close()
                        os.system('xyzsybyl '+path+'/'+file+'_lastshell.xyz')
                        CreatePolymer(path,file+'_lastshell',file+'_lastshell',num_mono=1,num_poly=1,lf=lx)
                else:
                    os.system('xyzsybyl '+path+'/'+tmpf+'>tmp.out')
                    CreatePolymer(path,tmpf,tmpf,num_mono=1,num_poly=1,lf=lx)
                    if not os.path.isfile(path+'/'+posfile):
                        shutil.copyfile(path+'/'+tmpf+'.pos',path+'/'+posfile)
                        os.system('rm '+path+'/'+tmpf+'.*')
                    else:
                        os.system('cat '+path+'/'+tmpf+'.pos >> '+path+'/'+posfile)
                        os.system('rm '+path+'/'+tmpf+'.*')
                    if infile.tell()>size:
                        os.system('rm '+path+'/'+file+'_last.*')
                        final=open(path+'/'+file+'_last.xyz','w')
                        final.write(xyzlines)
                        final.close()
                        os.system('xyzsybyl '+path+'/'+file+'_last.xyz')
                        CreatePolymer(path,file+'_last',file+'_last',num_mono=1,num_poly=1,lf=lx)
                j=0
                i+=1
            else:
                xyzlines+=line
                j+=1
                
if arc=='True':
    CreatePolymer(path,infile,infile,num_mono=1,num_poly=1,lf=65)
    #number of initial files
    dom = xml.dom.minidom.parse(path+'/'+infile+'.xml');

    # start by parsing the file
    hoomd_xml = dom.getElementsByTagName('hoomd_xml')[0];
    configuration = hoomd_xml.getElementsByTagName('configuration')[0];

    # read the box size
    box     = configuration.getElementsByTagName('box')[0];
    Lx = box.getAttribute('lx');
    l=float(Lx)
    arctopos(path,infile,shell,lx=l)
else:
    infile='D'+num+'_'+tornum+'_t_'
    outfile=infile
    tornum=sys.argv[6]
    for j in xrange(0,190,10):
        #os.system('xyzsybyl '+path+'/'+infile+'.001_'+str(j))
        os.system('analyze '+path+'/'+infile+str(j)+' '+path+'/../'+prmfile+' E > D'+num+'_'+tornum+'tor'+str(j)+'.out')
        #os.system('sybylxyz '+path+'/'+infile+str(j)+'_mp2')
        print path+'/'+infile+str(j)
        #CreatePolymer(path,infile+str(j),outfile+str(j),num_mono,num_poly,lf,tw=twist,latflag=lattice)
