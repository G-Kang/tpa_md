import os,sys
from solv_helper import *
from polymerize import CreatePolymer
def arctonrg(path,initpath,file,astr,period=10,archive=False,shell=False):
    polyfile=file.split('min')[0]+'min.xyz'
    num=int(file.split('_')[0][1:])
    #num=0
    runid=file.split('_')[1]
    #runid=0
    prmfile='_'.join(file.split('_')[0:4:2])+'.key'
    #prmfile=file.split('_')[0]+'.key'
    print path+'/'+file
    if shell==True:
        if 'cbsol' in file:
            sol_size=12
        elif 'cfsol' in file:
            sol_size=5
        shell_type='shell'
        solv_cut=1.0
        cutoff=8.0
    if archive==True:
        with open(path+'/'+file, 'r') as xyzf:
            numatoms=int(xyzf.readline().split()[0])
            linereset=(numatoms+1)*period
            lind=0
            i=0
            j=0
            xyzf.seek(-10,os.SEEK_END)
            size=xyzf.tell()
            print size
            xyzf.seek(0,0)
            Rg=[]
            En=[]
            Vdw=[]
            El=[]
            Ei=[]
            keylines=open(initpath+'/'+prmfile,'r').read()
            for i,line in enumerate(xyzf):
                if (i % linereset) == 0:
                    if i == 0:
                        j=0
                        tmpf='temp_'+str(num)+'_'+str(runid)
                        #os.system('rm '+path+'/'+tmpf+'*')
                        f=open(path+'/'+tmpf+'.xyz','w')
                        numatom=int(line.split()[0])
                        xyzlines=line
                        newfile=0
                    else:
                        #print i
                        print xyzf.tell()
                        f.write(xyzlines)
                        f.close()
                        if shell==True:
                            atom_array,tink_type,atom_type,connect,slines,natoms,cent_size=read_xyz(path,tmpf+'.xyz',polyfile,sol_size)
                            shell_atoms,all_atoms,cent_molecule=solv_shell(atom_array,connect,natoms,cent_size,sol_size,cutoff)
                            clean_solv_shell,clean_solv_all=deloop_solv(sol_size,shell_atoms,all_atoms,atom_array,cent_molecule,solv_cut)
                            write_xyz(path,tmpf+'.xyz',atom_array,atom_type,tink_type,connect,slines,cent_molecule,clean_solv_shell,shell_type)
                            infile=tmpf.split('.')[0]+'_shell'
                            s=open(path+'/'+infile+'.xyz','r')
                            numatoms=int(s.readline().split()[0])
                            s.close()
                            os.system('xyzsybyl '+path+'/'+infile+'.xyz > garb.out')
                            if num==1 or num==3 or num==6 or num==7 or num==8 or num==14 or num==15:
                                lf=65
                            else:
                                lf=100
                            CreatePolymer(path,infile,infile,num_mono=1,num_poly=1,lf=lf)
                        else:
                            infile=tmpf
                        xyzlines=line
                        numatom=int(line.split()[0])
                        for group in ['t','solv','poly']:
                            numpolyatms=int(open(initpath+'/'+polyfile).readline().split()[0])
                            if group=='solv':
                                open(initpath+'/'+prmfile+'.tmp','w').write(keylines+"\ninactive -1 "+str(numpolyatms))
                            elif group=='poly':
                                open(initpath+'/'+prmfile+'.tmp','w').write(keylines+"\ninactive -"+str(numpolyatms+1)+" "+str(numatoms))
                            else:
                                open(initpath+'/'+prmfile+'.tmp','w').write(keylines)
                            for item in astr.split():
                                os.system('analyze -k '+initpath+'/'+prmfile+'.tmp'+' '+path+'/'+infile+' '+item+'>'+tmpf+'mtmp.out')
                                if item=="M":
                                    Rglines=open(tmpf+'mtmp.out','r').readlines()
                                    for line in Rglines:
                                        if "Radius of Gyration" in line:
                                            Rg.append(float(line.split()[4]))
                                elif item=="E":
                                    Elines=open(tmpf+'mtmp.out','r').readlines()
                                    for line in Elines:
                                        if "Total Potential Energy" in line:
                                            if group=='solv':
                                                En_s=float(line.split()[4])
                                                num_s=float(numatoms-numpolyatms)
                                                En.append([En_s,num_s])
                                            elif group=='poly':
                                                En_p=float(line.split()[4])
                                                num_p=float(numpolyatms)
                                                En.append([En_p,num_p])
                                            else:
                                                En_t=float(line.split()[4])
                                                num_t=float(numatoms)
                                                En.append([En_t,num_t])
                                        elif "Van der Waals" in line:
                                            if group=='solv':
                                                Vdw_s=float(line.split()[3])
                                                num_s=float(line.split()[4])
                                                Vdw.append([Vdw_s,num_s])
                                            elif group=='poly':
                                                Vdw_p=float(line.split()[3])
                                                num_p=float(line.split()[4])
                                                Vdw.append([Vdw_p,num_p])
                                            else:
                                                Vdw_t=float(line.split()[3])
                                                num_t=float(line.split()[4])
                                                Vdw.append([Vdw_t,num_t])
                                        elif "Charge-Charge" in line:
                                            if group=='solv':
                                                El_s=float(line.split()[1])
                                                num_s=float(line.split()[2])
                                                El.append([El_s,num_s])
                                            elif group=='poly':
                                                El_p=float(line.split()[1])
                                                num_p=float(line.split()[2])
                                                El.append([El_p,num_p])
                                            else:
                                                El_t=float(line.split()[1])
                                                num_t=float(line.split()[2])
                                                El.append([El_t,num_t])
                                        elif "Intermolecular Energy" in line:
                                            if group=='solv':
                                                Ei_s=float(line.split()[3])
                                                num_s=float(numatoms-numpolyatms)
                                                Ei.append([Ei_s,num_s])
                                            elif group=='poly':
                                                Ei_p=float(line.split()[3])
                                                num_p=float(numpolyatms)
                                                Ei.append([Ei_p,num_p])
                                            else:
                                                Ei_t=float(line.split()[3])
                                                num_t=float(numatoms)
                                                Ei.append([Ei_t,num_t])
                        j=0
                        newfile=1
                        os.system('rm '+tmpf+'mtmp.out')
                        os.system('rm garb.out')
                        os.system('rm '+path+'/'+tmpf+'*')
                        if xyzf.tell()>size:
                            print "End of Archive"
                elif newfile==1:
                    xyzlines+=line
                    tmpf='temp_'+str(num)+'_'+str(runid)
                    #os.system('rm '+path+'/'+tmpf+'*')
                    f=open(path+'/'+tmpf+'.xyz','w')
                    newfile=0
                elif j<=numatom:
                    xyzlines+=line
                j += 1
    else:
        Rg=[]
        En=[]
        Vdw=[]
        El=[]
        Ei=[]
        xyzf=open(path+'/'+file, 'r')
        numatoms=int(xyzf.readline().split()[0])
        numpolyatms=int(open(initpath+'/'+polyfile).readline().split()[0])
        keylines=open(initpath+'/'+prmfile,'r').read()
        if shell==True:
            atom_array,tink_type,atom_type,connect,slines,natoms,cent_size=read_xyz(path,file,polyfile,sol_size)
            shell_atoms,all_atoms,cent_molecule=solv_shell(atom_array,connect,natoms,cent_size,sol_size,cutoff)
            clean_solv_shell,clean_solv_all=deloop_solv(sol_size,shell_atoms,all_atoms,atom_array,cent_molecule,solv_cut)
            write_xyz(path,file,atom_array,atom_type,tink_type,connect,slines,cent_molecule,clean_solv_shell,shell_type)
            infile=file.split('.')[0]+'_shell'
            s=open(path+'/'+infile+'.xyz','r')
            numatoms=int(s.readline().split()[0])
            s.close()
            os.system('xyzsybyl '+path+'/'+infile+'.xyz > garb.out')
            if num==1 or num==3 or num==6 or num==7 or num==8 or num==14 or num==15:
                lf=65
            else:
                lf=100
            CreatePolymer(path,infile,infile,num_mono=1,num_poly=1,lf=lf)
        else:
            infile=file
        for group in ['t','solv','poly']:
            if group=='solv':
                open(initpath+'/'+prmfile+'.tmp','w').write(keylines+"\ninactive -1 "+str(numpolyatms)+"\ndebug")
            elif group=='poly':
                open(initpath+'/'+prmfile+'.tmp','w').write(keylines+"\ninactive -"+str(numpolyatms+1)+" "+str(numatoms)+"\ndebug")
            else:
                open(initpath+'/'+prmfile+'.tmp','w').write(keylines+"\ndebug")
            for item in astr.split():
                print 'analyze -k '+initpath+'/'+prmfile+'.tmp'+' '+path+'/'+infile+' '+item
                os.system('analyze -k '+initpath+'/'+prmfile+'.tmp'+' '+path+'/'+infile+' '+item+'>'+file+group+'mtmp.out')
                if item=="M":
                    Rglines=open(file+group+'mtmp.out','r').readlines()
                    for line in Rglines:
                        if "Radius of Gyration" in line:
                            Rg.append(float(line.split()[4]))
                elif item=="E":
                    Elines=open(file+group+'mtmp.out','r').readlines()
                    for line in Elines:
                        if "Total Potential Energy" in line:
                            if group=='solv':
                                En_s=float(line.split()[4])
                                num_s=float(numatoms-numpolyatms)
                                En.append([En_s,num_s])
                            elif group=='poly':
                                En_p=float(line.split()[4])
                                num_p=float(numpolyatms)
                                En.append([En_p,num_p])
                            else:
                                En_t=float(line.split()[4])
                                num_t=float(numatoms)
                                En.append([En_t,num_t])
                        elif "Van der Waals" in line:
                            if group=='solv':
                                Vdw_s=float(line.split()[3])
                                num_s=float(line.split()[4])
                                Vdw.append([Vdw_s,num_s])
                            elif group=='poly':
                                Vdw_p=float(line.split()[3])
                                num_p=float(line.split()[4])
                                Vdw.append([Vdw_p,num_p])
                            else:
                                Vdw_t=float(line.split()[3])
                                num_t=float(line.split()[4])
                                Vdw.append([Vdw_t,num_t])
                        elif "Charge-Charge" in line:
                            if group=='solv':
                                El_s=float(line.split()[1])
                                num_s=float(line.split()[2])
                                El.append([El_s,num_s])
                            elif group=='poly':
                                El_p=float(line.split()[1])
                                num_p=float(line.split()[2])
                                El.append([El_p,num_p])
                            else:
                                El_t=float(line.split()[1])
                                num_t=float(line.split()[2])
                                El.append([El_t,num_t])
                        elif "Intermolecular Energy" in line:
                            if group=='solv':
                                Ei_s=float(line.split()[3])
                                num_s=float(numatoms-numpolyatms)
                                Ei.append([Ei_s,num_s])
                            elif group=='poly':
                                Ei_p=float(line.split()[3])
                                num_p=float(numpolyatms)
                                Ei.append([Ei_p,num_p])
                            else:
                                Ei_t=float(line.split()[3])
                                num_t=float(numatoms)
                                Ei.append([Ei_t,num_t])
        #os.system('rm '+file+group+'mtmp.out')
    return Rg,En,El,Vdw,Ei
