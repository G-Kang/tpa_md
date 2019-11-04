import os,sys,math,shutil
from natsort import natsorted
#path=sys.argv[1]
#polynum=sys.argv[2]
#solvtyp=sys.argv[3]
def lastStep(path,polynum,solvtyp,copyf=False,copydir=''):
    files=[f for f in os.listdir(path) if polynum in f and solvtyp in f and f.split('.')[1].isdigit()==True]
    files=natsorted(files)
    if len(files)<1:
        return 0
    else:
        print files[-1]
        if copyf=='True' and int(files[-1].split('.')[1])>=700:
            print copydir
            shutil.copyfile(path+'/'+files[-1],copydir+'/'+files[-1].split('.')[0]+'.xyz')
        return int(files[-1].split('.')[1])

if __name__=='__main__':
    path=sys.argv[1]
    polynum=sys.argv[2]
    solvtyp=sys.argv[3]
    copyf=sys.argv[4]
    copydir=sys.argv[5]
    lS=lastStep(path,polynum,solvtyp,copyf=copyf,copydir=copydir)
