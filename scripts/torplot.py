import os,sys,math
from natsort import natsorted
path=sys.argv[1]
tornum=sys.argv[2]
files=[f for f in os.listdir(path) if tornum+'tor' in f and '.out' in f and not 'full' in f]
files=natsorted(files)
print files

toren=[]
poten=[]
potred=[]
angle=[math.pi*a/180. for a in xrange(0,190,10)]
for file in files:
    lines=open(file,'r').readlines()
    for line in lines:
        if 'Torsional' in line:
            toren.append(float(line.split()[2]))
        elif 'Total' in line:
            poten.append(float(line.split()[4]))
ft=open(path+'/torplot'+tornum+'.txt','w')
fe=open(path+'/ntpotplot'+tornum+'.txt','w')
fte=open(path+'/tpotplot'+tornum+'.txt','w')
for t,e in zip(toren,poten):
    potred.append(e-t)
red_min=min(potred)
pot_min=min(poten)
tor_min=min(toren)
for t,e,te,a in zip(toren,potred,poten,angle):
    ft.write("%.2f\t%.2f\n" % (a,t-tor_min))
    fe.write("%.2f\t%.2f\n" % (a,e-red_min))
    fte.write("%.2f\t%.2f\n" % (a,te-pot_min))
ft.close()
fe.close()
fte.close()
