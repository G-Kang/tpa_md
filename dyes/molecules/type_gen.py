import os,sys

file=sys.argv[1]
pdb=open(file,'r').readlines()
istr=''
j=0
for i,l in enumerate(pdb):
    if 'CONECT' in l:
        break
    elif 'LIG' in l:
        j+=1
        istr+=str(j)+'\n'
f=open(file.split('.')[0]+'.type','w')
f.write(istr)
f.close()
