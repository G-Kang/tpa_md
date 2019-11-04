import os,sys
import collections
# param file
#prmf=sys.argv[1]
# num torsions of type 38
#aval=int(sys.argv[2])
# num torsions of type 74
#bval=int(sys.argv[3])
# A-type first
#cllist=list(sys.argv[4])

def prmread(prmf):
    file=open(prmf,'r')
    fread=file.read()
    file=open(prmf,'r')
    lines=file.readlines()
    #if "Thiophene" in fread:
    #    thiophene='True'
    #else:
    #    thiophene='False'
    prmstr=''
    natoms=0
    cltypes={}
    #Ntyp='47'
    #Styp='16'
    for k,l in enumerate(lines):
        if 'atom' in l:
            natoms+=1
            type=l.split()[1]
            cltyp=l.split()[2]
            #if cltyp in cldict:
            #    cltypes[cltyp]=100+int(cltyp)
            #else:
            #    cltypes[cltyp]=int(cltyp)
            itype=int(type)
            itype=1000+k
            nline=''
            for j,i in enumerate(l.split()):
                if j==1:
                    nline+=str(itype)+' '
                else:
                    nline+=i+' '
            #nline=nline.replace(cltyp,str(cltypes[cltyp]))
            prmstr+=nline+'\n'
            if 'vdw' in lines[k+1]:
                for numc in xrange(0,natoms):
                    prmstr+='charge '+str(itype-natoms+numc+1)+' c\n'
        elif 'vdw' in l:
            type=l.split()[1]
            itype=int(type)
            itype=1000+k-natoms
            nline=''
            for j,i in enumerate(l.split()):
                if j==1:
                    nline+=str(itype)+' '
                else:
                    nline+=i+' '
            prmstr+=nline+'\n'
        else:
            natoms=0
    return prmstr

#if Ntyp not in cltypes:
#    cltypes['47']=1
def getopls(prmstr,cltypes,cldict):
    print( cltypes)
    lines=open('oplsaa.prm','r').readlines()
    bl=''
    al=''
    dl=''
    ndih={'38':aval,'74':bval}
    if ndih['38']==0:
        typs=['74']
    elif ndih['74']==0:
        typs=['38']
    else:
        if abool=='yes':
            typs=['38','74']
        else:
            typs=['74','38']
    mult=[1,2,3,4]
    #mblist={1:[2],2:[2]}
    #malist={1:[3],2:[-1,1],3:[-2,-1,1]}
    #mdlist={1:[4],2:[-1,1],3:[-2,-1,1],4:[-3,-2,-1,1]}
    def findnth(source, target, n): 
        num = 0 
        start = -1 
        while num < n: 
            start = source.find(target, start+1) 
            if start == -1: return -1
            num += 1 
        return start 

    def replacenth(source, old, new, n): 
        p = findnth(source, old, n) 
        if n == -1: return source 
        return source[:p] + new + source[p+len(old):] 

    for m in mult:
        for l in lines:
            if 'bond' in l and m<3:
                bcl = l.split()
                for k,typ in enumerate(typs):
                    counter=collections.Counter(bcl[1:3])
                    if counter[typ]==m:
                        for j in xrange(0,ndih[typ]):
                            ncl=101+j+k
                            for i,at in enumerate(bcl[1:3]):
                                if at==typ and i==0 and all((k in cltypes) for k,v in counter.iteritems()):
                                    #if Ntyp in cltypes:
                                        #l=l.replace(Ntyp,Styp)
                                    bl+=l.replace(typ,str(ncl),1)
                                    #print l.replace(typ,str(ncl),1)
                                    #print l
                                elif at==typ and i==1 and bcl[1]!=typ and all((k in cltypes) for k,v in counter.iteritems()):
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    bl+=replacenth(l,typ,str(ncl),1)
                                    #print replacenth(l,typ,str(ncl),1)
                                    #print l
            elif 'angle' in l and m<4:
                acl = l.split()
                for k,typ in enumerate(typs):
                    counter=collections.Counter(acl[1:4])
                    if counter[typ]==m:
                        for j in xrange(0,ndih[typ]):
                            ncl=101+j+k
                            for i,at in enumerate(acl[1:4]):
                                if at==typ and i==0 and all((k in cltypes) for k,v in counter.iteritems()):
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    al+=l.replace(typ,str(ncl),1)
                                    #print l.replace(typ,str(ncl),1)
                                    #print l
                                elif at==typ and i==1 and all((k in cltypes) for k,v in counter.iteritems()):
                                    if acl[1]==typ:
                                        seq=2
                                    else:
                                        seq=1
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    al+=replacenth(l,typ,str(ncl),seq)
                                    #print replacenth(l,typ,str(ncl),seq)
                                    #print l
                                elif at==typ and i==2 and all((k in cltypes) for k,v in counter.iteritems()):
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    al+=replacenth(l,typ,str(ncl),m)
                                    #print replacenth(l,typ,str(ncl),m)
                                    #print l
            elif 'torsion' in l and m<5:
                dcl = l.split()
                for k,typ in enumerate(typs):
                    counter=collections.Counter(dcl[1:5])
                    if counter[typ]==m:
                        for j in xrange(0,ndih[typ]):
                            ncl=101+j+k
                            for i,at in enumerate(dcl[1:5]):
                                if at==typ and i==0 and all((k in cltypes) for k,v in counter.iteritems()):
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    dl+=l.replace(typ,str(ncl),1)
                                    #print l.replace(typ,str(ncl),1)
                                elif at==typ and i==1 and all((k in cltypes) for k,v in counter.iteritems()):
                                    if dcl[1]==typ:
                                        seq=2
                                    else:
                                        seq=1
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    dl+=replacenth(l,typ,str(ncl),seq)
                                    #print replacenth(l,typ,str(ncl),seq)
                                elif at==typ and i==2 and all((k in cltypes) for k,v in counter.iteritems()):
                                    if dcl[1]==typ and dcl[2]==typ:
                                        seq=3
                                    elif dcl[1]==typ or dcl[2]==typ:
                                        seq=2
                                    else:
                                        seq=1
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    dl+=replacenth(l,typ,str(ncl),seq)
                                    #print replacenth(l,typ,str(ncl),seq)
                                elif at==typ and i==3 and all((k in cltypes) for k,v in counter.iteritems()):
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    dl+=replacenth(l,typ,str(ncl),m)
                                    #print replacenth(l,typ,str(ncl),m)
    if (aval+bval)>1:
        for l in lines:
            if 'bond' in l:
                bcl = l.split()
                typ1=typs[0]
                typ2=typs[1]
                if typ1 in bcl[1:3] and typ2 in bcl[1:3]:
                    for j in xrange(0,ndih[typ1]):
                        ncl1=101+j
                        for k in xrange(j,ndih[typ2]+1):
                            if k==0:
                                continue
                            ncl2=ncl1+k
                            #if Ntyp in cltypes:
                            #    l=l.replace(Ntyp,Styp)
                            lh=l.replace(typ1,str(ncl1))
                            bl+=lh.replace(typ2,str(ncl2))
                            #print lh.replace(typ2,str(ncl2))
                            #print l
            elif 'angle' in l:
                acl = l.split()
                typ1=typs[0]
                typ2=typs[1]
                counter=collections.Counter(acl[1:4])
                if typ1 in acl[1:4] and typ2 in acl[1:4] and all((k in cltypes) for k,v in counter.iteritems()):
                    for j in xrange(0,ndih[typ1]):
                        ncl1=101+j
                        for k in xrange(j,ndih[typ2]+1):
                            if k==0:
                                continue
                            ncl2=ncl1+k
                            if counter[typ1]==2:
                                for seq in [1,2]:
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    lh=replacenth(l,typ1,str(ncl1),seq)
                                    al+=lh.replace(typ2,str(ncl2),1)
                                    #print lh.replace(typ2,str(ncl2),1)
                                    #print l
                            elif counter[typ2]==2:
                                for seq in [1,2]:
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    lh=replacenth(l,typ2,str(ncl2),seq)
                                    al+=lh.replace(typ1,str(ncl1),1)
                                    #print lh.replace(typ1,str(ncl1),1)
                                    #print l
                            else:
                                #if Ntyp in cltypes:
                                #    l=l.replace(Ntyp,Styp)
                                lh=l.replace(typ1,str(ncl1),1)
                                al+=lh.replace(typ2,str(ncl2),1)
                                #print lh.replace(typ2,str(ncl2),1)
                                #print l
            elif 'torsion' in l:
                dcl = l.split()
                typ1=typs[0]
                typ2=typs[1]
                counter=collections.Counter(dcl[1:4])
                if typ1 in dcl[1:5] and typ2 in dcl[1:5] and all((k in cltypes) for k,v in counter.iteritems()):
                    for j in xrange(0,ndih[typ1]):
                        ncl1=101+j
                        for k in xrange(j,ndih[typ2]+1):
                            if k==0:
                                continue
                            ncl2=ncl1+k
                            if counter[typ1]==3:
                                for seq in [1,2,3]:
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    lh=replacenth(l,typ1,str(ncl1),seq)
                                    dl+=lh.replace(typ2,str(ncl2),1)
                                    #print lh.replace(typ2,str(ncl2),1)
                                    #print l
                            elif counter[typ2]==3:
                                for seq in [1,2,3]:
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    lh=replacenth(l,typ2,str(ncl2),seq)
                                    dl+=lh.replace(typ1,str(ncl1),1)
                                    #print lh.replace(typ1,str(ncl1),1)
                                    #print l
                            elif counter[typ1]==2:
                                for seq in [1,2]:
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    lh=replacenth(l,typ1,str(ncl1),seq)
                                    dl+=lh.replace(typ2,str(ncl2),1)
                                    #print lh.replace(typ2,str(ncl2),1)
                                    #print l
                            elif counter[typ2]==2:
                                for seq in [1,2]:
                                    #if Ntyp in cltypes:
                                    #    l=l.replace(Ntyp,Styp)
                                    lh=replacenth(l,typ2,str(ncl2),seq)
                                    dl+=lh.replace(typ1,str(ncl1),1)
                                    #print lh.replace(typ1,str(ncl1),1)
                                    #print l
                            else:
                                #if Ntyp in cltypes:
                                #    l=l.replace(Ntyp,Styp)
                                lh=l.replace(typ1,str(ncl1),1)
                                dl+=lh.replace(typ2,str(ncl2),1)
                                #print lh.replace(typ2,str(ncl2),1)
                                #print l
    return prmstr,bl,al,dl
def writeprm(prmf,prmstr,bl='',al='',dl=''):
    prmnew=open(prmf.split('.')[0]+'new.prm','w')
    #if thiophene=='True':
    #    dl+='torsion\t101 16 102 103 0.000 0.0 1   7.250 180.0 2   0.000 0.0 3\n'
    #    dl+='torsion\t102 16 103 104 0.000 0.0 1   7.250 180.0 2   0.000 0.0 3\n'
    #    dl+='torsion\t103 16 104 105 0.000 0.0 1   7.250 180.0 2   0.000 0.0 3\n'
    #    dl+='torsion\t104 16 105 106 0.000 0.0 1   7.250 180.0 2   0.000 0.0 3\n'
    #    dl+='torsion\t103 16 102 101 0.000 0.0 1   7.250 180.0 2   0.000 0.0 3\n'
    #    dl+='torsion\t104 16 103 102 0.000 0.0 1   7.250 180.0 2   0.000 0.0 3\n'
    #    dl+='torsion\t105 16 104 103 0.000 0.0 1   7.250 180.0 2   0.000 0.0 3\n'
    #    dl+='torsion\t106 16 105 104 0.000 0.0 1   7.250 180.0 2   0.000 0.0 3\n'
    prmnew.write(prmstr)
    prmnew.write(bl)
    prmnew.write(al)
    prmnew.write(dl)
    prmnew.close()

def writeoprm(prmf):
    prmo=open(prmf.split('.')[0]+'o.prm','w')
    prmnew=open(prmf,'r').read()
    op=open('oplsaa.prm','r').read()
    prmo.write(op)
    prmo.write(prmnew)
    prmo.close()
