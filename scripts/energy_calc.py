from energetics import arctonrg

path='data/dynamic'
initpath='init/solvent'
poly=1
solv='cfsol'
shell=False
files=['D'+str(i)+'_0_poly4_cent_min_'+solv+'_all2.arc' for i in xrange(poly,poly+15)]
#files=[solv+'_box.arc']
print files
for file in files:
    Rg,En,El,Vdw,Ei=arctonrg(path,initpath,file,'M E',period=49,archive=True,shell=shell)
    print "calc done"
    if shell==True:
        outRg=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_shell_Rg.txt','w')
        outEn=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_shell_En.txt','w')
        outEl=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_shell_El.txt','w')
        outEi=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_shell_Ei.txt','w')
        outVdw=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_shell_Vdw.txt','w')
    else:
        outRg=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_Rg.txt','w')
        outEn=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_En.txt','w')
        outEl=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_El.txt','w')
        outEi=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_Ei.txt','w')
        outVdw=open(path+'/'+'_'.join(file.split('_')[0:2])+'_'+solv+'_Vdw.txt','w')
        #outRg=open(path+'/'+solv+'_Rg.txt','w')
        #outEn=open(path+'/'+solv+'_En.txt','w')
        #outEl=open(path+'/'+solv+'_El.txt','w')
        #outEi=open(path+'/'+solv+'_Ei.txt','w')
        #outVdw=open(path+'/'+solv+'_Vdw.txt','w')
    for itR,itEn,itEl,itV,itEi in zip(Rg,En,El,Vdw,Ei):
        outRg.write(str(itR)+'\n')
        outEn.write(str(itEn[0])+'\t'+str(itEn[1])+'\n')
        outEl.write(str(itEl[0])+'\t'+str(itEl[1])+'\n')
        outEi.write(str(itEi[0])+'\t'+str(itEi[1])+'\n')
        outVdw.write(str(itV[0])+'\t'+str(itV[1])+'\n')
    outRg.close()
    outEn.close()
    outEl.close()
    outEi.close()
    outVdw.close()

