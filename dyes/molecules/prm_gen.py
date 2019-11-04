import os,sys
import collections
from prmpredit import prmread,writeprm,writeoprm

cldict1={'16':1,'74':1,'77':1,'50':1,'71':1}
cldict2={}
cldict3={}
cldict4={}
cldict5={}
cldict6={}
cldict7={}
cldict8={}
cldict9={}
cldict10={}
cldict11={}
cldict12={}
cldict13={}
cldict14={}
cldict15={}
cldict=[eval('cldict'+str(i)) for i in range(1,16)]
print( cldict)
molecules=['diphenylethylene','diphenylstilbenen2','stilbenen2']
for i in molecules:
    file=i+'.prm'
    #prmstr=prmread(file)
    #writeprm(file,prmstr)
    writeoprm(file)
    os.system("prmedit "+ file.split('.')[0]+'o 5')
    os.system("cp parameter.prm "+ file.split('.')[0]+'o.prm')
    os.system("rm parameter.prm")
