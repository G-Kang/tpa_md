from hoomd_script import *
from hoomd_plugins.aniso_potentials import randomize_orientation
import os,sys
import numpy as np
out_dir = sys.argv[1]
overwrite = sys.argv[2]
for (path,dirs,files) in os.walk(out_dir):

    if overwrite=='yes':
        filesu = [os.path.splitext(f)[0] for f in os.listdir(path) if ".xml" in f]
    else:
        filesx = [os.path.splitext(f)[0] for f in os.listdir(path) if ".xml" in f]
        filesp = [os.path.splitext(f)[0] for f in os.listdir(path) if ".pos" in f]
        filesu = np.setdiff1d(filesx,filesp,assume_unique=True)
        print "skipping",filesp,'..'
    for file in filesu:
        print file+'.xml'
        randomize_orientation.xml_to_pos(path+'/'+file+'.xml',"AuShell")

