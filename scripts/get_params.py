import h5py
from paramInterface import GetMissingParams
from qchemsetup import StoreMissingParams
from qchemsetup import makeQchem
from shutil import copyfile
from os.path import isdir,isfile,join
from os import mkdir,listdir

hFile = 'testdbwithtop.h5'
initdir='molecules'
#xyzFile = 'squaraine-znap.xyz'
xyzlist=[f for f in listdir(initdir) if isfile(join(initdir,f)) and if 'xyz' in f]
#baseName = "squaraine-znap"
for xyzFile in xyzlist:
    baseName=xyzFile.split('.')[0]
    if isdir(baseName) != True:
        mkdir(baseName)
    missingParams, newAtoms, missingParams2, \
        fpDict, dihedralGroups, \
        fullGraph, frSets  = GetMissingParams(initdir+'/'+ xyzFile, hFile)

    StoreMissingParams(newAtoms, fpDict, missingParams2,
                       dihedralGroups, baseName)

    start = '0' # Interval start stop for potential scan
    end = '180'
    interval = '10'

    # You'll have to loop through a collection of $param to make all the files to get all the parameters
    # This will prep the basic DFT calculations
    for key, val in dihedralGroups.items():
        param = val[0] #dihedral Groups stores every duplicate of a particular 
                        #dihedral fingerprint it finds in a list, so you only need to parameterize one of them
        makeQchem(baseName, param, fullGraph, frSets, start, 
                    end, interval, charge=0, mult=1, template='qTemplate.in', makeXYZ=True)

        if isdir(baseName):
            titleString = ''
            for p in param:
                titleString += '_' + str(p)
            qchemtorpath=baseName+'/torsion'+titleString
            if isdir(qchemtorpath) != True:
                mkdir(qchemtorpath)
            if len(param) == 2:
                paramType = 'stre'
            elif len(param) == 3:
                paramType = 'bend'
            elif len(param) == 4:
                paramType = 'tors'
            newXYZFileName = baseName + '/' + baseName + '_' + paramType + titleString + '.xyz'
            qFileName = baseName + '/' + baseName + '_' + paramType + titleString + '.in'
            copyfile(newXYZFileName,qchemtorpath+'/input.xyz')
            copyfile(qFileName,qchemtorpath+'/input.in')
            copyfile('qchem.sh',qchemtorpath+'/qchem.sh')
