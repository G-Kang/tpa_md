import os, sys
import mmap
from filesort import lastStep

def drange(start, stop, step):
    r = start
    while r <= stop:
        yield r
        r += step

def tail(f, offset,window=1400):
    if offset==2:
        BUFSIZ = 1024
        f.seek(0,os.SEEK_END)
    else:
        BUFSIZ = 4048*20
        f.seek(offset,0)
    size = window
    block = -1
    data = ''
    exit = False
    while not exit:
        step = (block * BUFSIZ)
        if abs(step) >= f.tell():
            f.seek(0)
            newdata =f.read(BUFSIZ - (abs(step) - offset))
            exit = True
        else:
            if offset==2:
                f.seek(step,os.SEEK_END)
            else:
                f.seek(step,1)
            newdata = f.read(BUFSIZ)
        data = newdata + data
        if data.count('\n') >= window:
            break
        else:
            block -= 1
    return data.splitlines()[-window:]
def main():
    data = open('run_tinker.pbs-python', 'r').read()
    path = "init/solvent"
    restart=False
    lastTime=0
    i = 0
#    for xml_file in [f for f in os.listdir(path) if 'kayak' in f]:
#        xmlfilecat = xml_file.partition(".")[0]
    numstarts=1
    for poly in drange(1, 15, 1):
        for temp in drange(298, 298, 1):
            for pres in drange(0.,0.,.4):
                for size in [4]:
                    for id in drange(0,0,1):
                        if restart==True:
                            with open('D'+str(poly)+'_'+str(id)+'_'+str(size)+'dyncbsol.log', 'r') as logf:
                                numrun = 0
                                for line in logf:
                                    if 'Random' in line and numrun==numstarts:
                                        numrun += 1
                                        startsize = logf.tell()
                                        block = tail(logf,startsize)
                                    elif 'Random' in line:
                                        numrun += 1
                            times = []
                            for l in block:
                                if 'Instantaneous' in l:
                                    times.append(int(l.split()[6]))
                            if len(times)==0:
                                times.append(0)
                            lastTime = times[-1]
                            print lastTime
                            logf.close()
                            block = tail(open('D'+str(poly)+'_'+str(id)+'_'+str(size)+'dyncbsol.log'),2)
                            for l in block:
                                if 'Instantaneous' in l:
                                    times.append(int(l.split()[6]))
                            if len(times)==0:
                                times.append(0)
                            lastTime += times[-1]
                            print times[-1]
                        else:
                            #lastTime = lastStep('data/dynamic','D'+str(poly)+'_0','cbsol',copyf=True)
                            lastTime=0
                            print lastTime
                        step = 1000000 - lastTime
                        if step > 0:
                            print 'running ',step
                            dataf = (data % 
                                    {'initdir':path,'poly':poly,
                                    'temp':temp,'pres':pres,'size':size,'id':id,'ls':step})
                            open('run_tinker.pbs','w').write(dataf)
                            os.system("msub run_tinker.pbs")	
                        #elif lastTime > 500 and lastTime < 700:
                        #    print 'close enough to finish',500000-lastTime*1000
                        #    dataf = (data % 
                        #            {'initdir':path,'poly':poly,
                        #            'temp':temp,'pres':pres,'size':size,'id':id,'ls':step})
                        #    open('run_tinker.pbs','w').write(dataf)
                        #    #os.system("msub run_tinker.pbs")	
                        else:
                            print "Past 500 ps by ",lastTime/1000
                        print poly,id


if __name__ == '__main__':
	main()
