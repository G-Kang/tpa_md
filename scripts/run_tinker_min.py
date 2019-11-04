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
                        step=0
                        dataf = (data % 
                                {'initdir':path,'poly':poly,
                                'temp':temp,'pres':pres,'size':size,'id':id,'ls':step})
                        open('run_tinker.pbs','w').write(dataf)
                        os.system("msub run_tinker.pbs")	
                        print poly,id


if __name__ == '__main__':
	main()
