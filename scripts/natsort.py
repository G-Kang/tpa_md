def try_int(s):
    try: return int(s)
    except: return s

def natsort_key(s):
    import re
    return map(try_int,re.findall(r'(\d+|\D+)',s))

def natcmp(a,b):
    return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a,b):
    return natcmp(a.lower(), b.lower())

def natsort(seq,cmp=natcmp, reverse=False):
    if reverse:
        seq.sort(cmp, reverse=True)
    else:
        seq.sort(cmp)

def natsorted(seq, cmp=natcmp, reverse=False):
    import copy
    temp=copy.copy(seq)
    natsort(temp,cmp,reverse)
    return temp
