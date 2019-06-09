# Parse AR(1) equations given in string forms 
# into A*x = B*x(-1) + e form just by evaluating coefficients
# No derivatives taken, just evaluating equations
# Finding A and B matrices so far

def remblank(alist):
    alist = list(filter(lambda a: a != '', alist))
    return alist

import pdb
import numpy as np

fname = 'eqx.py'
with open(fname) as f:
    content = f.readlines()

content = [x.strip() for x in content]                  # Remove \n
content = remblank(content)     # Remove blank lines

blockset = ['var', 'varexo', 'parameters', 'paramval','model']


find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem] # find index of set of elements from a list
ind = find(content, blockset)
ind.sort(reverse=False)

mydict = {}
for i in range(len(ind)):
    if i == len(ind)-1:
        mydict[content[ind[i][0]]] = content[ind[i][0]+1:]
    else:
        mydict[content[ind[i][0]]] = content[ind[i][0]+1:ind[i+1][0]]

eqs = ' '.join(mydict['model'])
#eqs = ' '.join(eqs)

import re

eqs = re.split(';', eqs)
eqs = list(filter(lambda a: a != '', eqs))      # Remove blank lines
eqs = [x.strip() for x in eqs]                  # Remove spaces

a = eqs[0]

#Parameterize A and B matrix
A = np.zeros([len(eqs), len(mydict['var'])])
B = np.zeros([len(eqs), len(mydict['var'])])

for i in range(len(eqs)):
    a=eqs[i]
#re.split('((?<!\()[+\-*\/=])', b) # decomposing equation into variables
    c=re.split('((?<!\()[+\-=])', a) # decomposing equation into multiples
    for j in range(len(c)):
        d=re.split('([*\/])', c[j])
        pos, parval =[],[]
        for k in range(len(d)):
            myvar = d[k].strip()
            myvar = re.split('(\([\-+]\d+\))',myvar)
            if myvar == []:
                if d[k].strip() in mydict['var']:
                   pos = mydict['var'].index(d[k].strip())
                if d[k].strip() in mydict['parameters']:
                   parpos = mydict['parameters'].index(d[k].strip())
                   parval = re.findall('(?<==\s)\d+',mydict['paramval'][parpos])
                ind=0
            else: 
                if myvar[0] in mydict['var']:
                   pos = mydict['var'].index(myvar[0])
                if myvar[0] in mydict['parameters']:
                   parpos = mydict['parameters'].index(myvar[0])
                   parval = re.findall('(?<==\s)\d+',mydict['paramval'][parpos])
                if len(myvar)==1:
                    ind=1
                else:
                    ind=0
        if pos != [] and parval != [] and ind == 1:
            A[i, pos] = float(parval[0])
        if pos != [] and parval != [] and ind == 0:
            B[i, pos] = float(parval[0])