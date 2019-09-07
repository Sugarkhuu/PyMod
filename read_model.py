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

blockset = ['var', 'varexo', 'varobs','parameters', 'paramval','model', 'model_obs']


find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem] # find index of set of elements from a list
ind = find(content, blockset)
ind.sort(reverse=False)

mydict = {}
for i in range(len(ind)):
    if i == len(ind)-1:
        mydict[content[ind[i][0]]] = content[ind[i][0]+1:]
    else:
        mydict[content[ind[i][0]]] = content[ind[i][0]+1:ind[i+1][0]]

#eqs = ' '.join(mydict['model'])
eqs = ' '.join(mydict['model']+ mydict['model_obs'])

#eqs = ' '.join(eqs)

nobseq = len(mydict['model_obs'])


import re

eqs = re.split(';', eqs)
eqs = list(filter(lambda a: a != '', eqs))      # Remove blank lines
eqs = [x.strip() for x in eqs]                  # Remove spaces


# Replace LHS=RHS with LHS-(RHS)
eqs_keep = eqs

name_list = ['parameters', 'var', 'varexo', 'varobs']
names = []
nametypes = []
for i in range(len(name_list)):
    names = names + mydict[name_list[i]]
    nametypes = nametypes + [i]*len(mydict[name_list[i]])

list_eq = ['model', 'model_obs']
eq_all = []
eqtntypes = []
for i in range(len(list_eq)):
    eq_all = eq_all + mydict[list_eq[i]]
    eqtntypes = eqtntypes + [i]*len(mydict[list_eq[i]])

neq = len(eqs)
nnormeq = neq

patt = '^([^=]+)(=)?(.+)?$'
repl = r"\1-(\3)"

for i in range(neq):
    eqs[i] = re.sub(patt, repl, eqs[i])

a = eqs[0]


eqsrepl = eqs.copy()

# replace var and params with unicode characters, useful for getting lead and lag
pattchars = mydict['parameters'] + mydict['var'] + mydict['varexo'] +  mydict['varobs']
replchars = []
for i in range(len(pattchars)):
    replchars.append(chr(20480 + i))

for i in range(len(eqsrepl)):
    for j in range(len(pattchars)):
        eqsrepl[i] = re.sub(r'(?<!\w)' +pattchars[j]+r'(?!\w)', replchars[j], eqsrepl[i])
 
mystr = '(([' + str(chr(20480)) + '-' + str(chr(20500)) + '])' + '(\([\+\-]\d+\))?)'
#mystr = '([' + str(chr(20480)) + '-' + str(chr(20500)) + '])' + '(?:\([\+\-]\d+\))?'

re.findall(mystr, eqsrepl[0])

for m in re.finditer(mystr, eqsrepl[0]):
    print(m.group(2))
    print(m.group(3))

mint = 0
maxt = 0
key = []
keys = []

def my_replace(match):
    global key, mint, maxt
    if not not match.group(3):
        s = match.group(3)
        t = s[s.find("(")+1:s.find(")")]
        ind = 't' + t
        key.append(replchars.index(match.group(2))+int(t)*1j)
        mint = min(int(t), mint)
        maxt = max(int(t), maxt)
    else:
        ind = 't'
        key.append(replchars.index(match.group(2))+0j)
    
    return 'x[:,' + str(replchars.index(match.group(2)))  + ',' + ind + ']'

eqeq = eqsrepl
nvarpar = len(pattchars)
eqbaq = []

for i in range(neq):
    key = []
    eqbaq.append(re.sub(mystr, my_replace, eqeq[i]))
    keys.append(key)
'The quick @red2 fox jumps over the @lame4 brown dog.'

nt = maxt - mint + 1
t = maxt - mint + 1


occur = np.zeros((neq, nvarpar, nt))

for eq in range(neq):
    nocc = len(keys[eq])
    for j in range(nocc):
#        mar= np.ravel_multi_index([eq, int(keys[eq][j].real),int(keys[eq][j].imag+1-mint)-1], occur.shape)
        occur[eq, int(keys[eq][j].real),int(keys[eq][j].imag+1-mint)-1] = True
#        occur(mar) = True


fn = lambda x: x[:,1]**2 + 3*x[:,2] - 4

#yaya = eqbaq[3]
#yaya = yaya.replace('t','1')
#yaya = yaya.replace('=', '+')

#t=1
e=0
fms = []
#fm = lambda x, t: eval(eqbaq[0])
for i in range(neq):
#    fm = lambda x, t, e: eval(eqbaq[0])
    fms.append(lambda x, t, e: eval(eqbaq[i]))

fm = lambda i, x, t, e: eval(eqbaq[i])

occurS = occur.copy()
occurS = occurS.reshape(neq,nvarpar*nt, order='F')

from numpy import array
from scipy.sparse import csr_matrix
S = csr_matrix(occurS)
#print(S)

y = np.zeros([3,6])
a = {'b': fn}
a['b'](y)

npar = len(mydict['parameters'])
nvar = len(mydict['var'])
occur_all = occur.copy()
#occur = occur[:,-(nvar+nobs):,:]
occurT = occur[:-nobseq,npar:npar+nvar+nshock+nobs,:]
occurM = occur[-nobseq:,npar:npar+nvar+nshock+nobs,:]
occurTM = occur[:,npar:npar+nvar+nshock+nobs,:]

shift=np.zeros((2,nvar))


for i in range(nvar):
    shift[0,i] = min(min(np.where(occurT[:,i,:]==1)[1]) - nt + 2, shift[0,i])
    shift[1,i] = max(max(np.where(occurT[:,i,:]==1)[1]) - nt + 2, shift[1,i])
    if len(np.where(occurM[:,i,:]==1)[1]) > 0:
        print(i)
        shift[0,i] = min(min(np.where(occurM[:,i,:]==1)[1]) - nt + 1, shift[0,i])
        print(min(np.where(occurM[:,i,:]==1)[1]) - nt + 1)
    if shift[0,i] == shift[1,i]:
       shift[1,i] = 1


Min = int(min(shift[0]))
Max = int(max(shift[1]))

m = {}
m['systemid'] = []

for i in reversed(range(Min, Max+1)):
    print(i)
    aux = np.where((i>=shift[0]) & (i<shift[1]))
    auxb = list(map(lambda x: x + npar + i*1j, list(aux)))
    m['systemid'] = m['systemid'] + list(auxb[:][0][:])

nx = len(m['systemid'])
nu = sum(np.array(m['systemid']).imag >= 0);
npa = nx - nu;


nobs = len(mydict['varobs'])
mshocks = mydict['varexo']
nshock = len((mshocks))
#n = 10
#t=1
#t = m.tzero;
#n = sum(nname(1:3));


m['metaderiv'] = {}
m['metaderiv']['u'] = []
m['metaderiv']['uplus'] = []
m['metaderiv']['p'] = []
m['metaderiv']['e'] = []
m['metaderiv']['pplus'] = []
m['metasystem'] = {}
m['metasystem']['u'] = []
m['metasystem']['e'] = []
m['metasystem']['uplus'] = []
m['metasystem']['p'] = []
m['metasystem']['pplus'] = []
m['systemident'] = {}
m['systemident']['x'] = []
m['systemident']['xplus'] = []


m['metaderiv']['e'] = list(map(lambda x: x + 1 + (t-1)*(nobs+nvar+nshock) + nvar+nshock, list(range(nobs))))
m['metasystem']['e'] = list(range(nobs))

m['metadelete'] = [False]*nu
for i in range(nu):
   if any(m['systemid'][i]-1*1j == m['systemid'][nu+1:]):
       m['metadelete'] = True


for i in range(nu):
   id = m['systemid'][i]
   if id.imag == shift[0][int(id.real)-npar]:
      m['metaderiv']['u'].append(int((id.imag+t-1)*(nobs + nvar+nshock) + id.real))
      m['metasystem']['u'].append(i)
   m['metaderiv']['uplus'].append(int((id.imag+t+1-1)*(nobs + nvar+nshock) + id.real))
   m['metasystem']['uplus'].append(i)

for i in range(npa):
   id = m['systemid'][nu+i]
   if id.imag == shift[0][int(id.real)-npar]:
      m['metaderiv']['p'].append(int((id.imag+t-1)*(nobs + nvar+nshock) + id.real))
      m['metasystem']['p'].append(i)
   m['metaderiv']['pplus'].append(int((id.imag+t+1-1)*(nobs + nvar+nshock) + id.real))
   m['metasystem']['pplus'].append(i)


m['metaderiv']['y'] = list(map(lambda x: x + 1+ (t+ 1 -1)*(nobs + nvar+nshock), range(len(mshocks))))
m['metasystem']['y'] = []
m['metasystem']['y'] = list(range(len(mshocks)))

for i in range(nu+npa):
   id = m['systemid'][i]
   if id.imag != shift[0][int(id.real)-npar]:
      aux = np.zeros(nu+npa)
      aux[m['systemid'] == id-1*1j] = 1
      m['systemident']['xplus'].append(aux)
      aux = np.zeros(nu+npa)
      aux[i] = -1
      m['systemident']['x'].append(aux)


parvals = []
for i in range(len(mydict['paramval'])):
    parvals.append(float(re.findall('=\s*(\d*.\d*)', mydict['paramval'][i])[0]))

deriv = {}
deriv['f'] = np.zeros((neq,(nobs + nvar+nshock)*(nt+2)))

    
for ieq in range(neq):
    print(ieq)
#Soccurn = [19+tshift, 20+tshift, 21+tshift, 19+tshift]
#Soccurt = [0, 0, 0, 1]
    occurA = np.where(occurTM[ieq,:,:]==1)
    Soccurn = occurA[0]+npar
    Soccurt = occurA[1]
    
    #init = np.zeros((len(Soccurn),nvarpar, nt))
    init = np.zeros(nvarpar)
    init[:npar] = parvals
    init = np.broadcast_to(init[None,...,None],(len(Soccurn),) + init.shape+(nt,))
    h = np.ones(init.shape)

    plus  = init + h
    minus = init - h
    step = plus - minus

    fgrid = np.zeros((3,)+init.shape)
    fgrid[:] = init

    for i in range(len(Soccurn)):
        fgrid[0,i,Soccurn[i],Soccurt[i]] = minus[i,Soccurn[i],Soccurt[i]]
        fgrid[2,i,Soccurn[i],Soccurt[i]] = plus[i, Soccurn[i],Soccurt[i]]

    lower = fm(ieq, fgrid[0], t-maxt-1, 0)
    upper = fm(ieq, fgrid[2], t-maxt-1, 0)
    derivs = np.zeros(upper.shape)
    for i in range(len(Soccurn)):
        derivs[i] = (upper[i] - lower[i])/step[1, Soccurn[i], Soccurt[i]]
    index = (Soccurt-1+2)*(nvarpar-npar) + Soccurn
    deriv['f'][ieq][index] = np.array(derivs)


#######################################


nm = nobseq
aant = nnormeq

#mindex = find(eqselect(1:nm))
#tindex = find(eqselect(nm+1:end))


ny = len(m['metaderiv']['y'])
nx = len(m['systemid'])
ne = len(m['metaderiv']['e'])
nf = int(sum((np.array(m['systemid']).imag >= 0).astype('float')))
nb = nx - nf

m['system0']={}
m['system0']['K']=[]
m['system0']['A']=[]
m['system0']['B']=[]
m['system0']['E']=[]

m['system0']['K'].append(np.zeros([ny,1]))
m['system0']['K'].append(np.zeros([nx,1]))
m['system0']['A'].append(np.zeros([ny,ny]))
m['system0']['A'].append(np.zeros([nx,nx]))
m['system0']['B'].append(np.zeros([ny,nb]))
m['system0']['B'].append(np.zeros([nx,nx]))
m['system0']['E'].append(np.zeros([ny,ne]))
m['system0']['E'].append(np.zeros([nx,ne]))

system = m['system0']

# A1 y + B1 xb+ + E1 e + K1 = 0

system['K'][0, mindex] = deriv.c[mindex]
system['K'][0,tindex] = deriv.c(nm+tindex)



system.A{1}(mindex,m.metasystem.y) = deriv.f(mindex,m.metaderiv.y);
%system.B{1}(mindex,m.metasystem.u) = deriv.f(mindex,m.metaderiv.u);
system.B{1}(mindex,m.metasystem.pplus) = deriv.f(mindex,m.metaderiv.pplus);
system.E{1}(mindex,m.metasystem.e) = deriv.f(mindex,m.metaderiv.e);

system['A'][0][-nobseq:,m['metasystem']['y']] = deriv['f'][-nobseq:,m['metaderiv']['y']]
#system['B']{1}(mindex,m.metasystem.u) = deriv.f(mindex,m.metaderiv.u);
system['B'][0][-nobseq:,m['metasystem']['pplus']] = deriv['f'][-nobseq:,m['metaderiv']['pplus']]
system['E'][0][-nobseq:,m['metasystem']['e']] = deriv['f'][-nobseq:,m['metaderiv']['e']]

## A2 [xf+;xb+] + B2 [xf;xb] + E2 e + K2 = 0
#
system['A'][1][:neq-nobseq,m['metasystem']['uplus']] = deriv['f'][:neq-nobseq,m['metaderiv']['uplus']]
system['A'][1][:neq-nobseq,list(nf+np.array(m['metasystem']['pplus']))] = deriv['f'][:neq-nobseq,m['metaderiv']['pplus']]
system['B'][1][:neq-nobseq,m['metasystem']['u']] = deriv['f'][:neq-nobseq,m['metaderiv']['u']]
system['B'][1][:neq-nobseq,list(nf+np.array(m['metasystem']['p']))] = deriv['f'][:neq-nobseq,m['metaderiv']['p']]

system['E'][1][:neq-nobseq,m['metasystem']['e']] = deriv['f'][:neq-nobseq,m['metaderiv']['e']]
#
system['A'][1][neq-nobseq:,:] = m['systemident']['xplus']
system['B'][1][neq-nobseq:,:] = m['systemident']['x']

#Derivatives are almost done (shocks). Now turn them into the matrices.
# Add necessary leads and lags
# numbers on the left and right hand side