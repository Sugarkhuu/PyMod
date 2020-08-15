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

nobs = len(mydict['varobs'])
mshocks = mydict['varexo']
nshock = len((mshocks))

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


#nobs = len(mydict['varobs'])
#mshocks = mydict['varexo']
#nshock = len((mshocks))
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
   if any(m['systemid'][i]-1*1j == m['systemid'][nu:]):
       m['metadelete'][i] = True


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

#system['K'][0, mindex] = deriv.c[mindex]
#system['K'][0,tindex] = deriv.c(nm+tindex)


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

# Solve section!!!


import scipy

ny = nobs
nx = len(m['systemid'])
nb = sum(np.array(m['systemid']).imag < 0)
nf = nx - nb
ne = nshock
#fkeep = ~m.metadelete;
#nfkeep = sum(fkeep);
#nalt = size(m.assign,3);

AA = system['A'][1]
BB = system['B'][1]


#Line 379 of _decomp_qz.py of scipy.linalg
#    liwork = None
#    if liwork is None or liwork == -1:
#        result = tgsen(select, AA, BB, Q, Z, liwork=-1)
#        liwork = result[-2][0]
##    pdb.set_trace()
#    result = tgsen(select, AA, BB, Q, Z, lwork=lwork, liwork=liwork)
#    [SS,TT,alpha,beta,QQ,ZZ] = result[0], result[1], result[2], result[4], result[5], result[6]
#    sfunction = _select_function(sort2)
#    select = sfunction(alpha, beta)
#    result = tgsen(select, SS, TT, QQ, ZZ, lwork=lwork, liwork=liwork)

[SS,TT,QQ,ZZ] = scipy.linalg.qz(AA,BB)
#[SS,TT,alpha,beta,QQ,ZZ] = scipy.linalg.ordqz(AA,BB,'rhp')
[SS,TT,alpha,beta,QQ,ZZ] = scipy.linalg.ordqz(AA,BB,sort=sorter1,sort2=sorter2,output='real')
alpha/beta


AA =[[0,0,0,-1.0000,0.4000,0,0,0,0,0,0],
    [1.0000,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0],
         [0,         0,         0,         0,   -1.0000,         0,         0,         0,         0,         0,         0],
         [0,         0,         0,         0,    0.7000,         0,         0,         0,         0,         0,         0],
         [0,         0,         0,         0,         0,   -1.0000,         0,         0,         0,         0,         0],
         [0,         0,         0,         0,         0,         0,   -1.0000,         0,         0,         0,         0],
         [0,         0,         0,    1.0000,         0,         0,         0,         0,         0,         0,         0],
         [0,         0,         0,         0,         0,         0,         0,    1.0000,         0,         0,         0],
         [0,         0,         0,         0,         0,         0,         0,         0,    1.0000,         0,         0],
         [0,         0,         0,         0,         0,         0,         0,         0,         0,    1.0000,         0],
         [0,         0,         0,         0,         0,         0,         0,         0,         0,         0,    1.0000]]

BB =    [[0,    0.2000,         0,         0,         0,         0,    1.0000,         0,         0,         0,    0.8000],
         [0,   -1.0000,         0,         0,         0,         0,         0,         0,         0,         0,         0],
         [0,         0,         0,         0,    0.9000,         0,         0,         0,         0,         0,         0],
         [0,         0,   -1.0000,         0,         0,         0,         0,         0,         0,         0,         0],
         [0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0],
         [0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0],
   [-1.0000,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0],
         [0,         0,         0,   -1.0000,         0,         0,         0,         0,         0,         0,         0],
         [0,         0,         0,         0,         0,         0,         0,   -1.0000,         0,         0,         0],
         [0,         0,         0,         0,         0,         0,         0,         0,   -1.0000,         0,         0],
         [0,         0,         0,         0,         0,         0,         0,         0,         0,   -1.0000,         0]]

#matlab:
# [AA,BB,Q,Z] = qz(A,B)
# Q*A*Z = AA and Q*B*Z = BB
#python:
#AA, BB, Q, Z = linalg.qz(A, B)
#(A,B) = (Q*AA*Z', Q*BB*Z')
#Q'*A*Z = AA and Q'*B*Z = BB


#eigval = - alpha/beta
#abs(eigval)
#tolerance = 1e-10
#stable = abs(eigval) >= 1 + tolerance
#unit = abs(abs(eigval)-1) < tolerance
#clusters = np.zeros(eigval.shape)
#clusters[:] = False
#clusters[unit] = True
#clusters[stable] = True

def sorter(a, b):
    global clusters
    eigval = - alpha/beta
    tolerance = 1e-10
    stable = abs(eigval) >= 1 + tolerance
    unit = abs(abs(eigval)-1) < tolerance
    clusters = np.zeros(eigval.shape)
    clusters[unit] = 2
    clusters[stable] = 1

def sorter1(a,b):
    # this works well
      global clusters
      b[abs(b)== 0] = 0.001
      eigval = - a/b
      tolerance = 1e-5
      stable = abs(eigval) > 1 + tolerance
      # stable = abs(eigval) >= 1 + tolerance
      unit = abs(abs(eigval)-1) < tolerance
      clusters = np.zeros(eigval.shape)
      clusters[stable] = 1
      clusters[unit] = 1
      return clusters

def sorter2(a,b):
    # this works well
      global clusters
      b[abs(b)== 0] = 0.001
      eigval = - a/b
      tolerance = 1e-5
#      stable = abs(eigval) > 1 + tolerance
      # stable = abs(eigval) >= 1 + tolerance
      unit = abs(abs(eigval)-1) < tolerance
      clusters = np.zeros(eigval.shape)
#      clusters[stable] = 1
      clusters[unit] = 1
      return clusters
  
#def sorter(alhpa, beta):
#    return clusters
##    if abs(-a/b) >= 1 + tolerance:
##        c = 2
##    elif abs(-a/b) < 1 + tolerance:
##        c = 1
##    else:
##        c = 0
##    return c
#
#[SS,TT,alpha,beta,QQ,ZZ] = scipy.linalg.ordqz(SS,TT,sort=sorter)

#fkeep = ~m['metadelete']
import operator
fkeep = list(map(operator.not_, m['metadelete']))
nfkeep = sum(fkeep)
fkeep = fkeep + [False]*(ZZ.shape[0] - len(fkeep))

flag = True
C = np.dot(QQ,system['K'][1])
D = np.dot(QQ,system['E'][1])
S11 = SS[:nb,:nb]
S12 = SS[:nb,nb:]
S22 = SS[nb:,nb:]
T11 = TT[:nb,:nb]
T12 = TT[:nb,nb:]
T22 = TT[nb:,nb:]
Z11 = ZZ[fkeep,:nb]
Z12 = ZZ[fkeep,nb:]
Z21 = ZZ[nf:,:nb]
Z22 = ZZ[nf:,nb:]
C1 = C[:nb,0]
C2 = C[nb:,0]
D1 = D[:nb,:]
D2 = D[nb:,:]

U = Z21

G = -np.dot(np.linalg.inv(Z21),Z22)

# =============================================================================
#     % Unstable block.
# 
#   G = -Z21\Z22;
#   if any(isnan(G(:)))
#     flag = false;
#     return
#   end
#   Ru = -T22\D2;
#   if any(isnan(Ru(:)))
#     flag = false;
#     return
#   end
#   if m.linear
#     Ku = -(S22+T22)\C2;
#   else
#     Ku = zeros([nfkeep,1]);
#   end
#   if any(isnan(Ku(:)))
#     flag = false;
#     return
#   end
# 
#   % Transform stable block == transform backward-looking variables:
#   % a(t) = s(t) + G u(t+1).
# 
#   Ta = -S11\T11;
#   if any(isnan(Ta(:)))
#     flag = false;
#     return
#   end
#   Xa0 = S11\(T11*G + T12);
#   if any(isnan(Xa0(:)))
#     flag = false;
#     return
#   end
#   Ra = -Xa0*Ru - S11\D1;
#   if any(isnan(Ra(:)))
#     flag = false;
#     return
#   end
#   Xa1 = G + S11\S12;
#   if any(isnan(Xa1(:)))
#     flag = false;
#     return
#   end
#   if m.linear
#     Ka = -(Xa0 + Xa1)*Ku - S11\C1;
#   else
#     Ka = asstate(:,2) - Ta*asstate(:,1);
#   end
#   if any(isnan(Ka(:)))
#     flag = false;
#     return
#   end
# 
#   % Forward-looking variables.
# 
#   % Duplicit rows (metadelete) already deleted from Z11 and Z12.
#   Tf = Z11;
#   Xf = Z11*G + Z12;
#   Rf = Xf*Ru;
#   if m.linear
#     Kf = Xf*Ku;
#   else
#     Kf = xfsstate(:,2) - Tf*asstate(:,1);
#   end
#   if any(isnan(Kf(:)))
#     flag = false;
#     return
#   end
# 
#    % State-space form:
#    % [xf(t);a(t)] = T a(t-1) + K + R(L) e(t),
#    % U a(t) = xb(t).
#    T = [Tf;Ta];
#    K = [Kf;Ka];
#    R = [Rf;Ra];
# 
#    % Remove negligible entries from U.
#    sing = svd(U);
#    tol = size(U,1)*eps(sing(1))^(2/3);
#    U(abs(U) <= tol) = 0;
#   
#    % y(t) = Z a(t) + D + H e(t)
#    if ny > 0
#       Z = -full(system.A{1}\system.B{1});
#       if any(isnan(Z(:)))
#          flag = false;
#          return,
#       end
#       H = -full(system.A{1})\full(system.E{1}); % eye(ny);
#       if any(isnan(H(:)))
#          flag = false;
#          return
#       end
#       if m.linear
#          D = full(-system.A{1}\system.K{1});
#       else
#          D = ysstate - Z*xbsstate(:,2); % mirek told ondra how to fix this
#       end
#       if any(isnan(D(:)))
#          flag = false;
#          return
#       end
#       % Remove negligible entries from Z.
#       sing = svd(Z);
#       tol = size(Z,2)*eps(sing(1))^(2/3);
#       Z(abs(Z) <= tol) = 0;
#       Z = Z*U;
#    else
#       Z = zeros([0,nb]);
#       H = zeros([0,ne]);
#       D = zeros([0,1]);
#    end
# 
#   % Necessary initial conditions in xb vector.
#   m.icondix(1,:,ialt) = any(abs(T/U) > tolerance,1);
# 
#   % Forward expansion.
#   % a(t) <- -Xa J^(k-1) Ru e(t+k)
#   % xf(t) <- Xf J^k Ru e(t+k)
#   J = -T22\S22;
#   Xa = Xa1 + Xa0*J;
#   Jk = eye(size(J)); % highest computed power of J: e(t+k) requires J^k
# 
#   m.expand{1}(:,:,ialt) = Xa;
#   m.expand{2}(:,:,ialt) = Xf;
#   m.expand{3}(:,:,ialt) = Ru;
#   m.expand{4}(:,:,ialt) = J;
#   m.expand{5}(:,:,ialt) = Jk;
# 
#   m.solution{1}(:,:,ialt) = T;
#   m.solution{2}(:,:,ialt) = R;
#   m.solution{3}(:,:,ialt) = K;
#   m.solution{4}(:,:,ialt) = Z;
#   m.solution{5}(:,:,ialt) = H;
#   m.solution{6}(:,:,ialt) = D;
#   m.solution{7}(:,:,ialt) = U;
# =============================================================================
