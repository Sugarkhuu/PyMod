import os
os.chdir("/home/sugarkhuu/Documents/Documents/Cprojects/Numerical and quant econ/PyDynare/PyMod")

import scipy.io
mat = scipy.io.loadmat('mydata.mat') # Sample data from IRIS:kalman_:13 - algosmoother:284

#> In kalman_ (line 13)
#  In model/filter (line 56)
#  In algo_smoother>smoother (line 284)
#  In algo_smoother (line 197)
#  In run_action_bi201909_sm_function (line 34)
#  In run_action_bi201909_sm_script (line 1)
#  In run (line 91)
#  In ogiroot/execute_action (line 75)

# size_.m
import numpy as np
ny = sum(numpy.array(mat['m']['nametype'][0][0][0] == 1))
ne = sum(numpy.array(mat['m']['nametype'][0][0][0] == 3))
npb = sum(numpy.array(mat['m']['nametype'][0][0][0] == 4))

nxb = mat['m']['solution'][0][0][0][0].shape

nx = nxb[0]
nb = nxb[1]

if len(nxb) == 2:
    nalt = 1
nf = nx - nb

nargout = 9
dopredict = nargout > 7;
dosmooth = nargout > 8;
dostoreped = nargout > 2;


eps = np.spacing(1e7)

data = mat['data']
size_data = data.shape
if len(size_data) == 2:
    ndata = 1

data = np.c_[np.full([ny,1],np.nan), data]

range = mat['range']
range = np.r_[range[0][0]-1, range[0]]
nper = range.shape[0]
ahead = mat['options']['ahead'][0][0][0][0]

npout = 0 
nloop = max([ndata,nalt])

obj = [np.nan]*nloop
varscalar = [np.nan]*nloop

F = np.full([ny,ny,nper-1,2], np.nan)
pe = np.full([ny,nper-1], np.nan);


dopredict = 1

pred = {}
pred['mean_'] = [np.full([ny,nper], np.nan), 
                np.full([nx,nper], np.nan), 
                np.full([ne,nper], np.nan), 
                range]

pred['mse_'] = [np.full([ny,ny,nper], np.nan), 
                np.full([nx,nx,nper], np.nan), 
                np.full([ne,ne,nper], np.nan), 
                range]

slastsmooth = 1
nsmooth = nper - slastsmooth + 1

smooth = {}
smooth['mean_'] = [np.full([ny,nsmooth], np.nan), 
                   np.full([nx,nsmooth], np.nan), 
                   np.full([ne,nsmooth], np.nan), 
                   range[slastsmooth-1:]]

smooth['mse_'] = [np.full([ny,ny,nsmooth], np.nan), 
                  np.full([nx,nx,nsmooth], np.nan), 
                  np.full([ne,ne,nsmooth], np.nan), 
                  range[slastsmooth-1:]]

iloop = 1


# sspace_

alt = iloop 
nalt = alt

expand = 0

m = mat['m'][0]

T = m['solution'][0][0][0]
R = m['solution'][0][0][1]
K = m['solution'][0][0][2]
Z = m['solution'][0][0][3]
H = m['solution'][0][0][4]
D = m['solution'][0][0][5]
U = m['solution'][0][0][6]


ne = sum(sum(m['nametype'][0]==3))

R = R[:,:ne]

abseigval = abs(np.linalg.eig(T[2:,:])[0])
nunit = sum(abs(abseigval-1) <= eps)

Tf   = T[:nf,:]
Rf   = R[:nf,:ne]
Ta   = T[nf:,:]
Ra   = R[nf:,:ne]
Tat  = Ta.transpose()
Tft  = Tf.transpose()
Zt   = Z.transpose()
Rat  = Ra.transpose()
Rft  = Rf.transpose()
Ht   = H.transpose()

d = m['solution'][0][0][5]
k = m['solution'][0][0][2][:]
kf = k[:nf,:]
ka = k[nf:,:]

varvec = np.square(m['assign'][0][0][-3:]).transpose()

varvec[m['shocktype'][0][0]==0] = 0

varvecindex = np.ones(data.shape) 
varvecindex[:,0] = 0

varvec = np.tile(varvec,(nper,1))
varvec = varvec.transpose()

varveceq = np.r_[False,(varvec[:,1:] == varvec[:,:-1]).all(0)]

X = np.zeros([ny, npout, nper])
D = []

D = np.zeros([ny, nper])

y1 = data

if D.size:
      y1 = y1-D

yindex = ~np.isnan(y1)
jyeq = np.r_[False,(yindex[:,1:] == yindex[:,:-1]).all(0)]


ninit = 0
dostorepredict = 1




# pedl function
nb   = Ta.shape[0]
ny   = Z.shape[0]
nper = y1.shape[1]
npout = npout
ninit = ninit
nf = Rf.shape[0]
ne = Ra.shape[0]

Ta = Ta
Tat = Tat
ka = ka
d = d
Ra = Ra
Rat = Rat
Rf = Rf
H = H
Ht = Ht
jy = np.zeros([ny,1])
#Z = Z[jy,:]
Z= np.zeros([0,8])
X = X

varvec = varvec
varveceq = varveceq
y1 = y1
yindex = yindex



#mat['mult']['std_shock_dl_cpi_xf'][0][0][0][0][1]
#
#opt ={}
#opt['modelcode'] = 'bi'
#opt['description'] = 'BI analytical forecast delivered on 2019/09/19'
#opt['idr'] = 'bi201909'
#opt['ida'] = 'sm'
