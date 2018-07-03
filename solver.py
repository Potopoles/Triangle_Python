import numpy as np
import copy
import random
import time

from heightTendency import calcHeightTendency
from fluxTendency import calcWindTendency
from diffusion import numDif
from general import gaussianHill 
from IO import makePlot, plotIcosphere, plotIcosphere_V
from sphericalGrid import createIcoGrid, loadIcoGrid
from general import wind_reconstruct


from namelist import *
if i_loadGrid:
    F,V,FL = loadIcoGrid(nSubdiv)
else:
    F,V,FL = createIcoGrid(subdiv=nSubdiv)


from general import calc_curl
print('###########')
print('###########')
print('###########')

# INITIAL CONDITIONS
CURL = np.zeros(V.nv)
CURL = calc_curl(V,FL,np.full(FL.nfl,2),CURL)

H = np.zeros(F.nf)
H[:] = 1000

WIND = np.zeros(FL.nfl)

H = gaussianHill(F,H,mag=100,lon0=70,lat0=0,siglon=30,siglat=30)

u0 = 0
if u0 != 0:
    for flid in FL.ids:
        WIND[flid] = u0*np.dot(FL.vn[flid],[1,0])
    WIND_FACE = wind_reconstruct(F,FL,WIND)


#time_TTend = 0
#time_FTend = 0
#time_numDif = 0
#time_uv = 0
#time_0 = time.time()

#ind = 0
#F.FLdirs[264,ind] = F.FLdirs[264,ind]*-1
#ind = 1
#F.FLdirs[216,ind] = F.FLdirs[216,ind]*-1


outCount = -1
simTime = 0
for t in range(0,nts):
    print('###',t, ' ', str(np.sum(H)))
    simTime = simTime + dt

    #H_old = copy.copy(H)
    #WIND_old = copy.copy(WIND)

    CURL = calc_curl(V,FL,WIND,CURL)
    #now = time.time()
    dHdt = calcHeightTendency(F,FL,H,WIND)
    #time_TTend = time_TTend + (time.time() - now)
    #now = time.time()
    dWINDdt = calcWindTendency(F,FL,V,H,WIND,CURL,omega)
    #time_FTend = time_FTend + (time.time() - now)
    H = H + dHdt*dt
    WIND = WIND + dWINDdt*dt

    #CURL = calc_curl(V,FL,WIND,CURL)
    #dHdt = calcHeightTendency(F,FL,H,WIND)
    #dWINDdt = calcWindTendency(F,FL,V,H,WIND,CURL,omega)
    #H = H_old + dHdt*dt
    #WIND = WIND_old + dWINDdt*dt

    #inpRate = 0.002
    #outRate = 0.001*1E-3
    #for fid in F.ids:
    #    H[fid] = H[fid] - dt*outRate*H[fid] + inpRate*np.cos(F.lats[fid]/180*np.pi)*dt

    #now = time.time()
    #H,WIND = numDif(F,FL,H,WIND)
    #time_numDif = time_numDif + (time.time() - now)

    thresh = 0.004
    #print('>  ' + str(thresh) + '\n' + str(np.argwhere(dHdt > thresh)))
    print('>  ' + str(thresh))
    above = np.argwhere(dHdt > thresh)
    for ai in above:
        print(str(ai) + '  lats ' + str(F.lats[ai]) + '  lons ' + str(F.lons[ai]))
    print('< -' + str(thresh))
    below = np.argwhere(dHdt < -thresh)
    for bi in below:
        print(str(bi) + '  lats ' + str(F.lats[bi]) + '  lons ' + str(F.lons[bi]))


    i_save = 1
    if (t % save_nth_ts == 0) and (t > 0):
        print('OUTPUT')

        #now = time.time()
        WIND_FACE = wind_reconstruct(F,FL,WIND)
        #time_uv = time_uv + (time.time() - now)
        #coriolis = calcMeanVal(F,FL,coriolis)

        outCount = outCount + 1
        print('outcount ', outCount)
        plt,ax = plotIcosphere(F, V, WIND_FACE[:,0], 'U_WIND', outCount, simTime, i_save)
        #plt,ax = plotIcosphere(F, V, WIND_FACE[:,1], 'V_WIND', outCount, simTime, i_save)
        plt,ax = plotIcosphere(F, V, H, 'H', outCount, simTime, i_save)
        #plt,ax = plotIcosphere_V(F, V, CURL, 'CURL', outCount, simTime, i_save)
        #dHdt[291] = 10
        #dHdt[264] = 10
        plt,ax = plotIcosphere(F, V, dHdt, 'dHdt', outCount, simTime, i_save)

        #plt.hist(dHdt)
        #plt.show()

        #plt,ax = makePlot(F, V, WIND_FACE[:,0], 'U_WIND', outCount, simTime, i_save)
        #plt,ax = makePlot(F, V, WIND_FACE[:,1], 'V_WIND', outCount, simTime, i_save)
        #plt,ax = makePlot(F, V, H, 'H', outCount, simTime, i_save)

        #time_all = time_TTend + time_FTend + time_numDif + time_uv
        #print('TTend',int(time_TTend/time_all*100))
        #print('FTend',int(time_FTend/time_all*100))
        #print('numDif',int(time_numDif/time_all*100))
        #print('uv',int(time_uv/time_all*100))
        #print('all',int(time.time() - time_0))

