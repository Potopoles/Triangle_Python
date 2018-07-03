import numpy as np
import copy
from general import interp_rbf

#def numDif(P,F,H,difDiff):
def numDif(F,FL,H,WIND):
    # TODO only works for constant side lengths
    Hdiff = copy.copy(H)
    #for fid in F.ids:
    #    NFids = F.NFids[fid]
    #    add = 0
    #    for NFid in NFids: 
    #        add = add + 0.05*(H[NFid] - H[fid])

    #    Hdiff[fid] = H[fid] + add
    #difDiff = difDiff + (Hdiff - H)
    #return(Hdiff,difDiff)

    WINDdiff = copy.copy(WIND)
    # TODO only works for constant side lengths
    for flxid in FL.ids:
        Fids = FL.Fids[flxid]

        # WIND PROJECTION
        wind_faces = np.full((2,2), np.nan)
        for i in range(0,2):
            fcid = Fids[i]
            wind_faces[i,:] = interp_rbf(F.rbf_phi[fcid],F.rbf_phi0[fcid],WIND[F.FLids[fcid]])

        # TODO USE DISTANCES
        wind_proj = np.asarray([np.dot(wind_faces[0,:], FL.vn[flxid]),
                                np.dot(wind_faces[1,:], FL.vn[flxid])])                           

        # MOMENTUM ADVECTION
        delta = (wind_proj[1] - 2*WIND[flxid] + wind_proj[0])/4
        WINDdiff[flxid] = WIND[flxid] + 0.02*delta 

    WIND = WINDdiff


    return(Hdiff,WIND)
