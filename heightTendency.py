import numpy as np

def calcHeightTendency(F,FL,H,WIND):
    dHdt = np.zeros(F.nf)
    for fcid in F.ids:
        #fid = 4
        FLids = F.FLids[fcid]
        NFids = F.NFids[fcid]
        #print(FLids)
        Flux = np.zeros(3)
        for i in range(0,3):
            Hmean = 0.5*(H[fcid] + H[NFids[i]])
            #print(WIND[FLids[i]])
            #print(F.FLdirs[fid,i]*WIND[FLids[i]])
            Flux[i] = Hmean*F.FLdirs[fcid,i]*WIND[FLids[i]]*FL.d_edge[FLids[i]]
        #print(np.sum(Flux))
        dHdt[fcid] = np.sum(Flux)/F.A[fcid]
        #quit()
    return(dHdt)

    

