import numpy as np
from general import interp_rbf


def calcWindTendency(F,FL,V,H,WIND,CURL,omega):
    energy = np.zeros(FL.nfl)
    rotation = np.zeros(FL.nfl)

    for flxid in FL.ids:
        fcids = FL.Fids[flxid,:]
        
        # RADIAL BASIS FUNCTION INTERPOLATION OF FLUCXES ONTO FACE CENTRE POINTS
        # ALSO: PROJECTION ONTO DIRECTION OF FLUX (normal to edge)
        wind_i = np.full((2,2), np.nan)
        #u_i = np.full(2, np.nan)
        for c in range(0,2):
            wind_i[c,:] = interp_rbf(F.rbf_phi[fcids[c]],F.rbf_phi0[fcids[c]],
                                    WIND[F.FLids[fcids[c]]])
            #u_i[c] = np.dot(wind_i[c,:], FL.vn[flxid])

        ## interpolation of flux onto edge location
        fcwghts = FL.Fcwghts[flxid,:]
        wind_l = wind_i[0,:]*fcwghts[0] + wind_i[1,:]*fcwghts[1]

        ## projection onto flux normal and tangential
        #u_l = np.dot(wind_l, FL.vn[flxid])
        #v_l = np.dot(wind_l, FL.vp[flxid])

        ## ENERGY
        K_i = 0.5*np.sum(np.power(wind_i,2),axis=1)
        H_i = np.full(2, np.nan)
        H_i[0] =  H[fcids[0]]
        H_i[1] =  H[fcids[1]]
        E_i = 9.81*H_i + K_i
        energy[flxid] = (E_i[1] - E_i[0])/FL.d_points[flxid]

        ## ROTATION 

        ## projection onto flux tangential
        v_l = np.dot(wind_l, FL.vp[flxid])

        vertids = FL.Vertids[flxid,:]
        absVort_i = CURL[vertids] + 2*np.sin(V.lats[vertids]/180*np.pi)*omega
        absVort_l = np.mean(absVort_i)
        rotation[flxid] = absVort_l*v_l


    dWINDdt = - rotation - energy
    return(dWINDdt)


#def calcWindTendency(F,FL,V,H,WIND,omega):
#    dHds = np.zeros(FL.nfl)
#    FdFds = np.zeros(FL.nfl)
#    dWINDdt = np.zeros(FL.nfl)    
#    coriolis = np.zeros(FL.nfl)    
#
#    for flxid in FL.ids:
#        fcids = FL.Fids[flxid,:]
#        
#        # 1) PRESSURE GRADIENT TERM
#        dHds[flxid] = -(H[fcids[1]] - H[fcids[0]])/FL.d_points[flxid]
#
#        # RADIAL BASIS FUNCTION INTERPOLATION OF FLUCXES ONTO FACE CENTRE POINTS
#        # ALSO: PROJECTION ONTO DIRECTION OF FLUX (normal to edge)
#        wind_i = np.full((2,2), np.nan)
#        u_i = np.full(2, np.nan)
#        for c in range(0,2):
#            wind_i[c,:] = interp_rbf(F.rbf_phi[fcids[c]],F.rbf_phi0[fcids[c]],
#                                    WIND[F.FLids[fcids[c]]])
#            u_i[c] = np.dot(wind_i[c,:], FL.vn[flxid])
#
#
#        # 2) MOMENTUM ADVECTION
#        FdFds[flxid] = WIND[flxid]*(u_i[1] - u_i[0])/FL.d_points[flxid]
#
#
#        ## 3) CORIOLIS 
#        fcwghts = FL.Fcwghts[flxid,:]
#        ## interpolation of flux onto edge location
#        wind_l = wind_i[0,:]*fcwghts[0] + wind_i[1,:]*fcwghts[1]
#
#        ## projection onto flux tangential
#        v_l = np.dot(wind_l, FL.vp[flxid])
#
#        vertids = FL.Vertids[flxid,:]
#        absVort_vert0 = 2*np.sin(V.lats[vertids[0]]/180*np.pi)*omega
#        absVort_vert1 = 2*np.sin(V.lats[vertids[1]]/180*np.pi)*omega
#        
#        absVort_flx = 0.5*(absVort_vert0 + absVort_vert1)
#        # TODO this term should have a minus sign!!!
#        coriolis[flxid] = - absVort_flx*v_l
#
#
#        ## old approach
#        #rotMat = np.array([[0, -1],[1, 0]])
#        #if FL.lats[flxid] < 0:
#        #    rotMat = -rotMat
#        ## perpendicular normal vector
#        #vp = np.dot(rotMat, FL.vn[flxid])
#
#        #wind_l = np.mean(wind_i,0)
#        #v_l = np.dot(wind_l, vp)
#        #
#        #coriolis[flxid] = 2*np.abs(np.sin(FL.lats[flxid]/180*np.pi))*omega*v_l
#
#    #print('dHds',np.max(dHds))
#    #print('FdFds',np.max(FdFds))
#    #print('coriolis',np.max(coriolis))
#
#    #dWINDdt = dHds
#    #dWINDdt = dHds - FdFds
#    dWINDdt = dHds - FdFds + coriolis
#    return(dWINDdt)
    

#def calcFluxDivergence(P,F):
#    Fdiv = np.zeros(len(P.id))
#    for id in P.id:
#        Pid = P.Pid[id]
#        Fid = P.Fid[id]
#        nNP = len(Pid[Pid >= 0])
#        Flux = np.zeros(nNP)
#        for i in range(0,nNP):
#            Flux[i] = P.Fdir[id,i]*F.val[Fid[i]]*F.ds
#        Fdiv[id] = np.sum(Flux)/P.A[id]
#    return(Fdiv)


