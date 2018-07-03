import numpy as np

import Faces as Faces
import Vertices as Vertices
import Fluxes as Fluxes
from icosphere import icosphere
from general import find_midpoint, calc_rbf_phi, great_circle_dist, \
                        spherical_triangle_area, cart_to_polar
from general import spherical_pentagon_area
from IO import plotIcosphere

def saveIcoGrid(F,V,FL):
    path = 'grids/nSubdiv_'+str(F.nSubdiv)+'.p'
    import pickle
    f = open(path, 'wb')
    pickle.dump(F, f) 
    pickle.dump(V, f) 
    pickle.dump(FL, f) 

def loadIcoGrid(nSubdiv):
    path = 'grids/nSubdiv_'+str(nSubdiv)+'.p'
    import pickle
    f = open(path, 'rb')
    F = pickle.load(f)
    V = pickle.load(f)
    FL = pickle.load(f)
    return(F,V,FL)


def createIcoGrid(subdiv):

    ico = icosphere(scale=1,subdiv=subdiv)
    verts_cart = ico.verts
    faces = ico.faces


    rE = 6371000
    nf = len(faces) 
    nv = len(verts_cart)

    #######################
    #### PRIMARY GRID #####
    #######################

    # POINTS
    F = Faces.Faces(nf)
    F.rE = rE
    F.nSubdiv = subdiv

    for fi in range(0,len(faces)):
        xmean = 0
        ymean = 0
        zmean = 0
        for vi in faces[fi]:
            xmean = xmean + verts_cart[vi][0]
            ymean = ymean + verts_cart[vi][1]
            zmean = zmean + verts_cart[vi][2]
        xmean = xmean/3
        ymean = ymean/3
        zmean = zmean/3
        fact = 1/np.sqrt(xmean**2+ymean**2+zmean**2)
        xmean = xmean*fact
        ymean = ymean*fact
        zmean = zmean*fact
        # points_polar
        lat = np.arcsin(zmean)*(180/np.pi)
        lon = np.arctan2(ymean,xmean)*(180/np.pi) + 180
        cart_mean = [xmean,ymean,zmean]
        polar = cart_to_polar(cart_mean)
        
        # Fill Point Class Object
        F.lons[fi] = polar[0]
        F.lats[fi] = polar[1]
        F.carts[fi,:] = [xmean, ymean, zmean]
        F.Vids[fi,:] = faces[fi] 

    # VERTICES IN POLAR COORDINATES
    V = Vertices.Vertices(nv)
    for vi in range(0,nv):
        V.carts[vi,:] = verts_cart[vi]
        V.carts[vi,:] = V.carts[vi,:] 
        V.lats[vi] = np.arcsin(V.carts[vi,2])*(180/np.pi)
        V.lons[vi] = np.arctan2(V.carts[vi,1],V.carts[vi,0])*(180/np.pi) + 180

    # CALCULATE AREAS OF SPHERICAL TRIANGLES
    for fcid in F.ids:
        F.A[fcid] = rE**2*spherical_triangle_area(V.carts[F.Vids[fcid,0],:],
                                            V.carts[F.Vids[fcid,1],:],
                                            V.carts[F.Vids[fcid,2],:])

    #plt,ax = plotIcosphere(F, V, np.zeros(nf), 'H', 0, 0, 0)
    #plt.show()
    #quit()

    # FIND NEIGHBOUR POINTS
    for fi in F.ids:
        if fi % 100 == 0:
            print(fi/F.nf*100,'%')
        NFids = np.full(3,-1,np.int)
        c = 0
        for fii in F.ids[F.ids != fi]:
            if np.sum(np.in1d(F.Vids[fii], F.Vids[fi])) == 2:
                NFids[c] = fii
                c = c + 1
        F.NFids[fi,:] = NFids



    # CREATE FLUXES
    FL = Fluxes.Fluxes(F.nf)
    for fcid in F.ids:
        if fcid % 60 == 0:
            print(fcid/F.nf*100,'%')
        #print('####',fid)
        nfcids = F.NFids[fcid]
        j = 0

        for nfcid in nfcids:
            flid = FL.exists([fcid,nfcid])
            fldir = 1 # points inside
            if flid is None:
                nflon = F.lons[nfcid]
                nflat = F.lats[nfcid]
                if nflon - F.lons[fcid] < -300:
                    nflon = nflon + 360
                elif nflon - F.lons[fcid] > 300:
                    nflon = nflon - 360
                vnlon = nflon - F.lons[fcid]
                vnlat = nflat - F.lats[fcid]
                vn = [vnlon, vnlat]
                vn = vn/np.linalg.norm(vn)
                rotMat = np.array([[0, -1],[1, 0]])
                vp = np.dot(vn, rotMat)

                # vn_cart
                vn_cart = F.carts[nfcid] - F.carts[fcid]
                vn_cart = vn_cart/np.linalg.norm(vn_cart)

                # edge length (distance between vertices)
                vertids = F.Vids[fcid,np.in1d(F.Vids[fcid], F.Vids[nfcid])]
                vert_carts = V.carts[vertids,:]
                d_edge = rE*great_circle_dist(vert_carts[0,:], vert_carts[1,:])

                # NEW PERPENDICULAR VECTOR
                dlon = V.lons[vertids[1]] - V.lons[vertids[0]]
                dlat = V.lats[vertids[1]] - V.lats[vertids[0]]
                if dlon < -250:
                    dlon = dlon + 360
                elif dlon > 250:
                    dlon = dlon - 360
                vpn = [dlon, dlat]
                vpn = vpn/np.linalg.norm(vpn)

                # rotate those oriented in the wrong direction
                if np.cross(vn,vpn) > 0:
                    vpn = vpn*-1
                
                vp = vpn
                # rotate by 90 degrees to get vn
                rotMat = np.array([[0, 1],[-1, 0]])
                vn = np.dot(vp, rotMat)



                # distance between points
                d_points = rE*great_circle_dist(F.carts[fcid,:], F.carts[nfcid,:])

                # determine cartesian coordinates of flux
                flux_cart = find_midpoint(vert_carts[0,:], vert_carts[1,:])

                # determine polar coordinates of flux
                flux_polar = cart_to_polar(flux_cart)

                # weights of neighbor points 
                dist2 = great_circle_dist(F.carts[nfcid,:], flux_cart)
                dist1 = great_circle_dist(F.carts[fcid,:], flux_cart)
                totdist = dist1 + dist2
                dist1 = dist1/totdist
                dist2 = dist2/totdist
                # weights of neighbor points
                wght1 = 1 - dist1
                wght2 = 1 - dist2


                # CREATE FLUX
                flid = FL.add([fcid,nfcid],[wght1,wght2],
                        vertids,
                        flux_cart,flux_polar[0],flux_polar[1],
                        vn, vp, d_edge,d_points,
                        vn_cart)

                fldir = -1 # points oudside
            F.FLids[fcid,j] = flid
            F.FLdirs[fcid,j] = fldir
            j = j + 1

    ## FIND NEIGHBOR FLUXES IDS
    #for flid in FL.ids:
    #    Fids = FL.Fids[flid]
    #    NFLids = []
    #    NFLdirs = []
    #    for Fid in Fids:
    #        FLids,FLdirs = F.getFluxes(Fid)
    #        NFLids.extend(FLids[FLids != flid])
    #        NFLdirs.extend(FLdirs[FLids != flid])
    #    FL.NFLids[flid,0:len(NFLids)] = NFLids
    #    FL.NFLdirs[flid,0:len(NFLids)] = NFLdirs


    # STORE ARRAYS NECESSARY FOR RADIAL BASIS FUNCTION RECONSTRUCTION
    for fcid in F.ids:
        for j,flxid_j in enumerate(F.FLids[fcid]):
            F.rbf_phi0[fcid,j,:] = calc_rbf_phi(F.carts[fcid], FL.carts[flxid_j],
                    FL.vn[flxid_j]) 
            for k,flxid_k in enumerate(F.FLids[fcid]):
                F.rbf_phi[fcid,j,k] = calc_rbf_phi(FL.carts[flxid_k],
                        FL.carts[flxid_j], FL.vn[flxid_k], FL.vn[flxid_j]) 
        # invert matrix
        F.rbf_phi[fcid,:,:] = np.linalg.inv(F.rbf_phi[fcid,:,:]) 


    #######################
    ### SECONDARY GRID ####
    #######################

    # FIND FLUX IDS AND DIRECTION AS WELL AS FACE IDS
    for flxid in FL.ids:
        vrtids = FL.Vertids[flxid]
        fcids = FL.Fids[flxid]
        for vrtid in vrtids:
            V.FLids[vrtid][np.argwhere(np.isnan(V.FLids[vrtid]))[0]] = flxid
            for fcid in fcids:
                vrtfcids = V.Fids[vrtid]
                if fcid not in vrtfcids:
                    V.Fids[vrtid][np.argwhere(np.isnan(V.Fids[vrtid]))[0]] = fcid


    for vrtid in V.ids:
        vrt_area = 0
        flxids = V.FLids[vrtid,:]
        flxids = flxids[~np.isnan(flxids)].astype(np.int)
        for flxid in flxids:
            # neighbour vertice indices
            nvrtid = FL.Vertids[flxid,FL.Vertids[flxid,:] != vrtid]
            insertInd = np.argwhere(np.isnan(V.NVrtids[vrtid,:]))[0]
            # set neighbour vertice id
            V.NVrtids[vrtid,insertInd] = nvrtid 

            vrtToNvrt = np.zeros(2)
            vrtToNvrt[:] = [ V.lons[nvrtid] - V.lons[vrtid],
                        V.lats[nvrtid] - V.lats[vrtid] ]
            vrtToNvrt = vrtToNvrt/np.linalg.norm(vrtToNvrt)

            # vector pointing to neighbor vertice
            toN = V.carts[nvrtid] - V.carts[vrtid]
            toN = toN/np.linalg.norm(toN)

            # rotate it such that n_vl x t_vl = k_l
            # t_vl is the tangential component
            k_l = FL.carts[flxid]
            t_vl = - np.cross(toN,k_l)

            # find direction (contribution) of flux to vorticity
            # this means +1 = counterclockwise, -1 = clockwise
            direction = np.round(np.dot(t_vl, FL.vn_cart[flxid,:]),0)

            ## polar method of the same 
            #rotMat = np.array([[0, 1],[-1, 0]])
            #tangential = np.dot(vrtToNvrt, rotMat)

            ##direction = np.round(np.dot(tangential,FL.vn[flxid,:]))
            #direction = np.dot(tangential,FL.vn[flxid,:])
            ## TODO: is there a mistake? not all values are close to 1 or -1
            #direction = np.sign(direction)

            # insert direction to Vertice
            insertInd = np.argwhere(np.isnan(V.FLdirs[vrtid,:]))[0]
            V.FLdirs[vrtid,insertInd] = direction

            fcids = FL.Fids[flxid,:]
            fcCentres = F.carts[fcids,:]
            triang_area = spherical_triangle_area(V.carts[vrtid], fcCentres[0,:], fcCentres[1,:])
            vrt_area = vrt_area + triang_area

        V.A[vrtid] = vrt_area*rE**2




     
        



    saveIcoGrid(F,V,FL)
    return(F,V,FL)


if __name__ == '__main__':
    createIcoGrid(1)
