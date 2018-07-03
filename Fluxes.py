import numpy as np

class Fluxes:

    def __init__(self, nf):

        self.nfl = int(nf*3/2)


        self.ids = np.arange(0,self.nfl)
        self.lons = np.full(self.nfl,np.nan)
        self.lats = np.full(self.nfl,np.nan)
        self.carts = np.full((self.nfl,3),np.nan)
        # normal vector direction in polar
        self.vn = np.full((self.nfl,2),np.nan)
        # normal vector direction in cartesian
        self.vn_cart = np.full((self.nfl,3),np.nan)
        # perpendicular vector direction in polar
        self.vp = np.full((self.nfl,2),np.nan)
        # Face ids
        self.Fids = np.full((self.nfl,2),0,dtype=np.int)
        # Face weights
        self.Fcwghts = np.full((self.nfl,2),np.nan)
        # Vertice ids
        self.Vertids = np.full((self.nfl,2),0,dtype=np.int)
        ### Neighbor Flux ids
        #self.NFLids = np.full((self.nfl,4),0,dtype=np.int)
        ### Neighbor Flux directions
        #self.NFLdirs = np.full((self.nfl,4),0,dtype=np.int)

        # length of edge
        self.d_edge = np.full(self.nfl,np.nan)
        # distance between points
        self.d_points = np.full(self.nfl,np.nan)

        self.maxId = -1
        


    def exists(self,fids):
        
        inds = []
        for i in range(0,self.nfl):
            if np.array_equal(self.Fids[i,:], fids) or np.array_equal(self.Fids[i,:],
                                                        [fids[-1],fids[0]]):
                inds.append(i)

        if len(inds) > 0:
            return(inds[0])
        else:
            return(None)


    def add(self,fids,wghts,vertids,carts,lon,lat,vn,vp,d_edge,d_points,vn_cart):
        self.maxId = self.maxId + 1
        self.Fids[self.maxId,:] = fids
        self.Fcwghts[self.maxId,:] = wghts
        self.Vertids[self.maxId,:] = vertids
        self.carts[self.maxId,:] = carts 
        self.lons[self.maxId] = lon
        self.lats[self.maxId] = lat
        self.vn[self.maxId,:] = vn
        self.vp[self.maxId,:] = vp
        self.d_edge[self.maxId] = d_edge
        self.d_points[self.maxId] = d_points
        self.vn_cart[self.maxId,:] = vn_cart



        return(self.maxId)
