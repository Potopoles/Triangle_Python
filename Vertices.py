import numpy as np

class Vertices:

    def __init__(self, nv):

        self.nv = nv 

        self.A = np.full(nv,np.nan)

        self.ids = np.arange(0,nv)
        self.lons = np.full(nv,np.nan)
        self.lats = np.full(nv,np.nan)
        self.carts = np.full((nv,3),np.nan)


        # Neighbour Vertices ids
        self.NVrtids = np.full((nv,6),np.nan)
        # Flux ids
        self.FLids = np.full((nv,6),np.nan)
        ## Flux direction (1 positive, -1 negative)
        self.FLdirs = np.full((nv,6),np.nan)

        ## Face ids
        self.Fids = np.full((nv,6),np.nan)


