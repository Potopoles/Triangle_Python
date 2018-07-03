import numpy as np

class Faces:

    def __init__(self, nf):
        self.nf = nf

        self.A = np.full(nf,np.nan)
        #self.ds = None 
        #self.dP = None
        self.rE = None
        self.nSubdiv = None

        self.ids = np.arange(0,nf)
        self.lons = np.full(nf,np.nan)
        self.lats = np.full(nf,np.nan)
        self.carts = np.full((nf,3),np.nan)
        # Vertices ids
        self.Vids = np.full((nf,3),0,dtype=np.int)
        # Neighbour Faces ids
        self.NFids = np.full((nf,3),0,dtype=np.int)
        # Flux ids
        self.FLids = np.full((nf,3),0,dtype=np.int)
        # Flux direction (1 outward, -1 inward)
        self.FLdirs = np.full((nf,3),0,dtype=np.int)

        # radial basis function reconstruction
        self.rbf_phi = np.full((nf,3,3),np.nan)
        self.rbf_phi0 = np.full((nf,3,2),np.nan)


    def getFluxes(self,fid):
        return(self.FLids[fid,:],self.FLdirs[fid,:])

