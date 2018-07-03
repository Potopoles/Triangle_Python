import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm
import scipy as sp

elev_rotation = -90
azim_rotation = 0

def makePlot(F,V,field, varName, outCount, simTime, i_save):
    polygons = []
    colors = []

    maxVal = np.max(field)

    for fid in F.ids:
        #fid = 25
        pos = np.full((3,2),np.nan)
        for i in range(0,3):
            pos[i,0] = V.lons[F.Vids[fid,i]]
            pos[i,1] = V.lats[F.Vids[fid,i]]

        lons = pos[:,0]
        #print(np.std(lons)/F.nf)
        if F.nSubdiv == 0:
            thresh = 5
        elif F.nSubdiv == 1:
            thresh = 1
        elif F.nSubdiv == 2:
            thresh = 0.1
        elif F.nSubdiv == 3:
            thresh = 0.05
        elif F.nSubdiv == 4:
            thresh = 0.006
        if np.std(lons)/F.nf > thresh:
            deltas = F.lons[fid] - lons
            condition = np.argwhere((np.abs(deltas) > np.std(lons)) & (lons != 180))
            lons[condition] = lons[condition] + (
                    np.sign(deltas)*360)[condition]

        polygon = Polygon(pos, closed=True,linewidth=0)
        polygons.append(polygon)
        #colors.append(field[fid]/maxVal)
        colors.append(field[fid])
    #pc = PatchCollection(polygons, alpha=0.8, cmap='copper')
    pc = PatchCollection(polygons, alpha=0.8)
    pc.set_edgecolor('k')
    pc.set_array(np.array(colors))

    #fig,ax = plt.subplots(figsize=(20,12))
    fig,ax = plt.subplots(figsize=(10,6))
    ax.set_xlim(-90,540)
    ax.set_ylim(-90,90)
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax)
    count = str(outCount).zfill(4)
    plt.title(count + ' ' + varName + ' day: ' + str(np.round(simTime/3600/24,2)))
    if i_save:
        name = 'out/'+varName+'_'+count+'.png'
        plt.savefig(name)
        plt.close('all')
    else:
        plt.show()

    #plt.subplots_adjust(0,0,1,1)

    return(plt,ax)


def plotIcosphere(F, V, vals, varName, outCount, simTime, i_save):
    maxVal = np.max(vals)
    fig = plt.figure(figsize=(14,10))
    ax = fig.add_subplot(111, projection='3d')
    collection = []
    colors = []
    for fid in F.ids:
        #colors.append(vals[fid]/maxVal)
        colors.append(vals[fid])
        points = []
        for vid in F.Vids[fid]:
            points.append(tuple(V.carts[vid,:]))
        collection.append(points)
    #collection = Poly3DCollection(collection, cmap='summer')
    collection = Poly3DCollection(collection)
    collection.set_edgecolor('k')
    collection.set_array(np.array(colors))
    ax.add_collection3d(collection)
    ax.set_xlim((-1,1))
    ax.set_ylim((-1,1))
    ax.set_zlim((-1,1))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(elev_rotation, azim_rotation)
    fig.colorbar(collection, ax=ax)
    fig.tight_layout()
    count = str(outCount).zfill(4)
    plt.title(count + ' ' + varName + ' day: ' + str(np.round(simTime/3600/24,2)))
    if i_save:
        name = 'out/'+varName+'_'+count+'.png'
        plt.savefig(name)
        plt.close('all')
    else:
        plt.show()

    return(plt,ax)




def plotIcosphere_V(F, V, vals, varName, outCount, simTime, i_save):
    maxVal = np.max(vals)
    fig = plt.figure(figsize=(14,10))
    ax = fig.add_subplot(111, projection='3d')
    #plt.axis('equal')
    collection = []
    colors = []
    for vrtid in V.ids:
        #vrtid = 10
        colors.append(vals[vrtid])
        points = []
        polars = []
        fcids = V.Fids[vrtid,:]
        fcids = fcids[~np.isnan(fcids)].astype(np.int)
        for fcid in fcids:
            points.append(tuple(F.carts[fcid,:]))
            polar = F.carts[fcid] - V.carts[vrtid]
            polar = polar/np.linalg.norm(polar)
            polars.append(polar)

        #print(polars)
        points = np.asarray(points)
        order = [0]
        for i in range(1,len(polars)):
            o = 0
            crossp = np.cross(polars[order[o]],polars[i])
            norm = np.linalg.norm(V.carts[vrtid] + crossp)
            while (norm <= 1) and (o < len(order)-1):
                o = o + 1
                crossp = np.cross(polars[order[o]],polars[i])
                norm = np.linalg.norm(V.carts[vrtid] + crossp)
            if o == len(order)-1:
                o = o + 1
            order.insert(o,i)

        points = points[order]
        collection.append(points)

    collection = Poly3DCollection(collection)
    collection.set_edgecolor('k')
    collection.set_array(np.array(colors))
    ax.add_collection3d(collection)
    ax.set_xlim((-1,1))
    ax.set_ylim((-1,1))
    ax.set_zlim((-1,1))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(elev_rotation, azim_rotation)
    fig.colorbar(collection, ax=ax)
    fig.tight_layout()
    count = str(outCount).zfill(4)
    plt.title(count + ' ' + varName + ' day: ' + str(np.round(simTime/3600/24,2)))
    if i_save:
        name = 'out/'+varName+'_'+count+'.png'
        plt.savefig(name)
        plt.close('all')
    else:
        plt.show()

    return(plt,ax)
