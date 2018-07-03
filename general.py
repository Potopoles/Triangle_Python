import numpy as np

def gaussianHill(F,field, mag, lon0, lat0, siglon, siglat):
    for fid in F.ids:
        field[fid] = field[fid] + mag*np.exp(
                -((F.lons[fid] - lon0)**2/(2*siglon**2) +
                  (F.lats[fid] - lat0)**2/(2*siglat**2)))
                  

    return(field)


def cart_to_polar(carts):
    """
    converts cartesian coordinates to polar

    INPUT:
    carts: array with 3 cartesian (x,y,z) coordinates 
    OUTPUT:
    polar: array with 2 polar (lon,lat) coordinates
    """
    # points_polar
    lat = np.arcsin(carts[2])*(180/np.pi)
    lon = np.arctan2(carts[1],carts[0])*(180/np.pi) + 180
    return [lon,lat]




def find_midpoint(point1, point2):
    """
    calculates midpoint of the two input points on a sphere.

    INPUT:
    point1, point2: two points defined in cartesian coordinates (x,y,z)
    """

    # find midpoint on plain
    #half_dist = (point2 - point1)/2.
    #midpoint_plain = point1 + half_dist
    midpoint_plain = (point2 + point1)/2.

    # find projection on sphere
    absolute = np.linalg.norm(midpoint_plain)
    midpoint = midpoint_plain/absolute

    return midpoint

def dist(point1,point2):
    """
    calculate magnitude of distance between two points in a cartesian coordinate system.

    INPUT:
    point1, point2: two points defined in cartesian coordinates (x,y,z). numpy arrays

    OUTPUT:
    scalar of distance
    """

    distance = np.linalg.norm(point2 - point1)

    return distance



def great_circle_dist(point1,point2):
    """
    calculate great circle distance distance between two points on the sphere in a cartesian coordinate system.
    calculations assume a unit sphere

    INPUT:
    point1, point2: two points defined in cartesian coordinates (x,y,z). numpy arrays

    OUTPUT:
    scalar of distance
    """

    distance = dist(point1, point2)
    central_angle = 2*np.arcsin(distance/2)
    gr_cir_dist = central_angle
    return gr_cir_dist


def spherical_triangle_area(pointa,pointb,pointc):
    """
    calculate area of triangle on sphere between three points on the sphere in 
    a cartesian coordinate system.
    calculations assume a unit sphere

    INPUT:
    point1, point2, point3: three points defined in cartesian coordinates (x,y,z)
                            numpy arrays

    OUTPUT:
    scalar of area
    """

    angle_A = np.arcsin(np.linalg.norm(
                np.cross(np.cross(pointa,pointb),np.cross(pointa,pointc)))/(
                np.linalg.norm(np.cross(pointa,pointb))*
                np.linalg.norm(np.cross(pointa,pointc))))
    angle_B = np.arcsin(np.linalg.norm(
                np.cross(np.cross(pointb,pointc),np.cross(pointb,pointa)))/(
                np.linalg.norm(np.cross(pointb,pointc))*
                np.linalg.norm(np.cross(pointb,pointa))))
    angle_C = np.arcsin(np.linalg.norm(
                np.cross(np.cross(pointc,pointa),np.cross(pointc,pointb)))/(
                np.linalg.norm(np.cross(pointc,pointa))*
                np.linalg.norm(np.cross(pointc,pointb))))

    # calculate spherical triangle area
    area = angle_A + angle_B + angle_C - np.pi

    return area



def spherical_pentagon_area(pointa,pointb,pointc,pointd,pointe):
    """
    calculate area of pentagon on sphere between five points on the sphere
    in a cartesian coordinate system.
    calculations assume a unit sphere

    INPUT:
    pointa, pointb, ...: five points defined in cartesian coordinates (x,y,z)
                            numpy arrays

    OUTPUT:
    scalar of area
    """


    #print(angle_A)

    norm = np.cross(pointa,pointb)
    norm_ab = norm/np.linalg.norm(norm)
    norm = np.cross(pointb,pointc)
    norm_bc = norm/np.linalg.norm(norm)
    norm = np.cross(pointc,pointd)
    norm_cd = norm/np.linalg.norm(norm)
    norm = np.cross(pointd,pointe)
    norm_de = norm/np.linalg.norm(norm)
    norm = np.cross(pointe,pointa)
    norm_ea = norm/np.linalg.norm(norm)

    mat_a = np.asarray([pointe,pointa,pointb])
    mat_b = np.asarray([pointa,pointb,pointc])
    mat_c = np.asarray([pointb,pointc,pointd])
    mat_d = np.asarray([pointc,pointd,pointe])
    mat_e = np.asarray([pointd,pointe,pointa])
    print(np.sign(np.linalg.det(mat_a)))
    print(np.sign(np.linalg.det(mat_b)))
    print(np.sign(np.linalg.det(mat_c)))
    print(np.sign(np.linalg.det(mat_d)))
    print(np.sign(np.linalg.det(mat_e)))
    angle_A = np.sign(np.linalg.det(mat_a))*np.arccos(-np.dot(norm_ea,norm_ab))
    angle_B = np.sign(np.linalg.det(mat_b))*np.arccos(-np.dot(norm_ab,norm_bc))
    angle_C = np.sign(np.linalg.det(mat_c))*np.arccos(-np.dot(norm_bc,norm_cd))
    angle_D = np.sign(np.linalg.det(mat_d))*np.arccos(-np.dot(norm_cd,norm_de))
    angle_E = np.sign(np.linalg.det(mat_e))*np.arccos(-np.dot(norm_de,norm_ea))
    #print(angle_A)
    #print(angle_B)
    #print(angle_C)
    #print(angle_D)
    #print(angle_E)
    #print()

    #angle_A = np.arcsin(np.linalg.norm(
    #            np.cross(np.cross(pointa,pointb),np.cross(pointa,pointe)))/(
    #            np.linalg.norm(np.cross(pointa,pointb))*
    #            np.linalg.norm(np.cross(pointa,pointe))))
    #angle_B = np.arcsin(np.linalg.norm(
    #            np.cross(np.cross(pointb,pointc),np.cross(pointb,pointa)))/(
    #            np.linalg.norm(np.cross(pointb,pointc))*
    #            np.linalg.norm(np.cross(pointb,pointa))))
    #angle_C = np.arcsin(np.linalg.norm(
    #            np.cross(np.cross(pointc,pointd),np.cross(pointc,pointb)))/(
    #            np.linalg.norm(np.cross(pointc,pointd))*
    #            np.linalg.norm(np.cross(pointc,pointb))))
    #angle_D = np.arcsin(np.linalg.norm(
    #            np.cross(np.cross(pointd,pointe),np.cross(pointd,pointc)))/(
    #            np.linalg.norm(np.cross(pointd,pointe))*
    #            np.linalg.norm(np.cross(pointd,pointc))))
    #angle_E = np.arcsin(np.linalg.norm(
    #            np.cross(np.cross(pointe,pointa),np.cross(pointe,pointd)))/(
    #            np.linalg.norm(np.cross(pointe,pointa))*
    #            np.linalg.norm(np.cross(pointe,pointd))))

    #print(angle_A)
    #print(angle_B)
    #print(angle_C)
    #print(angle_D)
    #print(angle_E)
    #print(angle_A + angle_B + angle_C + angle_D + angle_E)
    #print(3*np.pi)
    #quit()
    area = angle_A + angle_B + angle_C + angle_D + angle_E - 3*np.pi
    print(area)

    return area




def radial_basis_func(radius):
    epsilon = 0.5 
    return np.exp(-(radius/epsilon)**2)


def calc_rbf_phi(face, flux, norm1, norm2=None):
    """
    equation 15 or 13 of radial basis function reconstruction (Ripodas et al 2009).  

    INPUT:
    face: cartesian coordinate (x,y,z) of face mid-point
    flux: cartesian coordinate (x,y,z) of 1 flux
    norm1: normal vector of flux in lat/lon coordinates.
    norm2: normal vector needed for equation 13. if None, eq 15 is applied. 

    OUTPUT:
    phi_j/phi_k,j fill value (2 elements)
    """
    if norm2 is None:
        result = radial_basis_func(dist(face, flux))*norm1
    else: 
        result = radial_basis_func(dist(face, flux))*np.dot(norm1,norm2)

    return result


def interp_rbf(phi,phi0,wind):
    wind_face = np.dot(np.transpose(np.dot(phi, wind)), phi0)
    return wind_face


def wind_reconstruct(F,FL,WIND):
    WIND_FACE = np.full((F.nf,2), np.nan)
    for fcid in F.ids:
        WIND_FACE[fcid,:] = interp_rbf(F.rbf_phi[fcid], F.rbf_phi0[fcid], WIND[F.FLids[fcid]]) 
    return WIND_FACE


def calc_curl(V,FL,WIND,CURL):
    #vrtid = 0
    for vrtid in V.ids:
        flxids = V.FLids[vrtid,:]
        flxdirs = V.FLdirs[vrtid,:]
        flxdirs = flxdirs[~np.isnan(flxids)].astype(np.int)
        flxids = flxids[~np.isnan(flxids)].astype(np.int)

        CURL[vrtid] = np.sum(WIND[flxids]*flxdirs*FL.d_points[flxids])/V.A[vrtid]
        #print(CURL[vrtid])
    return CURL





if __name__ == '__main__':
    
    point1 = np.array([1,2,3])/np.sqrt(1**2 + 2**2 + 3**2)
    point2 = np.array([2,2,2])/np.sqrt(12)
    #point1 = np.array([1,2,3])
    #point2 = np.array([2,2,2])
    point3 = np.array([2,3,2])/np.sqrt(17)
    #midpoint = find_midpoint(point1,point2)
    #distance = dist(point1,point2)
    #print(radial_basis_func(2))
    #print(midpoint)
    #print(distance)
    norm1 = np.asarray([0.5,0.3])
    norm2 = np.asarray([0.3,0.5])
    #print(calc_phi(point1, point2, norm1))
    point1 = [-0.5, 0.8, 0.]
    point1 = [-0.8, 0.5, 0.3]
    point1 = [-0.3, 0.8, 0.5]
    print(spherical_triangle_area(point1, point2, point3))

