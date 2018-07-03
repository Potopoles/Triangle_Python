from math import sqrt
import math
import numpy as np

# -----------------------------------------------------------------------------
# Settings

#scale = 1
#subdiv = 0

class icosphere:
    middle_point_cache = {}

    def __init__(self, scale, subdiv):
        self.scale = scale
        self.subdiv = subdiv

        # -----------------------------------------------------------------------------
        # Make the base icosahedron

        # Golden ratio
        PHI = (1 + sqrt(5)) / 2

        self.verts = [
                  self.vertex(-1,  PHI, 0),
                  self.vertex( 1,  PHI, 0),
                  self.vertex(-1, -PHI, 0),
                  self.vertex( 1, -PHI, 0),

                  self.vertex(0, -1, PHI),
                  self.vertex(0,  1, PHI),
                  self.vertex(0, -1, -PHI),
                  self.vertex(0,  1, -PHI),

                  self.vertex( PHI, 0, -1),
                  self.vertex( PHI, 0,  1),
                  self.vertex(-PHI, 0, -1),
                  self.vertex(-PHI, 0,  1),
                ]

        angle = 0.97*np.pi/3
        #angle = np.pi/3
        #print(angle/np.pi*180)
        rotMat = self.rotation_matrix([0,1,0],angle)
        for i,vert in enumerate(self.verts):
            vert = np.dot(rotMat,vert)
            self.verts[i] = vert


        self.faces = [
                 # 5 faces around point 0
                 [0, 11, 5],
                 [0, 5, 1],
                 [0, 1, 7],
                 [0, 7, 10],
                 [0, 10, 11],

                 # Adjacent faces
                 [1, 5, 9],
                 [5, 11, 4],
                 [11, 10, 2],
                 [10, 7, 6],
                 [7, 1, 8],

                 # 5 faces around 3
                 [3, 9, 4],
                 [3, 4, 2],
                 [3, 2, 6],
                 [3, 6, 8],
                 [3, 8, 9],

                 # Adjacent faces
                 [4, 9, 5],
                 [2, 4, 11],
                 [6, 2, 10],
                 [8, 6, 7],
                 [9, 8, 1],
                ]


        # -----------------------------------------------------------------------------
        # Subdivisions

        for i in range(subdiv):
            faces_subdiv = []

            for tri in self.faces:
                v1 = self.middle_point(tri[0], tri[1])
                v2 = self.middle_point(tri[1], tri[2])
                v3 = self.middle_point(tri[2], tri[0])

                faces_subdiv.append([tri[0], v1, v3])
                faces_subdiv.append([tri[1], v2, v1])
                faces_subdiv.append([tri[2], v3, v2])
                faces_subdiv.append([v1, v2, v3])

            self.faces = faces_subdiv



    # -----------------------------------------------------------------------------
    # Functions



    def vertex(self, x, y, z):
        """ Return vertex coordinates fixed to the unit sphere """

        length = sqrt(x**2 + y**2 + z**2)

        return [(i * self.scale) / length for i in (x,y,z)]


    def middle_point(self, point_1, point_2):
        """ Find a middle point and project to the unit sphere """

        # We check if we have already cut this edge first
        # to avoid duplicated verts
        smaller_index = min(point_1, point_2)
        greater_index = max(point_1, point_2)

        key = '{0}-{1}'.format(smaller_index, greater_index)

        if key in self.middle_point_cache:
            return self.middle_point_cache[key]

        # If it's not in cache, then we can cut it
        vert_1 = self.verts[point_1]
        vert_2 = self.verts[point_2]
        middle = [sum(i)/2 for i in zip(vert_1, vert_2)]

        self.verts.append(self.vertex(*middle))

        index = len(self.verts) - 1
        self.middle_point_cache[key] = index

        return index
		  
    def rotation_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


