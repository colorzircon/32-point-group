""" Geometry planes in three dimensions. 
"""
import numpy as np

#=== helper functions
def _det3x3(a) :
    """ determinant of a 3x3 matrix a.
    This can be also done using a numpy function numpy.linalg.det()."""
    return a[0][0] * (a[1][1]*a[2][2] - a[2][1]*a[1][2]) - a[0][1] * (
        a[1][0]*a[2][2] - a[2][0]*a[1][2]) + a[0][2] * (
        a[1][0]*a[2][1] - a[2][0]*a[1][1]) 


class Plane3D:
    """ Geometric plane represented by ax + by + cz + d = 0."""
    def __init__(self, coefficients):
        """ defines an equation of plane.
        
        coefficients is an array of [a, b, c, d] for the geometry equation
        ax + by + cz + d = 0
        Constraint: a, b, c cannot all be zero at the same time, 
            otherwise a ValueError is raised.
        Cartesian coordinate system.
        """
        a, b, c, d = coefficients[0:4]
        t = a*a + b*b + c*c
        if t == 0.0: raise ValueError
        t=np.sqrt(t)
        # normalize coef
        self.coeff = [a/t, b/t, c/t, d/t]

    def __repr__(self):
        return "Plane3D "+str(self.coeff)
    def __str__(self):
        return str(self.coeff)
        

    def distance_signed(self, p):
        """ returns the signed distance between a point and the plane.
        
        The sign is positive if the point p sits on the same side as the norm,
        and negative otherwise. Suppose p=[x, y, z], then the signed
        distance is:
        distance = (a*x + b*y + c*z + d) / sqrt(a*a + b*b + c*c)
        The plane-to-plane distance is to be implemented. 
        """
        a, b, c, d = self.coeff
        dist = a*p[0] + b*p[1] + c*p[2] + d
        return dist
    
    
    def distance(self, p):
        """ returns the distance between a point and the plane.
        
        Suppose p=[x, y, z], then the signed
        distance is
        distance = abs(a*x + b*y + c*z + d) / sqrt(a*a + b*b + c*c)
        The plane-plane distance is to be implemented. 
        """
        dist = self.distance_signed(p)
        if dist<0: 
            return -dist
        else: 
            return dist
    
    
    def as_list(self):
        """returns a list of the coefficient."""
        return self.coeff
        

    def is_parallel(self, another_plane):
        """ returns True if this plane is parallel to another, False otherwise.
        """
        t = np.cross(np.array(self.coeff[0:3]), 
                     np.array(another_plane.coeff[0:3]) )
        if t[0]*t[0] + t[1]*t[1] + t[2]*t[2] == 0.0:
            return True
        else: 
            return False
        

    def overlaps(self, another_plane):
        """ returns True if this plane is overlapping with another."""
        if self.coeff[0]==another_plane.coeff[0] and \
           self.coeff[1]==another_plane.coeff[1] and \
           self.coeff[2]==another_plane.coeff[2] and \
           self.coeff[3]==another_plane.coeff[3] :
           return True
        else:
            return False


    def set_distance(self, h):
        """ translate a plane so that the distance between origin is h.
        
        h: a number. The actual signed distance H = -h.
            If H>0 (that is h<0), the norm vector and the point are at the 
            same side of the plane. Otherwise, opposite sides.
        Note: since I want origin (that is, (0,0,0)) to be surrounded by 
        the plane, the signed distance between origin and the plane should 
        be negative. However, in most output files, the distance is a 
        positive value. So, I respect this convention. """
        #no need to reassure that a**2+b**2+c**2==1. They already are.
        #a, b, c = self.coeff[0:3]
        #self.coeff[3] = -h * np.sqrt(a*a + b*b + c*c)
        self.coeff[3] = -h
    

def by_normal_point(v, p):
    """ Defines a plane by its normal and a point lying on the plane.

    The normal (argument v) cannot be zero."""
    if v[0]==0 and v[1]==0 and v[2]==0: 
        raise ValueError("The normal vector is zero.")
    t = -np.dot(v,p)   # t is the signed origin-to-plane distance
    # let the normal vector always point out from the origin
    if (t<0.0):
        return Plane3D([v[0], v[1], v[2], t])
    else:
        return Plane3D([-v[0], -v[1], -v[2], -t])

    
def through_3points(p1, p2, p3):
    """ Returns a plane passing through three points.

    The three points must be non-collinear, and any of them must not be (0,0,0),
    otherwise a ValueError occur."""
    #The following will be checked in function by_normal_point()
    #D = _det3x3([p1,p2,p3])
    #if D==0: raise ValueError("The three points are collinear.")
    v1 = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
    v2 = [p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]]
    v = np.cross(v1,v2)
    return by_normal_point(v, p1)
    

def through_2points_1vector(p1, p2, v):
    """ Returns a plane passing through two points and one vector lying on it.

    The two points must not be parallel to the vector.
    otherwise a ValueError occur."""
    v2 = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
    n = np.cross(v2, v)
    return by_normal_point(n, p1)


def through_1point_2vectors(p1, v1, v2):
    """ Returns a plane passing through one point and two vectors lying on it.

    The two vectors must not be linearly dependent,
    otherwise a ValueError occur."""
    v = np.cross(v1, v2)
    return by_normal_point(v, p1)


# testing and debugging
if __name__=="__main__" :
    a = [1.0, 1.0, 1.0, -1.0]
    p1 = Plane3D(a)
    print p1

    p2 = through_2points_1vector([0.0, 1.0, 1.0], [2.0, 0.0, 0.0], 
        [1.0, -1.0, 0.0])
    print p2

    p3 = through_1point_2vectors([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 
        [0.0, 0.0, 1.0])
    print p3
