""" Defines basis vectors. This include a general class Basis, and functions 
    that generate Basis objects.
"""

import numpy as np
from copy import copy

#=== helper functions
def _distance(p1, p2) :
    """distance between two 2D or 3D points."""
    if len(p1)==2 and len(p2)==2 :
        return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
    elif len(p1)==3 and len(p2)==3 :
        return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2 )
    else :
        raise ValueError
    

#=== classes 
class Basis :
    """a basis set (a 3x3 array) of a crystal lattice. """
    
    def __init__(self, v1, v2=None, v3=None):
        """ initiate from 3 3-element vectors, or a 9-element or 3x3 matrix."""
        if v2!=None:      # three 3-element vectors
            self.basis = np.array([[v1[0], v2[0], v3[0]], 
                                   [v1[1], v2[1], v3[1]],
                                   [v1[2], v2[2], v3[2]]])
        elif len(v1)==9:  # one 9-elemet array
            self.basis = np.array([[v1[0], v1[3], v1[6]], 
                                   [v1[1], v1[4], v1[7]],
                                   [v1[2], v1[5], v1[8]]])
        elif len(v1)==3 and len(v1[0])==3 :  # one 3x3 matrix
            self.basis = np.array(v1)
        else :
            raise TypeError("Incompatible argument for Basis.")

        # If the 3 vectors are co-planar, this will fail.
        self.basisI = np.array(np.matrix(self.basis).I)


    def get_vec1(self):
        """ returns the first vector as a numpy array."""
        return np.array([self.basis[0][0], self.basis[1][0], self.basis[2][0]])
    

    def get_vec2(self):
        """ returns the second vector as a numpy array."""
        return np.array([self.basis[0][1], self.basis[1][1], self.basis[2][1]])
    

    def get_vec3(self):
        """ returns the second vector as a numpy array."""
        return np.array([self.basis[0][2], self.basis[1][2], self.basis[2][2]])
    
    
    def to_cartesian(self, xyz) :
        """Return the cartian coordinate of a fractional coordinate xyz.

	xyz is a single or an array of row-major coordinates."""
        return np.dot(self.basis, xyz)


    def to_fractional(self, xyz):
        """ Return the fractional coordinate of a cartesian coordinate xyz.
        
        xyz is a single or an array of row-major coordinates."""
        return np.dot(self.basisI, xyz)
    

    def distance(self, x, y):
        """returns the distance between two fractional coordinates."""
        cv = self.to_cartesian(x) - self.to_cartesian(y)
        return np.sqrt(cv[0]**2 + cv[1]**2 + cv[2]**2)
    

    def reciprocal(self) :
        """Return the reciprocal basis.
	
        The alogrithm is [b1 b2 b3]^T = [a1 a2 a3]^{-1}. """
        return Basis(np.array(np.matrix(self.basis).I.T))


    def inversion(self) :
        """Return the inversion basis."""
        return Basis(np.array(np.matrix(self.basis).I))

    
    def transform(self, s) :
        """ Return a transformed basis.
        
        s is an array-like object or another Basis object. This method
        is useful in defining super cells or primitive cells. For example,
        the primitive cell of a face-centered cubic convention cell is
	[[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]."""
        if type(s) is Basis :
            return Basis(np.dot(self.basis, s.basis))
        else :
            return Basis(np.dot(self.basis, s))


#=== frequently used Basis
def triclinic(a, b, c, alpha, beta, gamma):
    """ A triclinic basis. a is parallel to i, a-b are in i-j plane."""
    vec1 = [a, 0.0, 0.0]
    vec2 = [b*np.cos(gamma), b*np.sin(gamma), 0.0]
    vec3 = [c*np.cos(alpha), 
        c/np.sin(gamma) * (np.cos(alpha) - np.cos(beta)*np.cos(gamma)),
        c/np.sin(gamma) * np.sqrt(np.sin(gamma)**2 - np.cos(alpha)**2
        -np.cos(betall)**2 + 2.0*np.cos(alpha)*np.cos(beta)*np.cos(gamma))]
    return Basis(vec1, vec2, vec3)


def monoclinic_b_unique(a, b, c, beta):
    """ A monoclinic basis with unique b-axis. alpha=gamma=90 degree."""
    return Basis([a, 0.0, 0.0], [0.0, b, 0.0], 
                 [c*np.cos(beta), 0.0, c*np.sin(beta)])


def monoclinic_c_unique(a, b, c, gamma):
    """ A monoclinic basis with unique c-axis. alpha=beta=90 degree."""
    return Basis([a, 0.0, 0.0], [b*np.cos(gamma), b*np.sin(gamma), 0.0], 
                 [0.0, 0.0, c])


def monoclinic(a, b, c, theta):
    """ default setting is b-unique lattice."""
    return monoclinic_b_unique(a,b,c,theta)


def orthorhombic(a,b,c):
    """ A orthorhombic basis."""
    return Basis([a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c])


def tetragonal(a,c):
    """ A tetragonal basis."""
    return Basis([a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, c])


def hexagonal(a, c):
    """ A hexagonal basis."""
    return Basis([a, 0.0, 0.0], [-a*0.5, a*np.sqrt(0.75), 0.0], [0.0, 0.0, c])


def rhombohedral_hex(a,c):
    """ A rhombohedral basis using a hexagonal lattice parameters."""
    return hexagonal(a,c)


def rhombohedral(a, theta):
    """ A rhombohedral basis by an edge length and an angle."""
    t1 = a*np.sqrt((1.0-np.cos(theta))/6.0)
    t2 = t1*np.sqrt(3.0)
    t3 = a*np.sqrt((1.0+2.0*np.cos(theta))/3.0)
    return Basis([t2, t1, t3], [-t2, t1, t3], [0.0, -t1*2.0, t3])


def trigonal(a, theta):
    """ Synonym to rhombohedral. This is defined by edge and angle."""
    return rhombohedral(a, theta)


def trigonal_hex(a,c):
    """ Synonym to rhombohedral_hex. This is defined by long and shor axis."""
    return rhombohedral_hex(a,c)


def cubic(a):
    """ A cubic basis."""
    return Basis([a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a])


#=== Tests. Not exhausted.
if __name__ == "__main__" :
    b = Basis([[5, 4, 3], [8.0, 9, 10], [1, 7, 2]])
    br = b.reciprocal()
    bi = b.inversion()
    print "b is ", b.basis
    print "br is ", br.basis
    print "bi is ", bi.basis
    print "b * br is", np.dot(b.basis, br.basis.T)
    print "b * bi is", np.dot(b.basis, bi.basis)

    xyz = [3, 4, 5]
    print "xyz is", xyz
    print "in b, xyz is", b.to_cartesian(xyz)
    print "in br, xyz is", br.to_cartesian(xyz)
    print "in bi, xyz is", bi.to_cartesian(xyz)

