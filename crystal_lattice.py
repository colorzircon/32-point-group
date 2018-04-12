""" This module defines the crystallographic planes (short: crystal planes)
A major task is to convert a crystal plane to a geometry plane.
It requires modules: basis, plane3d.
"""
import numpy as np
import plane3d
import basis
import crystal_class

class CrystalLattice :
    """ A crystallographic lattice. It has a basis and crystal class.

    Pay special attention to delecate settings of monoclinic class (which has 
    two settings: the default b-unique, and a non-default c-unique)."""
    def __init__(self, crystal_basis, crystal_class):
        """ initiate from a basis and a crystal class.
        """
        self.basis = crystal_basis
        self.crystal_class = crystal_class

    def geometric_plane(self, hkl):
        """ Returns a geometry plane with a three-number index.

        hkl cannot be [0,0,0], otherwise a ValueError is raised."""
        a1 = self.basis.get_vec1()
        a2 = self.basis.get_vec2()
        a3 = self.basis.get_vec3()
        h,k,l = hkl
        if h==0 and k==0 and l==0: 
            raise ValueError("Crystal plane's index cannot be (000).")
        elif h==0 and k==0 and l!=0:
            return plane3d.through_1point_2vectors(a3/l, a1, a2)
        elif h==0 and k!=0 and l==0:
            return plane3d.through_1point_2vectors(a2/k, a3, a1)
        elif h!=0 and k==0 and l==0:
            return plane3d.through_1point_2vectors(a1/h, a2, a3)
        elif h==0 and k!=0 and l!=0:
            return plane3d.through_2points_1vector(a2/k, a3/l, a1)
        elif h!=0 and k==0 and l!=0:
            return plane3d.through_2points_1vector(a3/l, a1/h, a2)
        elif h!=0 and k!=0 and l==0:
            return plane3d.through_2points_1vector(a1/h, a2/k, a3)
        else:  # h!=0, k!=0, l!=0
            return plane3d.through_3points(a1/h, a2/k, a3/l)
            
        
    def geometric_plane_family(self, hkl):
        """ Returns a list of geometry planes in a crystal plane family.

        hkl cannot be [0,0,0], otherwise a ValueError is raised."""
        planes = []
        for i in self.crystal_class.plane_family(hkl):
            planes.append(self.geometric_plane(i))
        return planes

# Test and debug
if __name__=="__main__":
    cubic_cell = basis.cubic(3.0)
    cc1 = crystal_class.CrystalClass30()
    lattice1 = CrystalLattice(cubic_cell, cc1)

    for index in [[1,0,0], [1,1,0], [1,1,1], [2,1,0], [3,1,1], [3,2,1], ]:
        family = lattice1.geometric_plane_family(index)
        print index, "has", len(family), "faces, and they are:" 
        for f in family: print f

    orthorhombic_cell = basis.orthorhombic(10.0, 3.0, 5.0)
    cc2 = crystal_class.CrystalClass8()
    lattice2 = CrystalLattice(orthorhombic_cell, cc2)
    for index in [[1,0,0], [1,1,0], [1,1,1], [2,1,0], [3,1,1], [3,2,1],]:
        family = lattice2.geometric_plane_family(index)
        print index, "has", len(family), "faces, and they are:" 
        for f in family: print f

    
