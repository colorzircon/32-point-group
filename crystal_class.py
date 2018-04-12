""" This module defines 32 crystal classes. The base class is CrystalClassBase; 
other classes are derived classes, named as CrystalClass1, CrystalClass2, ...
"""
import numpy as np
import types
import plane3d as plane

def verify_plane_index(hkl):
    """ Verifys wheter hkl is a valid plane index, which is a 3-int list.
    
    returns 1 if hkl is a 3-int lit, otherwise raises TypeError.
    """
    if len(hkl)!=3:
        raise TypeError("Argument hkl's length is not 3.")
    h,k,l = hkl
    if isinstance(h,int) or isinstance(k,int) or isinstance(l,int):
        return 1
    else:
        raise TypeError("hkl's length is 3, but its elements are not int.")


class CrystalClassBase:
    """ a general 3-D crystal class."""

    def __init__(self, short_name, full_name, schoenflies, schoenflies_alt, order):
        self.short_name = short_name
        self.full_name = full_name
        self.schoenflies = schoenflies
        self.schoenflies_alt = schoenflies_alt
        self.order = order

    def plane_family(self, hkl):
        """ Returns a plane family of {h k l}.

        hkl: a three-integer list or tuple, [h, k, l].
        Returns:
        a list of three-integer lists: [ (hkl1), (hkl2), ... ].
        """
        pass

class CrystalClass1(CrystalClassBase):
    """ Point group 1."""
    def __init__(self):
        CrystalClassBase.__init__(self, "1", "1", "C_1", "", 1)

    def plane_family(self, hkl):
        """ For point group 1, returns the index itself, that is [hkl].
        """
        if verify_plane_index(hkl):
            return [tuple(hkl)]
        
class CrystalClass2(CrystalClassBase):
    """ Point group -1.""" 
    def __init__(self):
        CrystalClassBase.__init__(self, "-1", "-1", "C_i", "S_2", 2)

    def plane_family(self, hkl):
        """ For -1, a plane family hkl has 2 planes: (h,k,l), (-h, -k, -l)."""
        if verify_plane_index(hkl):
            h, k, l = hkl
            return [(h,k,l), (-h, -k, -l)]

class CrystalClass3(CrystalClassBase):
    """ Point group 2."""
    def __init__(self):
        CrystalClassBase.__init__(self, "2", "2", "C_2", "", 2)

    def plane_family_c_unique(self, hkl):
        """ Assume c-axis is unique, a plane family hkl has two planes: 
            (h, k, l), (-h, -k, l) if h and k are not zero at the same time.
        """
        if verify_plane_index(hkl):
            h, k, l = hkl
            if h==0 and k==0: return [(h,k,l)]
            else: return [(h,k,l), (-h, -k, l)]

    def plane_family_b_unique(self, hkl):
        """ Assume b-axis is unique, a plane family hkl has two planes: 
            (h, k, l), (-h, k, -l) if h and l are not zero at the same time.
        """
        if verify_plane_index(hkl):
            h, k, l = hkl
            if h==0 and l==0: return [(h,k,l)]
            else: return [(h,k,l), (-h, k, -l)]

class CrystalClass4(CrystalClassBase):
    """ Point group m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "m", "m", "C_s", "C_{1h}", 2)

    def plane_family_c_unique(self, hkl):
        """ Assume c-axis is unique, a plane family hkl has two planes: 
            (h, k, l), (h, k, -l) if l is not zero.
        """
        if verify_plane_index(hkl):
            h, k, l = hkl
            if l==0: return [(h,k,l)]
            else: return [(h, k, l), (h, k, -l)]
   
    def plane_family_b_unique(self, hkl):
        """ Assume b-axis is unique, a plane family hkl has two planes: 
            (h, k, l), (h, -k, l) if k is not zero.
        """
        if verify_plane_index(hkl):
            h, k, l = hkl
            if k==0: return [(h,k,l)]
            else: return [(h, k, l), (h, -k, l)]

class CrystalClass5(CrystalClassBase):
    """ Point group 2/m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "2/m", "2/m", "C_{2h}", "", 4)

    def plane_family_c_unique(self, hkl):
        """ Assume c-axis is unique, a plane family hkl has 4 planes:
            (h,k,l), (-h, -k, l), (h, k, -l), (-h, -k, -l) if h, k, l are not
            zero at the same time.
        """
        if verify_plane_index(hkl):
            h, k, l = hkl
            if l==0:
                return [(h, k, l), (-h, -k, l)]
            elif h==0 and k==0:
                return [(h, k, l), (h, k, -l)]
            else:
                return [(h, k, l), (h, k, -l), (-h, -k, l), (-h, -k, -l)]

    def plane_family_b_unique(self, hkl):
        """ Assume b-axis is unique, a plane family hkl has 4 planes:
            (h,k,l), (h, -k, l), (-h, k, -l), (-h, -k, -l) if h, k, l are not
            zero at the same time.
        """
        if verify_plane_index(hkl):
            h, k, l = hkl
            if k==0:
                return [(h, k, l), (-h, k, -l)]
            elif h==0 and l==0:
                return [(h, k, l), (h, -k, l)]
            else:
                return [(h, k, l), (-h, k, -l), (h, -k, l), (-h, -k, -l)]

class CrystalClass6(CrystalClassBase):
    """ Point group 222."""
    def __init__(self):
        CrystalClassBase.__init__(self, "2 2 2", "2 2 2", "D_2", "V", 4)

    def plane_family(self, hkl):
        """ a plane family hkl has 4 planes:
            (h,k,l), (h, -k, -l), (-h, k, -l), (-h, -k, l) if h,k,l are not 
            zero at the same time.
        """
        if verify_plane_index(hkl):
            h, k, l = hkl
            pf = [[h, k, l]]
            if k*l != 0: pf.append((h, -k, -l))
            if h*l != 0: pf.append((-h, k, -l))
            if h*k != 0: pf.append((-h, -k, l))
            return pf


class CrystalClass7(CrystalClassBase):
    """ Point group mm2."""
    def __init__(self):
        CrystalClassBase.__init__(self, "m m 2", "m m 2", "C_{2v}", "", 4)

    def plane_family(self, hkl):
        """ A plane family hkl has 4 planes:
        (h,k,l), (-h, k, l), (h, -k, l), (-h, -k, l) if h,k are not zero
        at the same time.
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            if h!=0 and k!=0: 
                return [(h,k,l), (-h, k, l), (h, -k, l), (-h, -k, l)]
            if h==0 and k!=0:
                return [(h,k,l), (h, -k, l)] 
            if h!=0 and k==0:
                return [(h,k,l), (-h, k, l)]

class CrystalClass8(CrystalClassBase):
    """ Point group mmm."""
    def __init__(self):
        CrystalClassBase.__init__(self, "m m m", "2/m 2/m 2/m", "D_{2h}", "V_h", 8)

    def plane_family(self, hkl):
        """ A plane family hkl has 8 planes, possibly overlapping:
        (h,k,l), (-h, -k, l), (h, k, -l), (-h, k, -l), (h, -k, l),
        (h, -k, -l), (-h, k, l), (-h, -k, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [ (h, k, l), (-h, -k, l), (h, k, -l), (-h, k, -l),
                (h, -k, l), (h, -k, -l), (-h, k, l), (-h, -k, -l)]
            pf = list(set(pf_buffer))   # remove redundant plane-indexes.
            return pf

class CrystalClass9(CrystalClassBase):
    """ Point group 4."""
    def __init__(self):
        CrystalClassBase.__init__(self, "4", "4", "C_4", "", 4)

    def plane_family(self, hkl):
        """ A plane family hkl has 4 planes, possibly overlapping:
        (h,k,l), (-k, h, l), (-h, -k, l), (k, -h, l).
        """
        if verify_plane_index(hkl): 
            h,k,l = hkl
            if h==0 and k==0: return [(h,k,l)]
            else: return [(h,k,l), (-k, h, l), (-h, -k, l), (k, -h, l)]

class CrystalClass10(CrystalClassBase):
    """ Point group -4."""
    def __init__(self):
        CrystalClassBase.__init__(self, "-4", "-4", "S_4", "", 4)

    def plane_family(self, hkl):
        """ A plane family hkl has 4 planes, possibly overlapping:
        (h,k,l), (-k, h, -l), (-h, -k, l), (k, -h, -l)."""
        if verify_plane_index(hkl):
            h,k,l = hkl
            if h==0 and k==0: return [(0,0,l), (0,0,-l)]
            else: return [(h,k,l), (-k, h, -l), (-h, -k, l), (k, -h, -l)]

class CrystalClass11(CrystalClassBase):
    """ Point group 4/m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "4/m", "4/m", "C_{4h}", "", 8)

    def plane_family(self, hkl):
        """ A plane family hkl has 8 planes, possibly overlapping:
        (h,k,l),  (-k, h, l),  (-h, -k, l),  (k, -h, l).
        (h,k,-l), (-k, h, -l), (-h, -k, -l), (k, -h, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [ (h,k,l), (-k, h, l), (-h, -k, l), (k, -h, l),
                (h,k,-l), (-k, h, -l), (-h, -k, -l), (k, -h, -l),]
            pf = list(set(pf_buffer))
            return pf

class CrystalClass12(CrystalClassBase):
    """ Point group 422."""
    def __init__(self):
        CrystalClassBase.__init__(self, "4 2 2", "4 2 2", "D_4", "", 8)

    def plane_family(self, hkl):
        """ A plane family hkl has 8 planes, possibly overlapping:
        (h,k,l), (-k, h, l), (-h, -k, l), (k, -h, l),
        (-h, k, -l), (k, h, -l), (h, -k, -l), (-k, -h, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-k, h, l), (-h, -k, l), (k, -h, l),
                (-h, k, -l), (k, h, -l), (h, -k, -l), (-k, -h, -l) ]
            pf = list(set(pf_buffer))
            return pf

class CrystalClass13(CrystalClassBase):
    """ Point group 4mm."""
    def __init__(self):
        CrystalClassBase.__init__(self, "4 m m", "4 m m", "C_{4v}", "", 8)

    def plane_family(self, hkl):
        """ A plane family hkl has 8 planes, possibly overlapping:
        (h, k, l), (-k, h, l), (-h, -k, l), (k, -h, l),
        (-h, k, l), (k, h, l), (h, -k, l), (-k, -h, l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-k, h, l), (-h, -k, l), (k, -h, l),
                (-h, k, l), (k, h, l), (h, -k, l), (-k, -h, l) ]
            pf = list(set(pf_buffer))
            return pf


class CrystalClass14(CrystalClassBase):
    """ Point group -42m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "-4 2 m", "-4 2 m", "D_{2d}", "V_d", 8)

    def plane_family(self, hkl):
        """ A plane family hkl has 8 planes, possibly overlapping:
        (h, k, l), (-k, h, -l), (-h, -k, l), (k, -h, -l),
        (-k, -h, l), (k, h, l), (h, -k, -l), (-h, k, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-k, h, -l), (-h, -k, l), (k, -h, -l),
                (-k, -h, l), (k, h, l), (h, -k, -l), (-h, k, -l) ]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None

class CrystalClass15(CrystalClassBase):
    """ Point group 4/m m m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "4/m m m", "4/m 2/m 2/m", "D_{4h}", "", 16)

    def plane_family(self, hkl):
        """ A plane family hkl has 16 planes, possibly overlapping:
        (h, k, l), (-k, h, l), (-h, -k, l), (k, -h, l),
        (h, k, -l), (-k, h, -l), (-h, -k, -l), (k, -h, -l),
        (k, h, l), (-k, -h, l), (h, -k, l), (-h, k, l),
        (k, h, -l), (-k, -h, -l), (h, -k, -l), (-h, k, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h, k,l), (-k, h, l), (-h, -k, l), (k, -h, l),
                (h, k, -l), (-k, h, -l), (-h, -k, -l), (k, -h, -l),
                (k, h, l), (-k, -h, l), (h, -k, l), (-h, k, l),
                (k, h, -l), (-k, -h, -l), (h, -k, -l), (-h, k, -l) ]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None


class CrystalClass16(CrystalClassBase):
    """ Point group 3."""
    def __init__(self):
        CrystalClassBase.__init__(self, "3", "3", "C_3", "", 3)

    def plane_family(self, hkl):
        """ A plane family hkl has 3 planes, possibly overlapping:
        (h,k,l), (-(h+k), h, l), (k, -(h+k), l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h, k, l), (-(h+k), h, l), (k, -(h+k), l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None


class CrystalClass17(CrystalClassBase):
    """ Point group -3."""
    def __init__(self):
        CrystalClassBase.__init__(self, "-3", "-3", "S_6", "C_{3i}", 6)

    def plane_family(self, hkl):
        """ A plane family hkl has 6 planes, possibly overlapping:
        (h,k,l), ((h+k), -h, -l), (k, -(h+k), l),
        (-h, -k, -l), (-(h+k), h, l), (-k, (h+k), -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), ((h+k), -h, -l), (k, -(h+k), l),
                (-h, -k, -l), (-(h+k), h, l), (-k, (h+k), -l) ]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None


   
class CrystalClass18(CrystalClassBase):
    """ Point group 32."""
    def __init__(self):
        CrystalClassBase.__init__(self, "3 2", "3 2", "D_3", "", 6)

    def plane_family(self, hkl):
        """ A plane family hkl has 6 planes, possibly overlapping:
        (h,k,l), (-(h+k), h, l), (k, -(h+k), l),
        (-(h+k), k, -l), (k, h, -l), (h, -(h+k), -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-(h+k), h, l), (k, -(h+k), l),
                (-(h+k), k, -l), (k, h, -l), (h, -(h+k), -l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None


class CrystalClass19(CrystalClassBase):
    """ Point group 3m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "3 m", "3 m", "D_{3v}", "", 6)

    def plane_family(self, hkl):
        """ A plane family hkl has 6 planes, possibly overlapping:
        (h,k,l), (-(h+k), h, l), (k, -(h+k), l),
        (-k, -h, l), (-h, (h+k), l), ((h+k), -k, l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-(h+k), h, l), (k, -(h+k), l),
                (-k, -h, l), (-h, (h+k), l), ((h+k), -k, l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None


class CrystalClass21(CrystalClassBase):
    """ Point group 6."""
    def __init__(self):
        CrystalClassBase.__init__(self, "6", "6", "C_6", "", 6)

    def plane_family(self, hkl):
        """ A plane family hkl has 6 planes, possibly overlapping:
        (h,k,l), (-k, (h+k), l), (-(h+k), h, l),
        (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-k, (h+k), l), (-(h+k), h, l),
                (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None


class CrystalClass22(CrystalClassBase):
    """ Point group -6."""
    def __init__(self):
        CrystalClassBase.__init__(self, "-6", "-6", "C_{3h}", "", 6)

    def plane_family(self, hkl):
        """ A plane family hkl has 6 planes, possibly overlapping:
        (h,k,l), (k, -(h+k), -l), (-(h+k), h, l),
        (h, k, -l), (k, -(h+k), l), (-(h+k), h, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (k, -(h+k), -l), (-(h+k), h, l),
                (h, k, -l), (k, -(h+k), l), (-(h+k), h, -l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None



class CrystalClass23(CrystalClassBase):
    """ Point group 6/m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "6/m", "6/m", "C_{6h}", "", 12)

    def plane_family(self, hkl):
        """ A plane family hkl has 12 planes, possibly overlapping:
            (h,k,l), (-k, (h+k), l), (-(h+k), h, l),
            (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l),
            (h,k,-l), (-k, (h+k), -l), (-(h+k), h, -l),
            (-h, -k, -l), (k, -(h+k), -l), ((h+k), -h, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-k, (h+k), l), (-(h+k), h, l),
                (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l),
                (h,k,-l), (-k, (h+k), -l), (-(h+k), h, -l),
                (-h, -k, -l), (k, -(h+k), -l), ((h+k), -h, -l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None



class CrystalClass24(CrystalClassBase):
    """ Point group 622."""
    def __init__(self):
        CrystalClassBase.__init__(self, "622", "622", "D_{6}", "", 12)

    def plane_family(self, hkl):
        """ A plane family hkl has 12 planes, possibly overlapping:
            (h,k,l), (-k, (h+k), l), (-(h+k), h, l),
            (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l),
            (h, -(h+k),-l),  (-k, -h, -l), (-(h+k), k, -l),
            (-h, (h+k), -l), (k, h, -l),   ((h+k), -k, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-k, (h+k), l), (-(h+k), h, l),
                (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l),
                (h, -(h+k),-l),  (-k, -h, -l), (-(h+k), k, -l),
                (-h, (h+k), -l), (k, h, -l),   ((h+k), -k, -l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None



class CrystalClass25(CrystalClassBase):
    """ Point group 6mm."""
    def __init__(self):
        CrystalClassBase.__init__(self, "6 m m", "6 m m", "C_{6v}", "", 12)

    def plane_family(self, hkl):
        """ A plane family hkl has 12 planes, possibly overlapping:
            (h,k,l), (-k, (h+k), l), (-(h+k), h, l),
            (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l),
            (-h, (h+k), l), (k, h, l), ((h+k), -k, l),
            (h, -(h+k), l), (-k, -h, l),   (-(h+k), k, l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-k, (h+k), l), (-(h+k), h, l),
                (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l),
                (-h, (h+k), l), (k, h, l), ((h+k), -k, l),
                (h, -(h+k), l), (-k, -h, l),   (-(h+k), k, l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None



class CrystalClass26(CrystalClassBase):
    """ Point group -6 m 2."""
    def __init__(self):
        CrystalClassBase.__init__(self, "-6 m 2", "-6 m 2", "D_{3h}", "", 12)

    def plane_family(self, hkl):
        """ A plane family hkl has 12 planes, possibly overlapping:
            (h, k, l),  (k, -(h+k), -l), (-(h+k), h, l),
            (h, k, -l), (k, -(h+k), l),  (-(h+k), h, -l),
            (-h, (h+k), l),  (-k, -h, -l), ((h+k), -k, l),
            (-h, (h+k), -l), (-k, -h, l),  ((h+k), -k, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h, k, l),  (k, -(h+k), -l), (-(h+k), h, l),
                (h, k, -l), (k, -(h+k), l),  (-(h+k), h, -l),
                (-h, (h+k), l),  (-k, -h, -l), ((h+k), -k, l),
                (-h, (h+k), -l), (-k, -h, l),  ((h+k), -k, -l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None



class CrystalClass27(CrystalClassBase):
    """ Point group 6/m m m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "6/m m m", "6/m 2/m 2/m", "D_{6h}", "", 24)

    def plane_family(self, hkl):
        """ A plane family hkl has 24 planes, possibly overlapping:
            (h,k,l), (-k, (h+k), l), (-(h+k), h, l),
            (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l),
            (-h, (h+k), l), (k, h, l), ((h+k), -k, l),
            (h, -(h+k), l), (-k, -h, l),   (-(h+k), k, l),
            (h,k,-l), (-k, (h+k), -l), (-(h+k), h, -l),
            (-h, -k, -l), (k, -(h+k), -l), ((h+k), -h, -l),
            (-h, (h+k), -l), (k, h, -l), ((h+k), -k, -l),
            (h, -(h+k), -l), (-k, -h, -l),   (-(h+k), k, -l).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h,k,l), (-k, (h+k), l), (-(h+k), h, l),
                (-h, -k, l), (k, -(h+k), l), ((h+k), -h, l),
                (-h, (h+k), l), (k, h, l), ((h+k), -k, l),
                (h, -(h+k), l), (-k, -h, l),   (-(h+k), k, l),
                (h,k,-l), (-k, (h+k), -l), (-(h+k), h, -l),
                (-h, -k, -l), (k, -(h+k), -l), ((h+k), -h, -l),
                (-h, (h+k), -l), (k, h, -l), ((h+k), -k, -l),
                (h, -(h+k), -l), (-k, -h, -l),   (-(h+k), k, -l)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None



class CrystalClass28(CrystalClassBase):
    """ Point group 2 3."""
    def __init__(self):
        CrystalClassBase.__init__(self, "2 3", "2 3", "T", "", 12)

    def plane_family(self, hkl):
        """ A plane family hkl has 12 planes, possibly overlapping:
            (h, k, l), (h, -k, -l), (-h, -k, l), (-h, k, -l), 
            (l, h, k), (l, -h, -k), (-l, -h, k), (-l, h, -k),
            (k, l, h), (k, -l, -h), (-k, -l, h), (-k, l, -h).  
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [
                (h, k, l), (h, -k, -l), (-h, -k, l), (-h, k, -l), 
                (l, h, k), (l, -h, -k), (-l, -h, k), (-l, h, -k),
                (k, l, h), (k, -l, -h), (-k, -l, h), (-k, l, -h) ]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None



class CrystalClass29(CrystalClassBase):
    """ Point group m -3."""
    def __init__(self):
        CrystalClassBase.__init__(self, "m -3", "2/m -3", "T_h", "", 24)

    def plane_family(self, hkl):
        """ A plane family hkl has 24 planes, possibly overlapping:
            (h, k, l), (h, -k, -l), (-h, -k, l), (-h, k, -l), 
            (l, h, k), (l, -h, -k), (-l, -h, k), (-l, h, -k),
            (k, l, h), (k, -l, -h), (-k, -l, h), (-k, l, -h),  

            (-h, -k, -l), (h, k, -l), (h, -k, l), (-h, k, l),
            (-l, -h, -k), (l, h, -k), (l, -h, k), (-l, h, k),
            (-k, -l, -h), (k, l, -h), (k, -l, h), (-k, l, h).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h, k, l), (h, -k, -l), (-h, k, -l), (-h, -k, l),
                (h, k, -l), (h, -k, l), (-h, k, l), (-h, -k, -l), 
                (l, h, k), (l, -h, -k), (-l, h, -k), (-l, -h, k),
                (l, h, -k), (l, -h, k), (-l, h, k), (-l, -h, -k), 
                (k, l, h), (k, -l, -h), (-k, l, -h), (-k, -l, h),
                (k, l, -h), (k, -l, h), (-k, l, h), (-k, -l, -h),]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None


class CrystalClass30(CrystalClassBase):
    """ Point group 4 3 2."""
    def __init__(self):
        CrystalClassBase.__init__(self, "4 3 2", "4 3 2", "O", "", 24)

    def plane_family(self, hkl):
        """ A plane family hkl has 24 planes, possibly overlapping:
            (h, k, l), (h, -k, -l), (-h, -k, l), (-h, k, -l), 
            (l, h, k), (l, -h, -k), (-l, -h, k), (-l, h, -k),
            (k, l, h), (k, -l, -h), (-k, -l, h), (-k, l, -h),  

            (-k, -h, -l), (k, h, -l), (k, -h, l), (-k, h, l),
            (-h, -l, -k), (h, l, -k), (h, -l, k), (-h, l, k),
            (-l, -k, -h), (l, k, -h), (l, -k, h), (-l, k, h).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h, k, l), (h, -k, -l), (-h, -k, l), (-h, k, -l), 
                (l, h, k), (l, -h, -k), (-l, -h, k), (-l, h, -k),
                (k, l, h), (k, -l, -h), (-k, -l, h), (-k, l, -h),
                (-k, -h, -l), (k, h, -l), (k, -h, l), (-k, h, l),
                (-h, -l, -k), (h, l, -k), (h, -l, k), (-h, l, k),
                (-l, -k, -h), (l, k, -h), (l, -k, h), (-l, k, h)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None


class CrystalClass31(CrystalClassBase):
    """ Point group -4 3 m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "-4 3 m", "-4 3 m", "T_d", "", 24)

    def plane_family(self, hkl):
        """ A plane family hkl has 24 planes, possibly overlapping:
            (h, k, l), (h, -k, -l), (-h, -k, l), (-h, k, -l), 
            (l, h, k), (l, -h, -k), (-l, -h, k), (-l, h, -k),
            (k, l, h), (k, -l, -h), (-k, -l, h), (-k, l, -h),  

            (k, h, l), (k, -h, -l), (-k, h, -l), (-k, -h, l),
            (h, l, k), (h, -l, -k), (-h, l, -k), (-h, -l, k),
            (l, k, h), (l, -k, -h), (-l, k, -h), (-l, -k, h).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [(h, k, l), (h, -k, -l), (-h, -k, l), (-h, k, -l), 
                (l, h, k), (l, -h, -k), (-l, -h, k), (-l, h, -k),
                (k, l, h), (k, -l, -h), (-k, -l, h), (-k, l, -h),  
                (k, h, l), (k, -h, -l), (-k, h, -l), (-k, -h, l),
                (h, l, k), (h, -l, -k), (-h, l, -k), (-h, -l, k),
                (l, k, h), (l, -k, -h), (-l, k, -h), (-l, -k, h)]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None



class CrystalClass32(CrystalClassBase):
    """ Point group m -3 m."""
    def __init__(self):
        CrystalClassBase.__init__(self, "m -3 m", "4/m -3 2/m", "O_h", "", 48)

    def plane_family(self, hkl):
        """ A plane family hkl has 24 planes, possibly overlapping:
            (h, k, l), (h, -k, -l), (-h, -k, l), (-h, k, -l), 
            (l, h, k), (l, -h, -k), (-l, -h, k), (-l, h, -k),
            (k, l, h), (k, -l, -h), (-k, -l, h), (-k, l, -h),  

            (k, h, l), (k, -h, -l), (-k, h, -l), (-k, -h, l),
            (h, l, k), (h, -l, -k), (-h, l, -k), (-h, -l, k),
            (l, k, h), (l, -k, -h), (-l, k, -h), (-l, -k, h),

            (-h, -k, -l), (h, k, -l), (h, -k, l), (-h, k, l),
            (-l, -h, -k), (l, h, -k), (l, -h, k), (-l, h, k),
            (-k, -l, -h), (k, l, -h), (k, -l, h), (-k, l, h),

            (-k, -h, -l), (k, h, -l), (k, -h, l), (-k, h, l),
            (-h, -l, -k), (h, l, -k), (h, -l, k), (-h, l, k),
            (-l, -k, -h), (l, k, -h), (l, -k, h), (-l, k, h).
        """
        if verify_plane_index(hkl):
            h,k,l = hkl
            pf_buffer = [
                (h, k, l), (h, -k, -l), (-h, -k, l), (-h, k, -l), 
                (l, h, k), (l, -h, -k), (-l, -h, k), (-l, h, -k),
                (k, l, h), (k, -l, -h), (-k, -l, h), (-k, l, -h),  
                (k, h, l), (k, -h, -l), (-k, -h, l), (-k, h, -l),
                (h, l, k), (h, -l, -k), (-h, -l, k), (-h, l, -k),
                (l, k, h), (l, -k, -h), (-l, -k, h), (-l, k, -h),
                (-h, -k, -l), (h, k, -l), (h, -k, l), (-h, k, l),
                (-l, -h, -k), (l, h, -k), (l, -h, k), (-l, h, k),
                (-k, -l, -h), (k, l, -h), (k, -l, h), (-k, l, h),
                (-k, -h, -l), (k, h, -l), (k, -h, l), (-k, h, l),
                (-h, -l, -k), (h, l, -k), (h, -l, k), (-h, l, k),
                (-l, -k, -h), (l, k, -h), (l, -k, h), (-l, k, h),]
            pf = list(set(pf_buffer))
            return pf
        else :
            return None


# Test and debug.
if __name__=="__main__":
    cc1 = CrystalClass32()
    for index in [[1,0,0], [1,1,0], [1,1,1], [2,1,0], [3,1,1], [3,2,1],]:
        family = cc1.plane_family(index)
        print index, "has", len(family), "faces, and they are:\n", family
    print "------------------------"
    cc2 = CrystalClass28()
    for index in [[1,0,0], [1,1,0], [1,1,1], [2,1,0], [3,1,1], [3,2,1],]:
        family = cc2.plane_family(index)
        print index, "has", len(family), "faces, and they are:\n", family

