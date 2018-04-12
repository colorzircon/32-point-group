""" Shape of a particle enclosed by a set of planes.
This module invokes function 'enclosure' written in C from cgeometry.so.
The pure-python version is particle_shape_py.py.
"""
import numpy as np
from cgeometry import enclosure    # function enclosure from cgeometry.so

def ps_setd_mingled(number_planes, planes, distance) :
    ''' set distances of planes in a plane set, whose subsets are mingled 
    into one big list.

    number_planes: numbers of planes in the subsets of the plane set.
    planes: planes of all the subsets of the plane set, in a consistant oder 
        as the numbers of planes in the subsets.
    distance: [h0, h1, ...], origin-plane distances negated; can be negative or positive.
    '''
    index_start = 0
    for i, subset_n in enumerate(number_planes):
        index_end = index_start + subset_n
        for p in planes[index_start:index_end]: p.set_distance(distance[i])
        index_start = index_end    


def ps_setd_subset(planeset, distance):
    """ set distances of planes in a plane set, organized into subsets.

    planeset: a list of plane subset; a two-dimensional list.
    distance: a list of numbers; its length must equal the length of 'planeset',
        otherwise a ValueError is raised.
    """
    if(len(planeset) != len(distance)):
        raise ValueError("In particle_shape.ps_setd_mingled(): lenghth of \
            'planeset' mismatches length of 'distance'")
    for i, d in enumerate(distance):  
        for p in planeset[i]: p.set_distance(d) 


def ps_mingle(planeset):
    """ mingle multiple plane sets into a single one.

    planeset:
        a list of plane subsets.
    returns:
    number_planes:
        a list of integers; numbers of planes of each plane subset.
    planes:
        a list of planes, in the order of the argument 'planeset'.
    """
    number_planes = []
    planes = []
    for ps in planeset:
        number_planes.append(len(ps))
        planes.extend(ps)
    return [number_planes, planes]


def ps_by_family(plane_family, latt):
    """ return plane set from plane-family indices and lattice.

    plane_family: plane-family indices,
    latt: lattice parameters and point group of the crystal.
    returns: [number_planes, planes] with the default distance -1.0.
    """
    cc = latt.crystal_class
    pf_int = cc.plane_family(plane_family[0])   # a list of 3-number tuples
    number_planes = [ len(pf_int) ]
    for p in plane_family[1:]:
        if  not (p in pf_int) :
            pf = cc.plane_family(p)
            pf_int.extend( pf )
            number_planes.append( len(pf) )
        else:
            print("Duplicate plane family found:", p)
    planes = [latt.geometric_plane(p) for p in pf_int]
    return [number_planes, planes]


# def qf(nf_fs, planes, volume):
def enclose_polyhedron(number_planes, planes, volume=0.0):
    ''' uses multiple planes to enclose a polyhedron, and calculates the 
    polyhedron's geometric parameters.
    
    number_planes:
        number of planes of each subset. [list of integers]
    planes: 
        planes of all subsets. Its length equals the sum of 'number_planes'.
        [list of Plane3D]
    volume: 
        a given volume to scale the calculated volume to this value.
        If it is zero, the unscaled volume is used.

    returns: [q, fractional_area, vertices, scale].
    
    Note:
    1. This function handles one volume at a time. For multiple volumes at 
       constant height-ratios, eg. for plotting q-d curves, use 
       enclose_polyhedron_multivol(), which is faster.
    '''    
    planes_c = [p.as_list() for p in planes]
    pr = enclosure(planes_c)   # call enclosure() which is implemented in C
    volume_calc = pr[0]  # the natural or unscaled volume of the polyhedron
    area_p = pr[1]       # areas of all the planes
    points_p = pr[2]     # interception points on the planes, grouped by plane
    
    # area of each plane family: sum of areas of surfaces of a plane family 
    area_ps = []
    index_start = 0
    for subset_n in number_planes:
        index_end = index_start + subset_n
        area_ps.append(np.sum(area_p[index_start : index_end]))
        index_start = index_end
    area_ps = np.array(area_ps)   # convert a list to a numpy.ndarray
    
    # total surface area
    area_sum = np.sum(area_ps)

    # weight of each surface family
    f_ps = area_ps / area_sum    # array operation: all elements are devided.

    # surface-area-to-volume ratio
    q = area_sum / volume_calc

    # scale the surface areas, point coordinates according to volume
    if volume!=0.0:
        r = np.power(volume / volume_calc, 1.0/3)
        # correct q, f_ps, points_p, (optionally, areal_ps, area_sum)
        q = q / r 
        # f_ps remains the same
        for p in points_p:
            # p is a list of points on one plane
            for xyz in p:
                xyz[0] = xyz[0]*r
                xyz[1] = xyz[1]*r
                xyz[2] = xyz[2]*r

        # area_ps = area_ps * (r*r)  # array operation: all elements are scaled.
        # area_sum = area_sum * (r*r)
        return (q, f_ps, points_p, volume)
    else:
        # r=1.0      # volume is the natural volume_calc.
        return (q, f_ps, points_p, volume_calc)
    

def _shape_f(dist, ratio_target, number_planes, planes, volume_ref):
    """ function object being used by shape_by_f.
    
    ratio_target: the target fractional surface area.
    number_planes,
    planes,
    volume_ref.
    
    returns: sum of square of difference in ratios of all surface families.
    """
    ps_setd_mingled(number_planes, planes, dist)   # set origin-to-plane distances

    qf = enclose_polyhedron(number_planes, planes, volume_ref)
    print "target", ratio_target, "real", qf[1], "dist", dist
    t = np.array(ratio_target) - np.array(qf[1])
    t2 = sum([x*x for x in t])
    return t2

def shape_by_f(ratio, lattice, plane_family, volume_ref, distance_initial):
    """ calculates a particle's shape from ratios of surface areas.

    Arguments:
    ratio: a list of ratios of surface areas. They are dimensionless.
        Each ratio corresponds to a plane family, and its value equals to
        the total surface area of all the planes in the plane family divided 
        by the total surface area of all the planes in all the plane families.
    lattice: a CrystalLattice object; it defines the basis and crystal class 
        (or crystallographic point group).
    plane_family: a list of plane families; its length must equal to that of
        ratio.
    volume_ref: the particle's volume.
    Returns:
    a list of [q, f_ps, points_p, volume_calc, distances], where
    q: the surface-area-to-volume ratio.
    f_ps: fractional surface area of each exposed plane.
    points_p: coordinates in each exposed plane.
    distances: origin-to-plane distance of each plane family.
    """
    # import the function fmin_cobyla()
    from scipy.optimize import fmin_cobyla

    # prepare for removing duplicate plane families, e.g., a user may accidentally
    # specify an errornous plane_family = [[1,1,1], [1,0,0], [1,1,1]].
    cc = lattice.crystal_class
    pf_int = cc.plane_family(plane_family[0])   # a list of 3-number tuples
    number_planes = [ len(pf_int) ]
    for p in plane_family[1:]:
        if  not (p in pf_int) :
            pf = cc.plane_family(p)
            pf_int.extend( pf )
            number_planes.append( len(pf) )
        else:
            print("Duplicate plane family found:", p)
    planes = [lattice.geometric_plane(p) for p in pf_int]

    # set params. This is the additional argument that will be passed on to
    # fmin_cobyla()
    params = (ratio, number_planes, planes, volume_ref)
    f = fmin_cobyla(_shape_f, distance_initial, args=params, cons=[], rhobeg=1.0e-1)
    return f
    

# testing and debugging (not exhaustive)
if __name__=="__main__" :
    import crystal_class as CC
    import basis
    from crystal_lattice import CrystalLattice
    b1 = basis.cubic(0.4)   # unit: nm
    cc1 = CC.CrystalClass28()
    latt = CrystalLattice(b1, cc1)
    planes1 = []
    number_planes1 = []
    for index in [[1,2,3], [1,1,0], [2,1,0], [3,3,1], [1,0,0], [1,1,1]]:
        pf = latt.geometric_plane_family(index)
        number_planes1.append(len(pf))
        planes1.extend(pf)

    #set distances of the planes
    ps_setd_mingled(number_planes1, planes1, [2.0, 2.0, 2.0, 2.0, 2.0, 2.0])

    diameter = 10.0     # unit: nm
    volume = 4.0/3.0*np.pi* (diameter*0.5)**3

    qfp = enclose_polyhedron(number_planes1, planes1, volume)
    print "q=", qfp[0]
    print "f=", qfp[1]
    print "numbers of points are", [len(p) for p in qfp[2]]

    latt2 = CrystalLattice(b1, CC.CrystalClass32())
    f = shape_by_f([0.25, 0.75], latt2, [[0,0,1], [1,1,1]], 10.0, [1.0, 1.0])
    print "shape by f", f
    ps = ps_by_family([[0,0,1], [1,1,1]], latt2)
    ps_setd_mingled(ps[0], ps[1], f)
    qfp2 = enclose_polyhedron(ps[0], ps[1], 10.0)
    print "q2=", qfp2[0]
    print "f2=", qfp2[1]
	
	



