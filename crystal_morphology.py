# -*- coding: utf-8 -*-
""" 
    This module defines the polygons.
    It requires modules: crystal_class, basis, CrystalLattice,
    ps_setd_mingled, mlab, tvtk and enclose_polyhedron.
"""
import numpy as np
import nlopt
from cgeometry import enclosure
import crystal_class as CC
import basis
from crystal_lattice import CrystalLattice
from crystal_shape import ps_setd_mingled
from mayavi import mlab
from tvtk.api import tvtk 
from crystal_shape import enclose_polyhedron 

def crystal_plane(hkl):
    '''
       return the geometric planes of the polyhedron, where the 32 point group is considered, 
	   and the lattice parameter is taken as 0.4nm.
	   hkl: Miller indices of crystallographic planes.
    ''' 
    b1     = basis.cubic(0.4)   # unit: nm
    cc1    = CC.CrystalClass32()
    latt   = CrystalLattice(b1, cc1)
    planes = []
    number_planes = []
    for index in hkl:
        pf = latt.geometric_plane_family(index)
        number_planes.append(len(pf))
        planes.extend(pf)
    return planes, number_planes

def surface_areas(d):
    '''
    input:
        origin-to-planes distances.
		
	    volume_calc: the natural or unscaled volume of the polyhedron.
		area_p: areas of all the planes.
		points_p: interception points on the planes, grouped by plane.   
		
	return: 
		area_ps: the area of each plane family (sum of areas of surfaces of a plane family). 
    ''' 
    ps_setd_mingled(number_planes, planes, d)
    planes_c = [p.as_list() for p in planes]
    pr       = enclosure(planes_c)   
    volume_calc = pr[0]              # the natural or unscaled volume of the polyhedron.
    area_p      = pr[1]              # areas of all the planes.
    points_p    = pr[2]              # interception points on the planes, grouped by plane.                                        
    area_ps     = []                 # area of each plane family: sum of areas of surfaces of a plane family.
    index_start = 0
    for subset_n in number_planes:
        index_end   = index_start + subset_n
        area_ps.append(np.sum(area_p[index_start : index_end]))
        index_start = index_end  
    diameter = 20.0 
    volume   = 4.0/3.0*np.pi* (diameter*0.5)**3 
    area_ps  = np.array(area_ps)*np.power(volume/volume_calc, 2.0/3) 
	
    return area_ps

def object_func(d,grad):
    '''
     	objective function of area fraction.
    '''
    areas      = surface_areas(d)
    areas_fra  = [areas[0]/(areas[0] + areas[1]), areas[1]/(areas[0] + areas[1])] # area fraction.
    areas_diff = np.array(area_target) - np.array(areas_fra)                      # the difference between the target area and the optimized area.
    areas_squ  = sum([x*x for x in areas_diff])                                   # the square of the area difference.
    print 'area percentage = ', areas_fra
	
    return areas_squ

	
def optimal_distance_local(Iter_max, re_xtol, abs_ftol, d0, dx):
    ''' The NLopt (Non-Linear Optimization) library is used for the optimization problem.
	    Local optimization methods include COBYLA, BOBYQA, PRAXIS, etc. Here we use COBYLA.	
    input:
        Iter_max: the maximum number of iterations.
        re_xtol；relative tolerance on optimization parameters.
        abs_ftol：absolute tolerance on function value.	
        d0：initial value.  		
        dx: initial steps.
	return: 
		optd: the optimized values of the optimization parameters.
		enumerated_constant:  enumerated constant which takes on 1, 2, 3, 4, 5, and 6.		
		minf: the function value corresponding to optd.
    '''  
    opt = nlopt.opt(nlopt.LN_COBYLA, len(d0))  
    opt.set_min_objective(object_func)    
    opt.set_maxeval(Iter_max) 
    opt.set_xtol_rel(re_xtol)
    opt.set_ftol_abs(abs_ftol)
    opt.set_initial_step(dx)
    optd = opt.optimize(d0)
    minf = opt.last_optimum_value()
    enumerated_constant = opt.last_optimize_result()	
    if enumerated_constant   == 1:
        print 'NLOPT_SUCCESS = 1, Generic success return value.'
    if enumerated_constant   == 2:
        print 'NLOPT_SUCCESS = 2, stopval (above) was reached.'
    if enumerated_constant   == 3:
        print 'NLOPT_SUCCESS = 3, ftol_rel or ftol_abs (above) was reached.'
    if enumerated_constant   == 4:
        print 'NLOPT_SUCCESS = 4, xtol_rel or xtol_abs (above) was reached.'
    if enumerated_constant   == 5:
        print 'NLOPT_SUCCESS = 5, maxeval (above) was reached.'
    elif enumerated_constant == 6:
        print 'NLOPT_SUCCESS = 6, maxtime (above) was reached.' 
 
    return optd, enumerated_constant, minf
	
def crystal_surface_family(Mil_ind):
    ''' 
        determine all the equivalent crystal planes of the crystal plane family in the 32 point group.  
		input: 
		Mil_ind: Miller indices of crystallographic planes.
		return:
		f: crystal plane family.
    '''  	
    cc1 = CC.CrystalClass32()
    f = []
    for index in Mil_ind:
        family = cc1.plane_family(index)
        f.append(family)
		
    return f

def plot3d_surfaces(f, val, hkl, planes, number_planes):
    ''' this function is used to construct polygons. 
	    {100} in red, {111}in blue, {122} in cyan, {110} in green, 
	    {012} in gray, {113} in brown, and {123} in purple. 

    number_planes: number of planes determined by the intersection.
    planes: planes that make up a polygon.
    '''		
    cells  = tvtk.CellArray()    # create a new CellArray object to assign the polys property.
    cells.set_cells(1, planes)   # the first parameter is the number of faces (here is 1), 	                              
    p1.polys = cells             # and the second parameter is an array describing the composition of each face.
    p1.point_data.scalars = np.linspace(0.0, 1.0, len(p1.points))  
    mlab.figure(number_planes, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    if i < val[0]  and (hkl[0][0], hkl[0][1], hkl[0][2]) in f[0]:
        mlab.pipeline.surface(p1, representation='surface', opacity = 1.0, color=(1.0, 0.0, 0.0))
    if i >= val[0] and (hkl[1][0], hkl[1][1], hkl[1][2]) in f[1]:
        mlab.pipeline.surface(p1, representation='surface', opacity = 1.0, color=(0.0, 0.0, 1.0))
    if i >= val[0] and (hkl[1][0], hkl[1][1], hkl[1][2]) in f[2]:
        mlab.pipeline.surface(p1, representation='surface', opacity = 1.0, color=(0.0, 0.5, 0.5))
    if i >= val[0] and (hkl[1][0], hkl[1][1], hkl[1][2]) in f[3]:
        mlab.pipeline.surface(p1, representation='surface', opacity = 1.0, color=(0.0, 1.0, 0.0))
    if i >= val[0] and (hkl[1][0], hkl[1][1], hkl[1][2]) in f[4]:
        mlab.pipeline.surface(p1, representation='surface', opacity = 1.0, color=(0.5, 0.5, 0.5))
    if i >= val[0] and (hkl[1][0], hkl[1][1], hkl[1][2]) in f[5]:
        mlab.pipeline.surface(p1, representation='surface', opacity = 1.0, color=(0.4, 0.4, 0.0))
    if i >= val[0] and (hkl[1][0], hkl[1][1], hkl[1][2]) in f[6]:
        mlab.pipeline.surface(p1, representation='surface', opacity = 1.0, color=(0.5, 0.0, 0.5))
    # else:
        # print hkl[1], 'plane is not defined.'	
		
    mlab.pipeline.surface(p1, representation='wireframe', opacity = 1.0, color=(1, 1, 1))
#    mlab.pipeline.glyph(p1, mode='point', color=(0, 1, 0), scale_factor=0, scale_mode='none') 

    return f

def array_to_list(qfp2):
    '''
       this function converts the intersection of faces, the number of faces into the form of a list.
    '''
    qfp2_list = []
    n_qfp2    = [len(p) for p in qfp2]
    for i in np.arange(len(n_qfp2)):
        for j in np.arange(n_qfp2[i]):
            qfp2_list.append(qfp2[i][j])
			
    return n_qfp2, qfp2_list

def number_of_planes(n_qfp2):
    '''
     	calculate the number of planes that make up the polyhedron.
    '''
    val = []  
    for item in n_qfp2:
        val.append(n_qfp2.count(item))
    val[0] = 6  # the cube has 6 equivalent (100) faces.
    return  val
		

if __name__ == "__main__":  
                      		
    hkl = [[1,0,0], [1,1,3]]                    
    planes, number_planes = crystal_plane(hkl)		
    area_target = [0.8, 0.2]                    # the sum of the area fractions is 1.
    Iter_max = 5000
    re_xtol  = 1e-12
    abs_ftol = 1e-12
    d0       = [1, 1]
    dx       = 0.0001	
    res      = optimal_distance_local(Iter_max, re_xtol, abs_ftol, d0, dx)  
    areas    = surface_areas(res[0])	
    areas_fraction = [areas[0]/(areas[0] + areas[1]), areas[1]/(areas[0] + areas[1])]
    print 'The optimal distance d is', res[0], ', area is ', areas, '.'
    print 'The area fraction is ', areas_fraction,  ', target area fraction is', area_target, '.'
	
    ps_setd_mingled(number_planes, planes, res[0])  # set distances of planes. 
    diameter = 10.0      
    volume   = 4.0/3.0*np.pi* (diameter*0.5)**3
    qfp      = enclose_polyhedron(number_planes, planes, volume)

    n_qfp2, qfp2_list = array_to_list(qfp[2])
    point_for_plane = open('point_for_plane.txt','w+') # open the file. 
    point_for_plane.truncate()                         # empty the file.
    np.savetxt('point_for_plane.txt', qfp2_list)        
	
    number_of_plane = open('number_of_plane.txt','w+') 
    number_of_plane.truncate()  
    np.savetxt('number_of_plane.txt', n_qfp2)

    Mil_ind  = [[1,0,0], [1,1,1], [1,2,2], [1,1,0], [0,1,2], [1,1,3], [1,2,3]]
    f        = crystal_surface_family(Mil_ind)		
    val      = number_of_planes(n_qfp2)
    p1       = tvtk.PolyData()  
    n_planes = len(qfp[2])  
    for i in range(n_planes):   
        print "plane", i
        p1.points = qfp[2][i]                    # the coordinates of the intersection of each plane.
        n_points  = len(qfp[2][i])
        if n_points == 3: 
            faces = [n_points,0,1,2]
            plot3d_surfaces(f,val, hkl, faces, 1)
        if n_points == 4:           
            faces = [n_points,0,1,2,3]
            plot3d_surfaces(f,val, hkl, faces, 1)   
        if n_points == 5: 
            faces = [n_points,0,1,2,3,4]
            plot3d_surfaces(f,val, hkl, faces, 1)
        if n_points == 6: 
            faces = [n_points,0,1,2,3,4,5]
            plot3d_surfaces(f,val, hkl, faces, 1)  
        if n_points == 7: 
            faces = [n_points,0,1,2,3,4,5,6]
            plot3d_surfaces(f,val, hkl, faces, 1)   
        if n_points == 8: 
            faces = [n_points,0,1,2,3,4,5,6,7]
            plot3d_surfaces(f,val, hkl, faces, 1)
        if n_points == 12:
            faces = [n_points,0,1,2,3,4,5,6,7,8,9,10,11]
            plot3d_surfaces(f,val, hkl, faces, 1)
            
        elif len(qfp[2][i]) > 12:
            print "The planes have more than 12 intersection points.", len(qfp[2][i])
	axe = tvtk.AxesActor(total_length=(3,3,3))
    mlab.show()
            
