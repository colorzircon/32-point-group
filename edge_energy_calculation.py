# -*- coding: utf-8 -*-
""" 
    This program is used to calculate the critical free energy ratios. 
    In the optimization part we use the local optimization method.
    It requires modules: enclosure, crystal_class, basis and crystal_shape.
"""
import numpy as np
from cgeometry import enclosure
import nlopt
import crystal_class as CC
import basis
from crystal_lattice import CrystalLattice
from crystal_shape import ps_setd_mingled
from scipy.optimize import minimize

def crystal_plane(hkl):
    '''
       return the geometric planes of the polyhedron, where the 32 point group is considered, 
	   and the lattice parameter is taken as 0.4nm.
	   hkl: Miller indices of crystallographic planes.
    ''' 
    b1     = basis.cubic(0.4)   # unit: m  0.4e-9
    cc1    = CC.CrystalClass32()
    # b1     = basis.hexagonal(0.4e-9,0.4e-9)   # unit: m
    # cc1    = CC.CrystalClass23()	
    latt   = CrystalLattice(b1, cc1)
    planes = []
    number_planes = []
    for index in hkl:
        pf = latt.geometric_plane_family(index)
        number_planes.append(len(pf))
        planes.extend(pf)
		
    return planes, number_planes

def one_dimension_vertex_list(planes_intersection):
    '''
       this function converts the intersection of faces, the number of faces into the form of a list.
    '''
    vertex_list = []
    n_vertex    = [len(p) for p in planes_intersection]
    for i in np.arange(len(n_vertex)):
        for j in np.arange(n_vertex[i]):
            vertex_list.append(planes_intersection[i][j])			
    return vertex_list
    
def intersection_point(number_planes, planes, d): 
    '''
       Calculate the intersection points of surfaces of a plane family.
    ''' 
    ps_setd_mingled(number_planes, planes, d)
    planes_c = [p.as_list() for p in planes]
    pr       = enclosure(planes_c)    
    points_p = pr[2]

    return points_p    
	
def vertices_number(number_planes, planes, d):
    '''
       Calculate the number of vertices of the polyhedron.
    ''' 
    points_p = intersection_point(number_planes, planes, d)
#    print 'points_p', points_p	
    for i in range(len(points_p)-1,-1,-1): # remove the empty set which repesents planes
        if points_p[i] == []:
            del points_p[i]
			
    vertex_list = one_dimension_vertex_list(points_p)	
    Pn = len(points_p)   # if calculate sum(number_planes) directly, empty sets will calculate together.
    vl = []	
    for i in vertex_list:    # remove the same vertices
      if not i in vl:
        vl.append(i)
	
    Vn = len(vl)   # number of vertices
    En = Pn+ Vn - 2    	   # number of edges
#    print 'number of planes is', Pn, ',number of vertices is', Vn, ',number of edges is', En 	
    return Vn, En, vl

def _dist_points(p1, p2):

    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)	
	
def edge_length_sort(number_planes, planes, d, vertices):

    '''
       Classify the length of edges of the polyhedron.
       1. find the find the different planes. take the last one in a set of planes.
    ''' 
    points_p = intersection_point(number_planes, planes, d)
    diff_planes = []    # the different planes.	
    plane_index = 0	
    for i in number_planes: # find the different planes, and take the last one in a set of planes.
        plane_index += i
        # if points_p[plane_index-1] != []: diff_planes.append(points_p[plane_index-1])  
        diff_planes.append(points_p[plane_index-1])   
		
    edge_length = []  # calculate the length of edge of different planes.
    for i in range(len(diff_planes)):
        if diff_planes[i] == []: # keep empty set.
            edge_length.extend([0])           
        elif diff_planes[i] != []:
            for j in range(len(diff_planes[i])-1):
                edge_length.extend([_dist_points(diff_planes[i][j], diff_planes[i][j+1])])  		         
            edge_length.extend([_dist_points(diff_planes[i][0], diff_planes[i][-1])])
        
    for i in range(len(edge_length)):  # reserved to four decimal places, because some numbers are the same, for example, 2.10123344 and 2.10123355.
        edge_length[i] = round(np.float(edge_length[i]), 4)
		
    edge_sort = []
    for i in edge_length:  # [edge_sort.append(i) for i in edge_length if not i in edge_sort]  
        if not i in edge_sort:  
            edge_sort.append(i)
           
    for i in range(len(edge_sort)): # convert [1,2] into [[1],[2]]
        edge_sort[i] = [edge_sort[i]]   
    # print edge_sort	
    
    for i in range(len(edge_sort)):   # calculate the distance between all vertices and classify them
        for j in range(len(vertices)-1):
            for k in range(j+1, len(vertices)):		
                dist = _dist_points(vertices[j], vertices[k]) 
                dist = round(dist, 4)				
                if dist == edge_sort[i][0]: 
                    edge_sort[i].append(dist)	
    # print 'edge_sort final', edge_sort                    
    return edge_sort 

	
def edge_length_fraction(edge_sort, volume, volume_calc):

    for i in range(len(edge_sort)):
        edge_sort[i]  = np.array(edge_sort[i])*np.power(volume/volume_calc, 1.0/3)
    # print 'edge_sort shrink is', edge_sort		
    edge_length_sum = 0	
    for i in range(len(edge_sort)): 
        edge_length_sum += sum(edge_sort[i][0:-1]) # edge_sort has a edge length, so take edge_sort[i][0:-1]
    # print 'edge_length_sum is', edge_length_sum
#    print len(edge_sort[0][0:-1]), len(edge_sort[1][0:-1])

    edge_length_f  = [] # fraction
    edge_length_o = [] # original	
    for i in range(len(edge_sort)):
        if list(edge_sort[i]) == list([0]): 
            edge_length_f.append(0)	
            edge_length_o.append(0)
        elif list(edge_sort[i]) != list([0]): 
            edge_length_f.append(sum(edge_sort[i][0:-1])/edge_length_sum)	
            edge_length_o.append(sum(edge_sort[i][0:-1]))
    
    return edge_length_f, edge_length_o				

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
    volume_calc = pr[0]                          # the natural or unscaled volume of the polyhedron.
    area_p      = pr[1]                          # areas of all the planes.
    points_p    = pr[2]                          # interception points on the planes, grouped by plane. 
#    print points_p                                                             	
    area_ps     = []                             # area of each plane family: sum of areas of surfaces of a plane family.               
    index_start = 0
    for subset_n in number_planes:
        index_end = index_start + subset_n
        area_ps.append(np.sum(area_p[index_start : index_end]))
        index_start = index_end  
    # diameter = 1.0e-9 
    # volume   = 4.0/3.0*np.pi* (diameter*0.5)**3   # this volume was used to scale the calculated volume.	
    area_ps  = np.array(area_ps)*np.power(volume/volume_calc, 2.0/3) 
	
    Vn, En, vl = vertices_number(number_planes, planes, d)		
    edge_sort  = edge_length_sort(number_planes, planes, d, vl)	
    edge_length_f, edge_length_o = edge_length_fraction(edge_sort, volume, volume_calc) 
 
    return area_ps, edge_length_f, edge_length_o

def area_fraction_fun(d,grad):
    '''
     	objective function of area fraction.
    '''
    print 'd is', d
    areas, edge_fraction, edge_length_o = surface_areas(d)
    print 'areas', areas, areas[0]/(areas[0]+areas[1]), areas[1]/(areas[0]+areas[1])
    print 'edge_fraction', edge_fraction
    print 'edge_length_o', edge_length_o
    areas_fra  = [areas[0]/(areas[0] + areas[1]), areas[1]/(areas[0] + areas[1])] # area fraction.
    areas_diff = np.array(area_target) - np.array(areas_fra)                      # the difference between the target area and the optimized area.
    areas_squ  = sum([x*x for x in areas_diff])                                   # the square of the area difference.
    print '-------------------------------------------------------------------'  
	
    return areas_squ
    
def area_edge_fraction_fun(d,grad):
    '''
     	objective function of area fraction.
    '''
    areas, edge_fraction, edge_length_o = surface_areas(d)
    areas_fra  = [areas[0]/(areas[0] + areas[1]), areas[1]/(areas[0] + areas[1])] # area fraction.
    areas_edge_fra = []
    areas_edge_fra = [areas_edge_fra.append(af) for af in areas_fra]
    print 'areas_edge_fra', areas_edge_fra
    areas_edge_fra = [areas_edge_fra.append(ef) for ef in edge_fraction]
    print 'areas_edge_fra', areas_edge_fra
    areas_edge_diff  = np.array(area_edge_target) - np.array(areas_edge_fra)       # the difference between the target area and the optimized area.
    areas_squ  = sum([x*x for x in areas_edge_diff])                                   # the square of the area difference.
    print '-------------------------------------------------------------------'  
	
    return areas_squ
	

def surface_energy_ratio(d, grad):
    '''
     	objective function of of critical energy ratio.
    '''
    # print 'd is', d
    # global surface_energy
    # print 'surface_energy0 is', surface_energy0
    areas, edge_fraction, edge_length_o = surface_areas(d)
    # print edge_fraction
    # areas_fraction  = [areas[0]/(areas[0] + areas[1]), areas[1]/(areas[0] + areas[1])] 
    # print 'areas_fraction is', 	areas_fraction
    # print d[0]/d[1]    
#    return surface_energy[0]*areas[0] + surface_energy[1]*areas[1] + 2.0e-10*edge_length_o[0] + 1.0e-10*edge_length_o[1]
#    return Gamma_ratio*areas[0] + areas[1] + 2.0e-10*edge_length_o[0] + 1.0e-10*edge_length_o[1]
#    return surface_energy0[0]*areas[0] + surface_energy0[1]*areas[1] + surface_energy0[2]*1.0e-9*edge_length_o[0] + surface_energy0[3]*1.0e-9*edge_length_o[1]
    # print 'areas is', areas
    # print 'edge_length_o is', edge_length_o
    if len(edge_length_o) == 2:
        return surface_energy0[0]*areas[0] + surface_energy0[1]*areas[1] + surface_energy0[2]*edge_length_o[0] + surface_energy0[3]*edge_length_o[1]
    elif len(edge_length_o) == 1:  
        return surface_energy0[0]*areas[0] + surface_energy0[1]*areas[1] + surface_energy0[3]*edge_length_o[0]
    
#    return surface_energy0[0]*areas[0] + areas[1] + surface_energy0[1]*edge_length_o[0] + surface_energy0[2]*edge_length_o[1]


def surface_energy_ratio_test(d, grad):
    '''
     	objective function of of critical energy ratio.
    '''

    areas, edge_fraction, edge_length_o = surface_areas(d)
    
    # return surface_edge_energy[0]*areas[0] + surface_edge_energy[1]*areas[1] + surface_edge_energy[2]*1.0e-9*edge_length_o[0] + surface_edge_energy[3]*1.0e-9*edge_length_o[1]
    return surface_edge_energy[0]*areas[0] + surface_edge_energy[1]*areas[1] + surface_edge_energy[2]*edge_length_o[0] + surface_edge_energy[3]*edge_length_o[1]

    
def area_edge_test(surface_energy,grad):
    '''
     	objective function of area fraction.
    '''
    print 'surface_energy is', surface_energy
    surface_energy0 = surface_energy
    global surface_energy0
    res1 = optimal_distance_local(surface_energy_ratio, 5000, 1e-6, 1e-12, [1.0, 1.0], 0.001)
    # res1 = optimal_distance_global(surface_energy_ratio, 5000, 1e-6, 1e-6, [1.0, 1.0], 0.001, [0.0, 0.0], [2.0,2.0])
    # print 'dist is', res1[0]     
    areas, edge_fraction, edge_length_o = surface_areas(res1[0]/res1[0][0]) # res1[0]  res1.x
    areas_fra = [areas[0]/(areas[0] + areas[1]), areas[1]/(areas[0] + areas[1])] # area fraction.
    print 'areas is', areas
    print 'edge length fraction is', edge_fraction 
    print 'edge length is', edge_length_o    
    areas_fra =  list(areas_fra)
    edge_fraction =  list(edge_fraction)
    if len(edge_fraction) == 2:  
        areas_fra.extend(edge_fraction)
    elif len(edge_fraction) == 1:
        areas_fra.extend([0,1])
    # print 'areas_fra', areas_fra  
    areas_edge_diff = np.array(area_edge_target) - np.array(areas_fra) 

    areas_squ  = sum([x*x for x in areas_edge_diff]) 
    print 'minf is', areas_squ    
#    print (edge_length_o[0]*surface_energy0[2]+edge_length_o[1]*surface_energy0[3])/(areas[0]*surface_energy0[0]+areas[1]*surface_energy0[1]+edge_length_o[0]*surface_energy0[2]+edge_length_o[1]*surface_energy0[3])    
    print '-------------------------------------------------------------------'  
	
    return areas_squ    
	
def optimal_distance_local(object_func, Iter_max, re_xtol, abs_ftol, d0, dx):
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
    opt = nlopt.opt(nlopt.LN_COBYLA, len(d0)) # PRAXIS  BOBYQA
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
	
def optimal_distance_global(object_func, Iter_max, re_xtol, abs_ftol, gamma0, dx, lb, ub):
    '''
     	This is the global optimization algorithm in NLopt, mainly including
		ESCH, ISRES, MLSL, etc.
    '''
    opt = nlopt.opt(nlopt.GN_ORIG_DIRECT, len(gamma0))   # G_MLSL_LDS, GN_ISRES, GN_ORIG_DIRECT, GN_ESCH, GN_CRS2_LM
    opt.set_min_objective(object_func)          
    opt.set_maxeval(Iter_max) 
    opt.set_xtol_rel(re_xtol)
    opt.set_ftol_abs(abs_ftol)
    #opt.set_ftol_rel(1e-10)
    opt.set_initial_step(dx)
	#a = opt.get_numevals()
    opt.set_lower_bounds(lb)
    opt.set_upper_bounds(ub)
    opt.set_population(100)
    opt.set_local_optimizer(nlopt.opt(nlopt.LN_BOBYQA, len(gamma0)))
    optd   = opt.optimize(gamma0)
    minf   = opt.last_optimum_value()
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
        print 'NLOPT_SUCCESS = 5, maxeval (above) was reached (global).'
    elif enumerated_constant == 6:
        print 'NLOPT_SUCCESS = 6, maxtime (above) was reached.' 
    
    return optd, enumerated_constant, minf  
    
def surface_energy_ratio_scipy(d):
    '''
     	objective function of of critical energy ratio.
    '''
    areas, edge_fraction, edge_length_o = surface_areas(d)
    
    return surface_energy0[0]*areas[0] + surface_energy0[1]*areas[1] + surface_energy0[2]*edge_length_o[0] + surface_energy0[3]*edge_length_o[1] 

def area_edge_test_scipy(surface_energy):
    '''
     	objective function of area fraction.
    '''
    print 'surface_energy is', surface_energy
    surface_energy0 = surface_energy
    global surface_energy0
    res1 = optimal_distance_local(surface_energy_ratio, 5000, 1e-6, 1e-12, [1.0, 1.0], 0.001)
    # res1 = optimal_distance_global(surface_energy_ratio, 5000, 1e-6, 1e-6, [1.0, 1.0], 0.001, [0.0, 0.0], [2.0,2.0])
    # print 'dist is', res1[0]   
    # res1 = minimize(surface_energy_ratio_scipy, [1.0, 1.0], method='nelder-mead', options={'xtol': 1e-12, 'disp': True})    
    areas, edge_fraction, edge_length_o = surface_areas(res1[0]/res1[0][0]) # res1[0]  res1.x
    areas_fra = [areas[0]/(areas[0] + areas[1]), areas[1]/(areas[0] + areas[1])] # area fraction.
    print 'areas is', areas
    print 'area fraction is', areas_fra
    print 'edge length fraction is', edge_fraction 
    print 'edge length is', edge_length_o    
    areas_fra =  list(areas_fra)
    edge_fraction =  list(edge_fraction)
    areas_fra.extend(edge_fraction)
    print 'areas_fra', areas_fra
    areas_edge_diff = np.array(area_edge_target) - np.array(areas_fra)
    print 'areas_edge_diff', areas_edge_diff     
    areas_squ  = sum([x*x for x in areas_edge_diff]) 
    print 'minf is', areas_squ    
#    print (edge_length_o[0]*surface_energy0[2]+edge_length_o[1]*surface_energy0[3])/(areas[0]*surface_energy0[0]+areas[1]*surface_energy0[1]+edge_length_o[0]*surface_energy0[2]+edge_length_o[1]*surface_energy0[3])
    print '-------------------------------------------------------------------'  
	
    return areas_squ  

if __name__=="__main__" :
	
    hkl = [[1,0,0], [1,1,1]]
    planes, number_planes = crystal_plane(hkl)
    diameter = 1.0 # 1.0e-9 
    volume   = 4.0/3.0*np.pi* (diameter*0.5)**3 
    
    # area_target = [21/(21+np.sqrt(3)), np.sqrt(3)/(21+np.sqrt(3))] 
    # res1 = optimal_distance_local(area_fraction_fun, 5000, 1e-6, 1e-6, [1.0, 1.0], 0.001)
    # print res1[0], res1[0][0]/res1[0][1], res1[0][1]/res1[0][0]
    
    # print '\n'    
    # print '============ using energy to calculate areas and edge length ==========='    
    # surface_edge_energy = [0.8, 1, 0.16, 0.12]   
    # res1 = optimal_distance_local(surface_energy_ratio_test, 5000, 1e-6, 1e-6, [1.0, 1.0], 0.001)
    # # res1 = optimal_distance_global(surface_energy_ratio_test, 8000, 1e-6, 1e-6, [1.0, 1.0], 0.1, [0.0, 0.0], [4.0, 4.0])
    # print 'surface_edge_energy is', surface_edge_energy
    # print 'd is', res1[0], res1[0][0]/res1[0][1] 	
    # area_ps, edge_length_f, edge_length_o = surface_areas(res1[0])	
    # print 'fraction of edge length is', edge_length_f 
    # print 'edge length is', edge_length_o  
    # print 'the minimum value of the function (minf) is',res1[2]     
    # print 'area fraction is', area_ps[0]/(area_ps[0]+area_ps[1]), area_ps[1]/(area_ps[0]+area_ps[1])
    # print 'area is', area_ps[0], area_ps[1]    
	
    print '\n'
    print '****************************** results *************************************'    
    lb = [0.0, 0.0, 0.0, 0.0]
    ub = [10.0, 10.0, 10.0, 10.0]
    area_edge_target = [21/(21+np.sqrt(3)), np.sqrt(3)/(21+np.sqrt(3)), 1.0/(1+np.sqrt(2)), np.sqrt(2)/(1+np.sqrt(2))]       # 
    # res2 = optimal_distance_local(area_edge_test, 5000, 1e-12, 1e-12, [0.5, 0.6, 0.5, 0.5], 0.01) 
    res2 = optimal_distance_global(area_edge_test, 5000, 1e-3, 1e-8, [5.0, 0.1, 7.0, 8.0], 0.01, lb, ub)  # 1.1295, 1.9691, 0.8906, 0.3253
    # res3 = optimal_distance_local(area_edge_test, 5000, 1e-8, 1e-14, res2[0], 0.001)
    # print res3[0], [res3[0][0]/res3[0][1], res3[0][1]/res3[0][0], res3[0][2]/res3[0][3], res3[0][3]/res3[0][2]]
    print 'surface area and edge length are', area_edge_target 
    print 'energy is', [res2[0][0], res2[0][1], res2[0][2], res2[0][3]]    
    print 'energy ratio is', [res2[0][0]/res2[0][1], res2[0][1]/res2[0][0], res2[0][2]/res2[0][3], res2[0][3]/res2[0][2]]
    print 'the minimum value of the function (minf) is', res2[2]

    # from scipy.optimize import fmin_cobyla
    # x0 = np.array([1.2, 1.0, 1.0, 1.0])
    # # res = minimize(area_edge_test_scipy, x0, method='nelder-mead', options={'xtol': 1e-10, 'disp': True})
    # res = fmin_cobyla(area_edge_test_scipy, x0, cons=[], rhobeg=1.0e-1)
    # print 'area_edge_target is', area_edge_target
    # print(res.x)
    # print res.x[0]/res.x[1], res.x[2]/res.x[3]  
    

 




    
    
