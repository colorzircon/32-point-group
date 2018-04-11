import numpy as np
from cgeometry import enclosure
import nlopt
import crystal_class as CC
import basis
from crystal_lattice import CrystalLattice
from crystal_shape import ps_setd_mingled

def crystal_plane(hkl):

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
		the area of each plane family (sum of areas of surfaces of a plane family). 
    '''  
    ps_setd_mingled(number_planes, planes, d)
    planes_c = [p.as_list() for p in planes]
    pr       = enclosure(planes_c)   # call enclosure() which is implemented in C
    volume_calc = pr[0]              # the natural or unscaled volume of the polyhedron
    area_p      = pr[1]              # areas of all the planes
    points_p    = pr[2]              # interception points on the planes, grouped by plane                                        
    area_ps     = []                 # area of each plane family: sum of areas of surfaces of a plane family 
    index_start = 0
    for subset_n in number_planes:
        index_end = index_start + subset_n
        area_ps.append(np.sum(area_p[index_start : index_end]))
        index_start = index_end  
    diameter = 20.0 
    volume   = 4.0/3.0*np.pi* (diameter*0.5)**3 
    area_ps  = np.array(area_ps)*np.power(volume/volume_calc, 2.0/3) 
	
    return area_ps
	

def target_func(d, grad):
    ''' target function of origin-to-planes distances.
    '''
    areas = surface_areas(d)
    
    return Gamma_ratio*areas[0] + areas[1]

	
def optimal_distance_local(Iter_max, re_xtol, abs_ftol, d0, dx):
    
    opt = nlopt.opt(nlopt.LN_COBYLA, len(d0))  # COBYLA BOBYQA PRAXIS NELDERMEAD
    opt.set_min_objective(target_func)    
    opt.set_maxeval(Iter_max) 
    opt.set_xtol_rel(re_xtol)
    opt.set_ftol_abs(abs_ftol)
    opt.set_initial_step(dx)
    optd = opt.optimize(d0)
    minf = opt.last_optimum_value()
    enumerated_constant = opt.last_optimize_result()	
    if enumerated_constant == 1:
        print 'NLOPT_SUCCESS = 1, Generic success return value.'
    if enumerated_constant == 2:
        print 'NLOPT_SUCCESS = 2, stopval (above) was reached.'
    if enumerated_constant == 3:
        print 'NLOPT_SUCCESS = 3, ftol_rel or ftol_abs (above) was reached.'
    if enumerated_constant == 4:
        print 'NLOPT_SUCCESS = 4, xtol_rel or xtol_abs (above) was reached.'
    if enumerated_constant == 5:
        print 'NLOPT_SUCCESS = 5, maxeval (above) was reached.'
    elif enumerated_constant == 6:
        print 'NLOPT_SUCCESS = 6, maxtime (above) was reached.' 
 
    return optd, enumerated_constant, minf, Gamma_ratio
	
def optimal_distance_global(Iter_max, re_xtol, abs_ftol, d0, dx, lb, ub):

    opt = nlopt.opt(nlopt.G_MLSL_LDS, len(d0))  # ESCH ISRES nlopt.G_MLSL GN_ISRES 
    opt.set_min_objective(target_func)          # G_MLSL_LDS  COBYLA
    opt.set_maxeval(Iter_max) 
    opt.set_xtol_rel(re_xtol)
    opt.set_ftol_abs(abs_ftol)
    #opt.set_ftol_rel(1e-10)
    opt.set_initial_step(dx)
	#a = opt.get_numevals()
    opt.set_lower_bounds(lb)
    opt.set_upper_bounds(ub)
    opt.set_population(10000)
    opt.set_local_optimizer(nlopt.opt(nlopt.LN_BOBYQA, 2))
    optd   = opt.optimize(d0)
    minf   = opt.last_optimum_value()
    enumerated_constant = opt.last_optimize_result()
    if enumerated_constant == 1:
        print 'NLOPT_SUCCESS = 1, Generic success return value.'
    if enumerated_constant == 2:
        print 'NLOPT_SUCCESS = 2, stopval (above) was reached.'
    if enumerated_constant == 3:
        print 'NLOPT_SUCCESS = 3, ftol_rel or ftol_abs (above) was reached.'
    if enumerated_constant == 4:
        print 'NLOPT_SUCCESS = 4, xtol_rel or xtol_abs (above) was reached.'
    if enumerated_constant == 5:
        print 'NLOPT_SUCCESS = 5, maxeval (above) was reached.'
    elif enumerated_constant == 6:
        print 'NLOPT_SUCCESS = 6, maxtime (above) was reached.' 
    
    return optd, enumerated_constant, minf, Gamma_ratio

# def opt_ratio(Iter_max, re_xtol, abs_ftol, d0, dx, lb, ub, n):
   
    # res   = optimal_distance(Iter_max, re_xtol, abs_ftol, d0, dx)
    # areas = surface_areas(res[0])   
    # print areas
		
    # if areas[0] > 0:		
        # lb = Gamma_ratio
    # elif areas[0] == 0:	
        # ub = Gamma_ratio
    # print 'iteration number = ', k 
    # mid_ub_lb = round((ub + lb)/2.0, 4) 
    # print mid_ub_lb 
    # print '--------------------------------------------'	
        


if __name__=="__main__" :
	
    # lb  = 0.0
    # ub  = 1.0
    # n   = 11
    hkl = [[1,2,3], [1,1,3]]
    planes, number_planes = crystal_plane(hkl)
	
    # for Gamma_ratio in np.linspace(lb, ub, n):
        # res   = optimal_distance(5000, 1e-12, 1e-12, [1.0, 1.0], 0.0001)
        # areas = surface_areas(res[0])
        # Gamma_ratio = round(Gamma_ratio, 5)
        # print 'Gamma_ratio = ', Gamma_ratio, ', areas = ', areas             
        # print 'optd = ', res[0]
        # print 'minf = ', res[2]
        # print '--------------------------------------------'

    lb_ub    = [[0.0, 1.0], [1.0, 2.0]]	
    Iter_max = 5000 
    re_xtol  = 1e-6
    abs_ftol = 1e-6
    d0 = [1, 1]      
    dx = 0.0001
    k  = 0	
    crit_lb_ub = np.zeros(len(lb_ub))
    opt_dist   = np.empty((len(lb_ub),len(lb_ub)), dtype=np.float)     
    opt_area   = np.empty((len(lb_ub),len(lb_ub)), dtype=np.float)
    for i in np.arange(len(lb_ub)):
        lb = lb_ub[i][0]
        ub = lb_ub[i][1]
        while abs(ub - lb)/ub > 1e-6:		   # abs(ub - lb)/ub > 1e-4  Iter_max = 5000  
            Gamma_ratio = (lb + ub)/2.0          
            k = k + 1
            res   = optimal_distance_local(Iter_max, re_xtol, abs_ftol, d0, dx)
            areas = surface_areas(res[0]) 			
            print areas
            if i == 0:                             # Calculate the lower bound of energy ratio interval.
                if areas[1] == 0:                  # if areas[1] == 0 and areas[0]>0:		
                    lb = Gamma_ratio
                elif areas[1] > 0:	               # elif areas[1] > 0  and areas[0]==0:	
                    ub = Gamma_ratio
                print 'iteration number = ', k 
                crit_lb_ub[i] = round((ub + lb)/2.0, 3) 
                opt_dist[i]   = res[0]
                opt_area[i]   = areas				
                print crit_lb_ub 				
                print '--------------------------------------------'  
   				
            elif i == 1:                           # Calculate the upper bound of energy ratio interval.
                if areas[0] > 0:		
                    lb = Gamma_ratio
                elif areas[0] == 0:	
                    ub = Gamma_ratio
                print 'iteration number = ', k 
                crit_lb_ub[i] = round((ub + lb)/2.0, 3) 
                opt_dist[i]   = res[0]
                opt_area[i]   = areas	
                print crit_lb_ub 
                print '--------------------------------------------'  
    print '****************** Critical surface energy ratios, distances and areas are as follows. *******************************'  			
    print 'The critical energy ratio interval is ', crit_lb_ub, ', lower bound is', crit_lb_ub[0], ', upper bound is', crit_lb_ub[1], '.' 
    print 'When the critical optimal distance is ', opt_dist[0], ', the area is', opt_area[0], '.'
    print 'When the critical optimal distance is ', opt_dist[1], ', the area is', opt_area[1], '.'
		
		
