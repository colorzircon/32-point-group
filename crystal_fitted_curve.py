# -*- coding: utf-8 -*-
""" 
    This program is mainly used to calculate the functional relationship between energy ratio and area fraction.
    It mainly requires modules: enclosure, crystal_class, basis, crystal_lattice, matplotlib.pyplot, and crystal_shape.
"""
import numpy as np
from cgeometry import enclosure
import nlopt
import crystal_class as CC
import basis
from crystal_lattice import CrystalLattice
from crystal_shape import ps_setd_mingled
import matplotlib.pyplot as plt

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
	
def surface_energy_ratio(d, grad):
    '''
       the objective function of critical energy ratio.
    '''    
    areas = surface_areas(d)
    
    return Gamma_ratio*areas[0] + areas[1]
	
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
    opt = nlopt.opt(nlopt.LN_COBYLA, len(d0))  
    opt.set_min_objective(object_func)    
    opt.set_maxeval(Iter_max) 
    opt.set_xtol_rel(re_xtol)
    opt.set_ftol_abs(abs_ftol)
    opt.set_initial_step(dx)
    optd = opt.optimize(d0)
    minf = opt.last_optimum_value()
    enumerated_constant = opt.last_optimize_result()	
    # if enumerated_constant   == 1:
        # print 'NLOPT_SUCCESS = 1, Generic success return value.'
    # if enumerated_constant   == 2:
        # print 'NLOPT_SUCCESS = 2, stopval (above) was reached.'
    # if enumerated_constant   == 3:
        # print 'NLOPT_SUCCESS = 3, ftol_rel or ftol_abs (above) was reached.'
    # if enumerated_constant   == 4:
        # print 'NLOPT_SUCCESS = 4, xtol_rel or xtol_abs (above) was reached.'
    # if enumerated_constant   == 5:
        # print 'NLOPT_SUCCESS = 5, maxeval (above) was reached.'
    # elif enumerated_constant == 6:
        # print 'NLOPT_SUCCESS = 6, maxtime (above) was reached.' 
 
    return optd, enumerated_constant, minf

def area_fraction(Gamma_ratio):
    '''
       this function uses the distance obtained by the optimization method 
	   to calculate the area, and then gives the area fraction.
    input: energy ratio (Gamma_ratio).
	return: area fraction and areas.
    ''' 
    res   = optimal_distance_local(surface_energy_ratio, Iter_max, re_xtol, abs_ftol, d0, dx)
    areas = surface_areas(res[0])
    area_f    = np.zeros(2)
    area_f[0] = areas[0]/(areas[0] + areas[1])
    area_f[1] = areas[1]/(areas[0] + areas[1])
	
    return area_f, areas
	
def target_areas(lb1, lb2, ub1, ub2, n_point):
    '''
	area fraction of different surfaces.
    input:
        lb1, lb2, ub1, ub2: the upper and lower bounds of the area where the area fraction is located.
		n_point: number of divisions of the interval.   
		
	return: 
		ta: the given area fraction. 
    ''' 	
    ta0 = np.linspace(lb1,ub1,n_point)  
    ta1 = np.linspace(ub2,lb2,n_point)					
    ta = []	
    for i in range(len(ta0)):	
        ta.extend([[ta0[i], ta1[i]]])
    
    return ta, ta0
	
def array_to_list(qfp2):

    qfp2_list = []
    n_qfp2    = [len(p) for p in qfp2]
    for i in np.arange(len(n_qfp2)):
        for j in np.arange(n_qfp2[i]):
            qfp2_list.append(qfp2[i][j])
			
    return n_qfp2, qfp2_list 

def fitted_fun(a, x):
    '''
        this function is used to fit the energy ratios and the area fractions.
    '''
    # (1-x)*(a[0]-a[1]*np.sqrt(x/((a[0]-1)*x + 1)))
    # return 1.0/(a[0]-a[1]*np.sqrt((1-x)/((a[0]-1)*(1-x) + 1)))    
    # return (np.exp(-x) + a[0] + a[1]*np.exp(x-a[2]) - np.exp(11*x-10*a[2]-1))/(np.exp((x-a[2])*10)+1)
    return (1-x)*(a[0]-a[1]*np.sqrt(x/((a[0]-1)*x + 1))) + x*(1.0/(a[2]-a[3]*np.sqrt((1-x)/((a[2]-1)*(1-x) + 1))))
    		   
def _residual(a, x_data, y_data):

    return y_data - fitted_fun(a, x_data)

def _residual_sq(a, grad):

    r = y_data - fitted_fun(a, x_data)
	
    return sum(r*r)


if __name__=="__main__" :
    
    hkl = [[1,0,0], [0,1,2]]
    planes, number_planes = crystal_plane(hkl)
    Iter_max  = 5000 
    re_xtol   = 1e-12
    abs_ftol  = 1e-12
    lb1 = 0.00001
    ub1 = 1.0 
    lb2 = 0.0
    ub2 = 0.99999
    n_point = 21
    ta, ta0 = target_areas(lb1, lb2, ub1, ub2, n_point)
    energy_ratio = []	
	
    for i in range(len(ta)):
        print "i = ", i
        print ta[i][0], ta[i][1]
        k     = 0	
        d0 = [1, 1]      
        dx = 0.001	
        lb = 0
        ub = 2.0
        print 'lower bound = ', lb
        print 'upper bound = ', ub
		
        while (ub - lb)/ub > 1e-6:		
            Gamma_ratio = (lb + ub)/2.0
            k = k + 1
            af = area_fraction(Gamma_ratio) 
            if af[0][0] > ta[i][0]:		
                lb = Gamma_ratio
            if af[0][0] == 0:		
                ub = Gamma_ratio
            if af[0][1] == 0:		
                lb = Gamma_ratio					
            elif af[0][1] > ta[i][1]:	
                ub = Gamma_ratio
        print 'iteration number = ', k            
        mid = round((ub + lb)/2.0, 4)         
        energy_ratio.extend([mid])
	
	print energy_ratio	
	
    qfp2 = list([ta0, energy_ratio])		
    n_qfp2, qfp2_list = array_to_list(qfp2)
	
    print n_qfp2
    print qfp2_list

    data = open('data012.txt','w+') # open the file 
    data.truncate()                 # empty the file
    np.savetxt('data012.txt', qfp2_list) 
	
    x_data   = ta0  
    y_data   = np.array(energy_ratio)  # the fitted function unsupported operand type(s) for -: 'int' and 'list'.
    Iter_max = 50000
    re_xtol  = 1e-20
    re_ftol  = 1e-20
    d0 = [1, 1, 1, 1]
    dx = 0.01
    fp = optimal_distance_local(_residual_sq, Iter_max, re_xtol, re_ftol, d0, dx)	
    print 'The solution is', fp[0]
    print 'NLOPT_SUCCESS = ', fp[1]
    
    x1  = min(x_data)  
    x2  = max(x_data)
    xr  = x2-x1
    y1  = min(y_data)
    y2  = max(y_data)
    yr  = y2-y1
    x   = np.linspace(x1, x2, 500)
    y   = fitted_fun(fp[0], x)
    l1, = plt.plot(x_data, y_data, 'ro')
    l2, = plt.plot(x, y, 'b-', linewidth=3.0)

    plt.xlim(0, 1)
    plt.ylim(0, 2)  
    plt.title("{110}")
    plt.xlabel("Area Fraction, $A_{\{100\}}/(A_{\{100\}}+A_{\{110\}})$")
    plt.ylabel("Energy Ratio, $\gamma_{\{100\}}/\gamma_{\{122\}}$")
    plt.legend(handles = [l1,l2], labels = ['Actual values.', 'Fitting curve.'], loc = 'best')

#    plt.figure(figsize=(4,3))
#    plt.savefig("110.png", dpi=200)
    plt.show()


