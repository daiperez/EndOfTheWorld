""" rootfind.py -- library of rootfinding routines
     
    Language: Python 3
    Mark A. Caprio
    University of Notre Dame
    Written for Computational Methods in Physics, Spring 2014.
"""

def bisection(f,interval,tolerance,verbose=False):
    """ Find root by bisection.

    The 'approximation' x_i at each iteration is defined by the
    midpoint of the interval.
    
    The 'error' x_i-x_(i-1) is defined by the change in midpoint from
    the midpoint of the last interval.  (Of course, for bisection,
    that is always half the width of the new interval.)

    Returns None if the sign of the function does not change on the
    given interval.  Otherwise, returns final midpoint x_i when
    termination condition is reached.

    f: function for rootfinding
    interval: tuple containing initial interval endpoints (xa,xb)
    tolerance: difference x_i-x_(i-1) at which search should terminate
    verbose (optional): whether or not to print iteration log
    """

    # set up initial bracketing interval
    #   Note: Sign of function *must* change in this interval for method to work.
    (xa,xb) = interval
    fxa = f(xa)
    fxb = f(xb)
    if (fxa*fxb >=0):
        # no sign change in interval
        return None

    # set up for first iteration
    xm = (xb + xa)/2
    error = (xb - xa)/2
    iteration_count = 0

    # bisect until tolerance reached
    while (abs(error) > tolerance):

        # increment iteration count
        iteration_count += 1
        
        # evaluate function
        fxa = f(xa)
        fxb = f(xb)
        fxm = f(xm)

        # find which subinterval contains root
        if (fxm == 0):
            # accidentally landed on root (often occurs for "toy" test intervals)
            xa = xm
            xb = xm
        elif ((fxa * fxm) < 0):
            # sign change is in left half of interval
            xb = xm
        else:
            # sign change is in right half of interval
            xa = xm

        # find new midpoint (and change in midpoint)
        xm_old = xm
        xm = (xb + xa)/2
        error = xm - xm_old

        # verbose iteration log
        if (verbose):
            print("iteration", iteration_count, "(bisection):",
                  "interval", (xa, xb), "root", xm)
            
    return xm


"""
ODE Solve function
Written by Jason
"""

def solve(f,y0,interval,steps,order=1):
    """ Solve ODE by Euler or Runge-Kutta methods, with fixed number
    of steps.

    In contrast to the examples of Newman Chapter 8, which build up a
    list, point by point, 
    
    f: function giving ODE as y'=f(x,y)
    y0: initial value
    interval: tuple region (a,b) on which to solve ODE
    steps: number of steps
    order: order of solution method (1 for Euler, 2 or 4 for Runge-Kutta)
    
    Returns (x,y) points, as (steps+1)x2 numpy array.
    """
    
    from math import sin
    from numpy import arange
    from pylab import plot, xlabel,ylabel, show
    
    a,b = interval
    h = (b - a)/steps
    x = 0.0
    t_points = arange(a,b,h)
    x_points=[]
    
    # solve Euler method - source code found in Newman section 8.11
    if order == 1:
      
        
        for t in t_points:
            x_points.append(x)
            x += h*f(x,t)

        
       #print(t_points,x_points)
    
    # solve range-kutta 2nd order
    if order == 2:
        for t in t_points:
            x_points.append(x)
            k1 = h*f(x,t)
            k2 = h*f(x+0.5*k1, t+0.5*h)
            x += k2
    
    # solve range-kutta 4th order
    if order == 4:
        for t in t_points:
            x_points.append(x)
            k1 = h*f(x,t)
            k2 = h*f(x+0.5*k1,t+0.5*h)
            k3 = h*f(x+0.5*k2,t+0.5*h)
            k4 = h*f(x+k3,t+h)
            x += (k1+2*k2+2*k3+k4)/6
            
    plot(t_points,x_points)
    xlabel("t")
    ylabel("x(t)")
