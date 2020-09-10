
# Dictionary for Central Body
import numpy as np 

# Altitudes and Relative Densities
atm = np.array([
            [63.096,2.059e-4],
            [251.189,5.909e-11],
            [1000,3.561e-15]
            ])

# Constants For Earth
earth = {

    'name':'Earth',
    'mass':5.972e24,
    'mu':398600.4418,
    'radius':6378,
    'J2':1.08262668e-3,
    'zs' : atm[:,0],
    'rhos':atm[:,1]*10**8,
    'Earth_rot': np.array([0.0,0.0,72.9211e-6])

}