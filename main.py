
import numpy as np 
from scipy.integrate import ode
from scipy.integrate import odeint
import dictionary as dy 
import math
from scipy.optimize import minimize, Bounds
from scipy.interpolate import interp1d
import toolbox as tb
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from orbit_propagator import hohm
from orbit_propagator import propagator



def spacecraft():
    return {

        'thrust':29400,
        'Isp':245,
        'CD':0.1,
        'area':30*(1*10**(-3))**2
    }

# craft spec and central body
sp = spacecraft()
body = dy.earth

def TwoBody(t,y,u):
    ''' 
    'integrate_thrust.m' algorithm from Orbital Mechanics for Engineering Students
    '''
    rx,ry,rz,vx,vy,vz,m = y

    r = np.array([rx,ry,rz]) 
    v = np.array([vx,vy,vz])

    norm_r = np.linalg.norm(r)
    norm_v = np.linalg.norm(v)

    time = 3600*4


    throttle = np.array(u[0:6])
    alpha = np.array(u[6:12])
    sample_time = [0,time/6, 2*time/6,3*time/6, 4*time/6, 5*time/6]


    nthrottle = interp1d(sample_time,throttle,bounds_error=False,fill_value='extrapolate')
    nalpha = interp1d(sample_time,alpha,bounds_error=False,fill_value='extrapolate')
    

    # Two body relative motion
    a = -r*body['mu']/norm_r**3

    # Adding perturbations for J2, thrust and atmospheric drag
    z2 = r[2]**2
    r_2 = norm_r**2
    t_x = r[0]/norm_r*(5*z2/r_2-1) 
    t_y = r[1]/norm_r*(5*z2/r_2-1)
    t_z = r[2]/norm_r*(5*z2/r_2-3)

    # J2 Acceleration
    acc_j2 = 1.5*body['J2']*body['mu']*body['radius']**2/norm_r**4*np.array([t_x,t_y,t_z])
    a+=acc_j2

    
    gamma = math.asin(np.dot(r,v)/(norm_r*norm_v))


    T = nthrottle(t)*sp['thrust']*np.array([math.cos(gamma+nalpha(t)*180/math.pi),
                                            math.sin(gamma+nalpha(t)*180/math.pi),
                                            0])

    
    thrust_acc =(T/m)*(v/norm_v)/1000 #km/s^2
    a+=thrust_acc

    mdot = -np.linalg.norm(T)/9.87/sp['Isp']
    
    # atmospheric drag
    
    z = norm_r - body['radius']

    dens = tb.calc_rho(z)
    v_rel = v - np.cross(body['Earth_rot'],r)

    atmos_drag = -v_rel*0.5*dens*(np.linalg.norm(v_rel))*sp['CD']*sp['area']/m

    a+=atmos_drag
  
    func = [vx,vy,vz,a[0],a[1],a[2],mdot]
        
    return func

def solver(u,t,dt,z0):

    step = 0
    nsteps= int(np.ceil(t/dt))
   
    ts = np.zeros((nsteps+1,1))
    ys = np.zeros((nsteps+1,7))


    ts[0] = 0
    ys[0,:] = z0
    
    solver = ode(TwoBody)
    solver.set_integrator('dopri5',nsteps=t/dt)
    solver.set_initial_value(ys[0,:],0)
    solver.set_f_params(u)

    while solver.successful() and step < nsteps:
        solver.integrate(solver.t+dt)

        ys[step] = solver.y
        ts[step] = solver.t

        step+=1

    

    # time
    ts = ts[:step]
    # masses
    ms = ys[:step,-1]
 
    return ys,ms,ts
   
def objective(u,t,dt,y):
        
    sol = solver(u,t,dt,y)
    m = y[-1]
    mf = sol[1][-1]
    dv = sp['Isp']*9.81*np.log(m/mf)

    return dv

def constraint(u,t,dt,z0,keptarget):

    sol = solver(u,t,dt,z0)
    
    a,inc,e = keptarget

    rs = sol[0][:,:3]
    vs = sol[0][:,3:6]
    '''
    kep_elem = []

    kep_elem.append(tb.RV2COE(rs, vs, mu=body['mu']))
    '''
    
    kep_elem = np.zeros((len(rs)-1,4))

    for i in range(len(rs)-1):
    
        kep_elem[i,:] = tb.RV2COE(rs[i,:], vs[i,:], mu=body['mu'])
    
    wa = 1000
    wi = 0
    we = 0.1

    diffa = abs(kep_elem[-1][1]-a)/(wa*90)

    diffe = abs(kep_elem[-1][0]-e)/(we*10)


    if inc==0:

        diffi = 0
    else: 
        diffa = abs(kep_elem[-1][1]-a)/(wa*33.3)

        diffe = abs(kep_elem[-1][0]-e)/(we*33.3)

        diffi = abs(kep_elem[-1][2]-inc)/(wi*33.3)

    difftotal = diffa + diffe + diffi
    
    print(kep_elem[-1][1],kep_elem[-1][0],difftotal)
    return difftotal

def constraint2(u,t,dt,z0,kep_target):

    #sol=solver(u,t,dt,z0)

    m_dry = 2000
    m_final = sol[-1]

    return m_final- m_dry
def calc_elem(y):

    rs = y[0][:,:3]
    vs = y[0][:,3:6]

    kep_elem = np.zeros((len(rs)-1,4))

    for i in range(len(rs)-1):
    
        kep_elem[i,:] = tb.RV2COE(rs[i,:], vs[i,:], mu=body['mu'])

    return kep_elem

def plot_kep_elements(ts,kep_elem,figsize=(16,8),title="Keplerian Elements"):

    # Plotting the keplerian elements over time using the 'calculate_elements' function
    #Create multiple Axes and title
    fig,ax = plt.subplots(nrows = 2, ncols=2,figsize=figsize)
    fig.suptitle(title,fontsize=20)
    sol = solver(u,t,dt,y)
    fp_angle = tb.flight_path_graph(sol)
    xlabel='time (s)'
        

    # Eccentricity
    ax[0,0].plot(ts,kep_elem[:,0])
    #ax[0,0].set_title('Eccentricity / Time')
    ax[0,0].grid(True)
    ax[0,0].set_ylabel('e')
    ax[0,0].set_xlabel(xlabel)

        # Semi-Major Axis
    ax[0,1].plot(ts,kep_elem[:,1])
    #ax[0,1].set_title('Semi-Major Axis / Time')
    ax[0,1].grid(True)
    ax[0,1].set_ylabel('Semi-Major Axis (km)')
    ax[0,1].set_xlabel(xlabel)
    
    '''  # Inclination
    ax[1,0].plot(ts,kep_elem[:,2])
    ax[1,0].grid(True)
    ax[1,0].set_ylabel('IFuel')
    ax[1,0].set_xlabel(xlabel)
    '''
        # Flight path angle
    ax[1,0].plot(ts,fp_angle[:,0],color='r')
    ax[1,0].grid(True)
    ax[1,0].set_ylabel('Flight Path Angle (deg)')
    ax[1,0].set_xlabel(xlabel)
        
    
    ax[1,1].plot(ts,sol[1],color='k')
    ax[1,1].grid(True)
    ax[1,1].set_ylabel('Fuel Consumption (kg)')
    ax[1,1].set_xlabel(xlabel)
    
    # Calculating Perigee and Apogee radius
    
    rp = (kep_elem[:,1]*(1-kep_elem[:,0]**2))/(1+kep_elem[:,0]*math.cos(0))
    ra = (kep_elem[:,1]*(1-kep_elem[:,0]**2))/(1+kep_elem[:,0]*math.cos(math.pi))
    
    fig2, ax2 = plt.subplots(nrows = 2, ncols=2,figsize=figsize)
    fig2.suptitle(title,fontsize=20)

    ax2[0,0].plot(ts,rp,label='Periapsis')
    ax2[0,0].plot(ts,ra, label='Apoapsis')
    ax2[0,0].grid(True)
    ax2[0,0].set_xlabel('Time (s)')
    ax2[0,0].set_ylabel('Radius (km)')
    ax2[0,0].legend()


    plt.show()
    
if __name__ == '__main__':

# target elements [semi-major axis,inclination,eccentricity]
    kep_target=[9000,0,0]

# circular orbit
    r = 6371+200
    v = math.sqrt(body['mu']/r)
    a_trans = (r+kep_target[0])/2


    rf,vf = tb.COE2RV(body['mu'],a_trans,kep_target[2],0,0,0,0)

    y = [r,0,0,0,4,0,30000]
    yf = [rf[0],0,0,0,vf[1],0]
    t = 3600*4
    dt = 100


# constraint
    cons = ({ 'type':'eq','fun': constraint,'args':(t,dt,y,kep_target)})
        #{'type':'ineq','fun': constraint2,'args':(t,dt,y,kep_target)}

# control
    u = np.array(
       [0.001, .001, 0.001, 0.001,
       0.001, 0.001, 0.1, .1,
       .1, 9.76128628e-01, 1.79814754e-01, 2.06618409e-01])




# minimize function

    res = minimize(
                objective,u,
                bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1), # throttle bounds
                        (-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1), ), # alpha bounds
                args=(t,dt,y),
                constraints=cons,
                options={'maxiter':100},
                tol=0.001,
                method='SLSQP'
                
                )
 
    print(res)

# reassigning new control parameters
    u = res.x

    sol=solver(u,t,dt,y)



# minimized fuel consumption
#min_fuel = objective(u,t,dt,y)
#print(min_fuel)
# create figure
    rs = sol[0]
    rs = rs[0:len(sol[0])-1]

# plot trjacectory and elements

# plotting initial and final orbits
    op = propagator(y[0:6],2*math.pi*math.sqrt((r**3/body['mu'])),dt)

# stop propagation and plot final orbit
    u1 = (0,0,0,0,0,0,
    0,0,0,0,0,0)
    sol2 = solver(u1,t,dt,rs[-1])
    rs2 = sol2[0]
    rs2 = rs2[0:len(sol2[0])-1]

    tb.plot_n_orbit([op.rs,rs,rs2],labels=['0','1','2'])

# plot kep elements
    kep_elem = calc_elem(sol)
    plot_kep_elements(sol[2],kep_elem)

    tb.control_graph(u,t)









