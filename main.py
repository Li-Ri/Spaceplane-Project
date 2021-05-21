
import numpy as np 
from scipy.integrate import ode
from scipy.integrate import odeint
import central_body_constants as dy
import math
from scipy.optimize import minimize, Bounds
from scipy.interpolate import interp1d
import toolbox as tb
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Orbit Class

class Orbit_Optimizer:

    def __init__(self,u,t,dt,y,body=dy.earth):

        self.u = u
        self.t = t
        self.dt = dt
        self.y = y    
        self.body =body

        self.sol = self.solver(self.u,self.t,self.dt,self.y)
    @staticmethod
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

        # Flight Path Angle (0 if circular orbit)
        gamma = math.asin(np.dot(r,v)/(norm_r*norm_v))


        # Thrust with parameter gamma and nalpha(t) controlling angle of attack and nthrottle managing power

        T = nthrottle(t)*sp['thrust']*np.array([math.cos(gamma+nalpha(t)*180/math.pi),
                                                math.sin(gamma+nalpha(t)*180/math.pi),
                                                0])

        
        thrust_acc =(T/m)*(v/norm_v)/1000 #km/s^2
        a+=thrust_acc

        # Mass Flow Rate
        mdot = -np.linalg.norm(T)/9.87/sp['Isp']
        
        # atmospheric drag
        
        z = norm_r - body['radius']

        dens = tb.calc_rho(z)
        v_rel = v - np.cross(body['Earth_rot'],r)

        atmos_drag = -v_rel*0.5*dens*(np.linalg.norm(v_rel))*sp['CD']*sp['area']/m

        a+=atmos_drag
    
        func = [vx,vy,vz,a[0],a[1],a[2],mdot]
            
        return func
 
    def solver(self,u,t,dt,z0):

        self.step = 0
        self.nsteps= int(np.ceil(self.t/self.dt))
    
        self.ts = np.zeros((self.nsteps+1,1))
        self.ys = np.zeros((self.nsteps+1,7))


        self.ts[0] = 0
        self.ys[0,:] = z0
        
        self.solve = ode(self.TwoBody)
        self.solve.set_integrator('dopri5',nsteps=self.t/self.dt)
        self.solve.set_initial_value(self.ys[0,:],0)
        self.solve.set_f_params(self.u)

        while self.solve.successful() and self.step < self.nsteps:
            self.solve.integrate(self.solve.t+self.dt)

            self.ys[self.step] = self.solve.y
            self.ts[self.step] = self.solve.t

            self.step+=1

        

        # time
        self.ts = self.ts[:self.step]
        # masses
        self.ms = self.ys[:self.step,-1]
    
        return self.ys,self.ms,self.ts
    
    def objective(self,u,t,dt,z0):

        self.sol = self.solver(self.u,self.t,self.dt,self.y)
            
        m = self.y[-1]
        mf = self.sol[1][-1]
        dv = sp['Isp']*9.81*np.log(m/mf)

        return dv

    def constraint(self,u,t,dt,z0,keptarget):
        self.a,self.inc,self.e = keptarget

        self.sol = self.solver(self.u,self.t,self.dt,self.y)

        self.rs = self.sol[0][:,:3]
        self.vs = self.sol[0][:,3:6]
        '''
        kep_elem = []

        kep_elem.append(tb.RV2COE(rs, vs, mu=body['mu']))
        '''
        
        self.kep_elem = np.zeros((len(self.rs)-1,4))

        for i in range(len(self.rs)-1):
        
            self.kep_elem[i,:] = tb.RV2COE(self.rs[i,:], self.vs[i,:], mu=body['mu'])
        
        self.wa = 1000
        self.wi = 0
        self.we = 0.1

        diffa = abs(self.kep_elem[-1][1]-self.a)/(self.wa*90)

        diffe = abs(self.kep_elem[-1][0]-self.e)/(self.we*10)


        if self.inc==0:

            diffi = 0
        else: 
            diffa = abs(self.kep_elem[-1][1]-self.a)/(self.wa*33.3)

            diffe = abs(self.kep_elem[-1][0]-self.e)/(self.we*33.3)

            diffi = abs(self.kep_elem[-1][2]-self.inc)/(self.wi*33.3)

        difftotal = diffa + diffe + diffi
        
        # Show Elements
        print(self.kep_elem[-1][1],self.kep_elem[-1][0],difftotal)

        return difftotal

    def calc_elem(self):
  
        self.rs = self.sol[0][:,:3]
        self.vs = self.sol[0][:,3:6]

        self.kep_elem = np.zeros((len(self.rs)-1,4))

        for i in range(len(self.rs)-1):
        
            self.kep_elem[i,:] = tb.RV2COE(self.rs[i,:], self.vs[i,:], mu=body['mu'])

        return self.kep_elem, self.sol[2]

    def plot_kep_elements(self,ts,kep_elem,figsize=(16,8),title="Keplerian Elements"):

        # Plotting the keplerian elements over time using the 'calculate_elements' function
        #Create multiple Axes and title
        fig,ax = plt.subplots(nrows = 2, ncols=2,figsize=figsize)
        fig.suptitle(title,fontsize=20)
        self.sol = self.solver(self.u,self.t,self.dt,self.y)
        fp_angle = tb.flight_path_graph(self.sol)
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
            
        
        ax[1,1].plot(ts,self.sol[1],color='k')
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
        


# Space Craft and Main Function
def spacecraft():
    return {

        'thrust':30000,
        'Isp':245,
        'CD':0.1,
        'area':30*(1*10**(-3))**2
    }
def main():

    # target elements [semi-major axis,inclination,eccentricity]
    kep_target=[9000,0,0]

    # circular orbit
    r = 6371+200
    v = math.sqrt(body['mu']/r)
    a_trans = (r+kep_target[0])/2


    rf,vf = tb.COE2RV(body['mu'],a_trans,kep_target[2],0,0,0,0)

    y = [r,0,0,0,4,0,300]
    yf = [rf[0],0,0,0,vf[1],0]
    t = 3600*4
    dt = 100

    u = np.array(
    [0.00, .00, 0.00, 0.00,
    0.001, 0.001, 0.1, .1,
    .1, 9.76128628e-01, 1.79814754e-01, 2.06618409e-01])

    propagator = Orbit_Optimizer(u,t,dt,y)
    # constraint
    cons = ({ 'type':'eq','fun': propagator.constraint, 'args':(t,dt,y,kep_target)})
        #{'type':'ineq','fun': constraint2,'args':(t,dt,y,kep_target)}


    # # minimize function

    res = minimize(
                propagator.objective,u,
                bounds=((0,1),(0,1),(0,1),(0,1),(0,1),(0,1), # throttle bounds
                        (-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1), ), # alpha bounds
                args=(t,dt,y),
                constraints=cons,
                options={'maxiter':100},
                tol=0.001,
                method='SLSQP'
                )
    
    return res.x





if __name__ == '__main__':
    body = dy.earth
    kep_target=[9000,0,0]

    # circular orbit
    r = 6371+200
    v = math.sqrt(body['mu']/r)
    a_trans = (r+kep_target[0])/2


    rf,vf = tb.COE2RV(body['mu'],a_trans,kep_target[2],0,0,0,0)

    y = [r,0,0,0,4,0,300]
    yf = [rf[0],0,0,0,vf[1],0]
    t = 3600*4
    dt = 100

    u = np.array(
    [0.00, .00, 0.00, 0.00,
    0.001, 0.001, 0.1, .1,
    .1, 9.76128628e-01, 1.79814754e-01, 2.06618409e-01])
    # Initiat Spacecraft
    sp = spacecraft()
    #Initiate Earths Gravity and J2
    

    propagator = Orbit_Optimizer(u,t,dt,y)
    kep, ts = propagator.calc_elem()
    propagator.plot_kep_elements(ts,kep)


# reassigning new control parameters
    
    

    












