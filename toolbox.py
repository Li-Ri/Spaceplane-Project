import numpy as np 
import dictionary as dy
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

d2r = math.pi/180
cb = dy.earth

def plot_n_orbit(rs,labels,cb=dy.earth, show_plot=False,save_plot=False, Title= 'Many Orbits',colour = ['b','.-r','g'],linewidth=[1.5,.5,1.5]):

        fig = plt.figure(figsize=(12,8))
        ax = fig.add_subplot(111,projection='3d')
        #plot trjacectory
        legend = [
            'Final Orbit',
            'Trajectory',
            'Final Orbit'
        ]
        i=0
        
        for r in rs:
            ax.plot(r[:,0],r[:,1],r[:,2],colour[i], label=labels[i],linewidth=linewidth[i])
            ax.legend(legend)
            ax.set_xlabel('x-axis (km)')
            ax.set_ylabel('y-axis (km)')
            ax.set_zlabel('z-axis (km)')
            i+=1
            
        

        if show_plot:
            plt.show()

def COE2RV(mu,a, ecc, i, RA, w, theta):

    cos = math.cos 
    sin = math.sin
    sqrt = math.sqrt
    rad = math.radians
    p = a*(1-ecc**2)
    # Position vector components in x,y,z

    r_x = (p * cos(rad(theta))) / (1 + ecc * cos(rad(theta)))
    r_y = (p * sin(rad(theta))) / (1 + ecc * cos(rad(theta)))
    r_z = 0

    # Velocity vector components in x,y,z

    v_x = -sqrt(mu / p) * sin(rad(theta))
    v_y = sqrt(mu / p) * (ecc + cos(rad(theta)))
    v_z = 0

    # Rotation Matrix

    phi_11 = cos(rad(w))*cos(rad(RA)) - sin(rad(w))*cos(rad(i))*sin(rad(RA))
    phi_12 = -sin(rad(w))*cos(rad(RA)) - cos(rad(w))*cos(rad(i))*sin(rad(RA))
    phi_13 = sin(rad(i))*sin(rad(RA))
    phi_21 = cos(rad(w))*sin(rad(RA)) + sin(rad(w))*cos(rad(i))*cos(rad(RA))
    phi_22 = -sin(rad(w))*sin(rad(RA)) + cos(rad(w))*cos(rad(i))*cos(rad(RA))
    phi_23 = -sin(rad(i))*cos(rad(RA))
    phi_31 = sin(rad(w))*sin(rad(i))
    phi_32 = cos(rad(w))*sin(rad(i))
    phi_33 = cos(rad(i))

    phi_matrix = np.array(

                        [[phi_11, phi_12, phi_13],

                        [phi_21, phi_22, phi_23],

                        [phi_31, phi_32, phi_33]]
                                
                        )

    # Position and Velocities transformed into a 1D array

    pos = np.array([
                    r_x,
                    r_y,
                    r_z
                    ])

    vel = np.array([
                    v_x,
                    v_y,
                    v_z
                    ])

    # 3D Transformed Positions and Velocities
    Position_ijk = np.matmul(phi_matrix,np.transpose(pos))
    Velocity_ijk = np.matmul(phi_matrix,np.transpose(vel))

    return Position_ijk, Velocity_ijk
      
def ecc_anom(e,M_anom):
    pi = math.pi
    sin = math.sin
    cos = math.cos

    margin = 1*10**(-8)
    ratio = 1

    if M_anom < pi:
        E = M_anom + e/2
    else:
        E = M_anom - e/2

    while abs(ratio) > margin:
        ratio = (E - e*sin(E) - M_anom)/(1-e*cos(E))
        E -= ratio    
    return E

def rel_rho(z,zs=cb['zs'],rhos=cb['rhos']):
    if not 1.0<z<1000:
        return [[0.0,0.0],[0.0,0.0]]


    for n in range(len(rhos)-1):

        if zs[n]<z<zs[n+1]:

            return [[rhos[n],rhos[n+1]],[zs[n],zs[n+1]]]

    return [[0.0,0.0],[0.0,0.0]]

def calc_rho(z):

    dens,zs = rel_rho(z)

    if dens[0] == 0:
        return 0
    
    h0 = -(zs[1]-zs[0])/math.log(dens[1]/dens[0]) 

    return dens[0]*math.exp(-(z-zs[0])/h0)





def RV2COE(r_ijk,v_ijk,mu):


    # Magnitudes of vectors r and v

    mag_r = np.linalg.norm(r_ijk)
    mag_v = np.linalg.norm(v_ijk)

    # Angular Momentum
    h = np.cross(r_ijk,v_ijk)

    mag_h = np.linalg.norm(h)

    '''
    K = [0,0,1]

    
    # Node Vector
    node_v = np.cross(K,h)

   
    # Node Vector Magnitude
    mag_node = np.linalg.norm(node_v)
    '''

    # eccentricity vector
    #ecc_v = (1 / mu) * ((((mag_v **2) - (mu/mag_r))*r_ijk) - (np.dot(r_ijk,v_ijk)*v_ijk))
    ecc_v = ((mag_v**2 -cb['mu']/mag_r)*r_ijk-np.dot(r_ijk,v_ijk)*v_ijk)/cb['mu']
    
    # eccentricity
    e = np.linalg.norm(ecc_v)

    # energy equation
    eps = ((mag_v**2) / 2) - (cb['mu']/mag_r)

    # Semi-major axis

    a = - (cb['mu']/(2 * eps))

    # Semi-Latus Rectum
    #p = mag_h**2/cb['mu']

    # Inclination
    i = math.acos(h[2]/mag_h)

    # Right Ascension of the Ascending Node
    '''
    RA = math.acos(node_v[0] / mag_node)

    if node_v[1] < 0 :

         RA = 2*math.pi - RA
 

         

    # Argument of Perigee

    w = math.acos(np.dot(node_v,ecc_v) / mag_node/e)

    if ecc_v[2] < 0:

        w = 2*math.pi - w
    '''
       

    # True Anomaly

    v_anom = math.acos(np.dot(ecc_v,r_ijk) / (e * mag_r))
    
    if np.dot(r_ijk,v_ijk) < 0:

        v_anom = 2*math.pi - math.acos(np.dot(ecc_v,r_ijk) / (e * mag_r))

    return [e,a,i,v_anom]


def flight_path_graph(ys):

    r = ys[0][:,:3]
    y = ys[0][:,3:6]

    gamma = np.zeros((len(r)-1,1))
   

    for i in range(len(r)-1):
    
        gamma[i,:] = math.asin((np.dot(r[i,:],y[i,:]))/(np.linalg.norm(r[i,:])*np.linalg.norm(y[i,:])))

    return gamma*180/math.pi

def control_graph(u,t):
    
    
    throttle = np.array(u[0:6])
    alpha = np.array(u[6:12])
    sample_time = [0,t/6, 2*t/6,3*t/6, 4*t/6, 5*t/6]


    nthrottle = interp1d(sample_time,throttle,bounds_error=False,fill_value='extrapolate')
    nalpha = interp1d(sample_time,alpha,bounds_error=False,fill_value='extrapolate')
    time_span = np.linspace(0,t,len(sample_time))

    plt.plot(time_span,nalpha(time_span),'g')
    plt.grid()
    plt.xlabel('Time (s)')
    plt.ylabel('Angle of Attack (rad)')
    plt.show()


    plt.plot(time_span,nthrottle(time_span), 'g')
    plt.grid()
    plt.xlabel('Time (s)')
    plt.ylabel('Throttle')
    plt.show()

u =[ 1.        ,  1.        ,  1.        ,  1.        ,  0.69976436,
        0.77309885,  1.        ,  1.        ,  1.        ,  1.        ,
       -1.        ,  1.        ]

t = 24*35*3600
control_graph(u,t)
