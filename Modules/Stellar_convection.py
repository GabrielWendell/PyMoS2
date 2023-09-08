"""
This is the Stellar_model module from the PyMoS4 package whose objective 
is to simulate hydrodynamically the convection in the stellar interior

Rocha, G. W. C. (2021)
"""

from matplotlib import pyplot as plt
from matplotlib import animation
import FVis3 as FVis 
import numpy as np
# import smplotlib



# ---------------------------------------------------------------------------------------------------------------------
def opt_plot():
    """
    This function helps to improve the aesthetics of a single graphic!
    """
    # plt.style.use('dark_background')
    plt.grid(True, linestyle = ':', color = '0.50')
    plt.minorticks_on()
    plt.tick_params(axis = 'both', which = 'minor', direction = "in",
                    top = True, right = True, length = 5, width = 1,labelsize = 15)
    plt.tick_params(axis='both', which='major', direction = "in",
                    top = True, right = True, length = 8, width = 1, labelsize = 15)

def opt_plot_mod(ax):
    """
    This function helps to improve the aesthetics of a graph tied to a certain axis!
    """
    # plt.style.use('dark_background')
    ax.grid(True, linestyle = ':', color = '0.50')
    ax.minorticks_on()
    ax.tick_params(axis = 'both', which = 'minor', direction = "in",
                        top = True,right = True, length = 5,width = 1, labelsize = 15)
    ax.tick_params(axis = 'both', which = 'major', direction = "in",
                        top = True, right = True, length = 8,width = 1, labelsize = 15) 
    
def style():
    """
    This function applies this style to all graphics generated by the code.
    """
    from matplotlib import rcParams
    import matplotlib.font_manager as font_manager
    
    # Uncomment the lines below to use custom fonts by specifying the path to them.
    
    # font_dirs = 'C:/Python36/Lib/site-packages/matplotlib/mpl-data/fonts/ttf'
    # font_files = font_manager.findSystemFonts(fontpaths = font_dirs)
    # font_list = font_manager.createFontList(font_files)
    # font_manager.fontManager.ttflist.extend(font_list)
    rcParams['font.family'] = 'monospace'
    rcParams['font.size'] = 18
    rcParams['mathtext.fontset'] = 'dejavuserif'
    rcParams['axes.linewidth'] = 1.2

style()
    
    


# ---------------------------------------------------------------------------------------------------------------------



# 

class twoD_convection(object):

   # Set variables
    def __init__(self): 

        # Constants and parameters used
        self.G = 6.672e-11                            # [G] = N * m^2 * kg^-2
        self.M_Sun = 1.989e30                         # [M_Sun] = kg
        self.R_Sun = 6.96e8                           # [R_Sun] = m
        self.g = self.G*self.M_Sun/(self.R_Sun**2)    # [g] = m * s^(-2)
        self.k_B = 1.382e-23                          # [k_B] = m^2 * kg * s^-2 * K^-1
        self.m_u = 1.6605e-27                         # [m_u] = kg
        self.mu = 0.61                                # Mean Molecular Mass
        self.gamma = 5./3.                            # Adiabatic index (monoatomic ideal gas)
        
        self.nx = 300                                 # Number of cells in horizontal direction
        self.ny = 100                                 # Number of cells in vertical direction
        self.X = 12.e6                                # [X] = m (Width of the box)
        self.Y = 4.e6                                 # [Y] = m (Height of the box)
        self.Delta_x = self.X/self.nx                 # Size of a cell, x-direction
        self.Delta_y = self.Y/self.ny                 # Size of a cell, y-direction

        # Arrays with values inside the box - these will be animated, 
        # first index specifies width and the second index - height
        self.rho_box = np.zeros((self.nx,self.ny))    # [rho_box] = kg * m^(-3) : Density
        self.u_box = np.zeros((self.nx,self.ny))      # [u_box] = m * s^(-1) : Horizontal velocity
        self.w_box = np.zeros((self.nx,self.ny))      # [w_box] = m * s^(-1) : Vertical velocity
        self.e_box = np.zeros((self.nx,self.ny))      # [e_box] = J * m^(-3) : Internal energy density
        self.P_box = np.zeros((self.nx,self.ny))      # [P_box] = Pa : Pressure
        self.T_box = np.zeros((self.nx,self.ny))      # [T_box] = K : Temperature
        
        # Store all time steps (used for flow animation)
        self.dt_array = [] 
        # Store convective flux along y-direction for each time-step
        self.flux_array = [] 
        
    def initialize(self, perturbation = True):
       #initialize T, P, e, rho and turn off the perturbation if needed
        T_top = 5778.                                 # [T_top] = K : Temperature at the top of the box (Solar Photosphere)
        P_top = 1.8e8                                 # [P_top] = Pa : Pressure at the top of the box (Solar Photosphere)
        Gradient = 2./5. + 0.001                      # Double logarithmic gradient
        # Initialize temperature, dependent on height
        for j in range(self.ny): 
           self.T_box[:,j] = T_top + (Gradient*self.mu*self.m_u*(self.g)/self.k_B)*self.Delta_y*(self.ny-j-1)

        # Initialize P, dependent on the unperturbed temperature
        self.P_box = P_top*((self.T_box/T_top)**(1./Gradient)) 

        # Add perturbation to the temperature if needed
        if perturbation == True: 
           x_array = np.linspace( 0.,self.X-self.Delta_x, self.nx)
           y_array =  np.linspace( 0.,self.Y-self.Delta_y, self.ny)
           sigma = 1.*1e6                             # Standard Deviation
           x,y = np.meshgrid(x_array, y_array)
           # Transpose this matrix to match the shape of T-array
           gaussian_perturbation = (6000.*np.exp(-((x-self.X/2.)**2 + (y-self.Y/2.)**2 )/(2.*sigma**2))).T 
           self.T_box = self.T_box + gaussian_perturbation
        
        # Initialize rho and e 
        self.rho_box = (self.mu*self.m_u*self.P_box)/(self.k_B*self.T_box)
        self.e_box = (self.rho_box*self.T_box)*self.k_B/((self.gamma -1.)*self.mu*self.m_u)
        '''
        A different way to initialize e, and which produces almost the same results
        '''
        # self.e_box = self.P_box/(self.gamma-1.) 
        
    def timestep(self, rho_time_derivative, e_time_derivative, u, w): #calculate timestep
        # Exclude stationary points and very small values from velocity fields
        u_without_singularities = []
        w_without_singularities = []
        for i in range(self.nx-2):
           for j in range(self.ny-2):
              if abs(u[i,j]) > 1e-5:
                 u_without_singularities.append(u[i,j])
              if abs(w[i,j]) > 1e-5:
                 w_without_singularities.append(w[i,j])
        u_without_singularities = np.array(u_without_singularities)
        w_without_singularities = np.array(w_without_singularities)
        
        if len(u_without_singularities) == 0:
           max_relative_x = 0.
        else:
           relative_x = abs(u_without_singularities/self.Delta_x)
           # Largest relative change on the grid for x-position 
           max_relative_x = np.amax(relative_x) 
        if len(w_without_singularities) == 0:
           max_relative_y = 0.
        else:
           relative_y = abs(w_without_singularities/self.Delta_y)
           # Largest relative change on the grid for y-position
           max_relative_y = np.amax(relative_y) 
        
        relative_rho = abs(rho_time_derivative*(1./self.rho_box[1:-1,1:-1]))
        # Largest relative change on the grid for rho
        max_relative_rho = np.amax(relative_rho)
        relative_e = abs(e_time_derivative*(1./self.e_box[1:-1,1:-1]))
        # Largest relative change on the grid for e
        max_relative_e = np.amax(relative_e) 
        
        # Largest relative change per time-step
        delta = np.max([max_relative_rho, max_relative_e, max_relative_x, max_relative_y]) 
        # Avoid very big time-steps
        if delta < 1.: 
           dt = 0.1
        else:
           # Some small number    
           p = 0.1 
           dt = p/delta
        
        return dt    

    def boundary_conditions(self):
       #apply boundary conditions for e, rho, u and w
        
        # Horizontal boundaries
        self.e_box[0,:] = np.copy(self.e_box[-2,:])
        self.e_box[-1,:] = np.copy(self.e_box[1,:])
        self.rho_box[0,:] = np.copy(self.rho_box[-2,:])
        self.rho_box[-1,:] = np.copy(self.rho_box[1,:])
        self.u_box[0,:] = np.copy(self.u_box[-2,:])
        self.u_box[-1,:] = np.copy(self.u_box[1,:])
        self.w_box[0,:] = np.copy(self.w_box[-2,:])
        self.w_box[-1,:] = np.copy(self.w_box[1,:])
       
        # Vertical boundaries
        self.w_box[:,0] = 0.
        self.w_box[:,-1] = 0.
        self.u_box[:,0] = (4.*self.u_box[:,1] - self.u_box[:,2])/3.
        self.u_box[:,-1] = (4.*self.u_box[:,-2] - self.u_box[:,-3])/3.
        e_T_1 = (2.*self.Delta_y*self.mu*self.m_u*(-self.g))/(self.k_B*self.T_box[:,0])
        self.e_box[:,0] = (4.*self.e_box[:,1] - self.e_box[:,2])/(e_T_1 + 3.)
        e_T_2 = (2.*self.Delta_y*self.mu*self.m_u*(-self.g))/(self.k_B*self.T_box[:,-1])
        self.e_box[:,-1] = (self.e_box[:,-3] -4.*self.e_box[:,-2])/(e_T_2 - 3.)
        self.rho_box[:,0] = ((self.gamma-1.)*self.mu*self.m_u*self.e_box[:,0])/(self.k_B*self.T_box[:,0])
        self.rho_box[:,-1] = ((self.gamma-1.)*self.mu*self.m_u*self.e_box[:,-1])/(self.k_B*self.T_box[:,-1])
        
        
    def central_x(self,func): #central difference scheme in x-direction
        derivative = (np.roll(func,-1, axis=0) - np.roll(func, 1, axis=0))/(2.*self.Delta_x)
        return derivative[1:-1,1:-1]

    def central_y(self,func): #central difference scheme in y-direction
        derivative = (np.roll(func,-1, axis=1) - np.roll(func, 1, axis=1))/(2.*self.Delta_y)
        return derivative[1:-1,1:-1]
   
    def upwind_x(self,func,wind): #upwind difference scheme in x-direction
        derivative = np.zeros((self.nx,self.ny))
        derivative[np.where(wind >= 0.)] = ((func - np.roll(func,1, axis=0))[np.where(wind >= 0.)])/self.Delta_x
        derivative[np.where(wind < 0.)] = ((np.roll(func,-1, axis=0) - func)[np.where(wind < 0.)])/self.Delta_x
        return derivative[1:-1,1:-1]

    def upwind_y(self,func,wind): #upwind difference scheme in y-direction
        derivative = np.zeros((self.nx,self.ny))
        derivative[np.where(wind >= 0.)] = ((func - np.roll(func,1, axis=1))[np.where(wind >= 0.)])/self.Delta_y
        derivative[np.where(wind < 0.)] = ((np.roll(func,-1, axis=1) - func)[np.where(wind < 0.)])/self.Delta_y
        return derivative[1:-1,1:-1]
        
        
           
    def hydro_solver(self): #hydrodynamic equations solver
        
        # Variable-arrays before update
        rho, u, w, e, P = self.rho_box[:], self.u_box[:], self.w_box[:], self.e_box[:], self.P_box[:]
        
        # Calculate time-derivatives inside the box, for elements with indices 1 to -2
        drho_dt = -self.rho_box[1:-1,1:-1]*(self.central_x(u) + self.central_y(w)) -self.u_box[1:-1,1:-1]*self.upwind_x(rho,u) - self.w_box[1:-1,1:-1]*self.upwind_y(rho,w)
        drhou_dt = -(self.rho_box*self.u_box)[1:-1,1:-1]*(self.upwind_x(u,u) + self.upwind_y(w,u))-self.u_box[1:-1,1:-1]*self.upwind_x((rho*u),u)-self.w_box[1:-1,1:-1]*self.upwind_y((rho*u),w)-self.central_x(P) 
        # Ideal equilibrium (for 'w') requires 'np.round' before hydro_equilibrium-factor to avoid numerical uncertainites
        hydro_eq = (-self.central_y(P)-self.rho_box[1:-1,1:-1]*self.g) 
        drhow_dt = -(self.rho_box*self.w_box)[1:-1,1:-1]*(self.upwind_x(u,w) + self.upwind_y(w,w))-self.w_box[1:-1,1:-1]*self.upwind_y((rho*w),w)-self.u_box[1:-1,1:-1]*self.upwind_x((rho*w),u) + hydro_eq
        de_dt = -(self.e_box + self.P_box)[1:-1,1:-1]*(self.central_x(u) + self.central_y(w)) - self.u_box[1:-1,1:-1]*self.upwind_x(e,u) - self.w_box[1:-1,1:-1]*self.upwind_y(e,w)
        
        # Find suitable time-step
        dt = self.timestep(drho_dt, de_dt, self.u_box[1:-1,1:-1], self.w_box[1:-1,1:-1]) 
        # print('Present time-step:', dt, 's')
        
        # Update elements with indices 1 to -2:
        self.rho_box[1:-1,1:-1] = rho[1:-1,1:-1] + drho_dt*dt
        self.u_box[1:-1,1:-1] = ((u*rho)[1:-1,1:-1] + drhou_dt*dt)/self.rho_box[1:-1,1:-1]   # Here division by rho already progressed in time
        self.w_box[1:-1,1:-1] = ((rho*w)[1:-1,1:-1] + drhow_dt*dt)/self.rho_box[1:-1,1:-1]   # Here division by rho already progressed in time
        self.e_box[1:-1,1:-1] = e[1:-1,1:-1] + de_dt*dt

        # Update the bundary points
        self.boundary_conditions() 

        # Update whole arrays of secondary variables
        self.P_box[:] = (self.gamma - 1.)*self.e_box
        self.T_box[:] = ((self.gamma - 1.)*self.mu*self.m_u*self.e_box)/(self.rho_box*self.k_B)
        
        # Store flux and time-steps needed for animation of the convective flux
        # Stores the total convective energy flux along y-direction
        ev_sum_now = []             # [ev_sum_now] = W * m^(-2)
        for j in range(self.ny):
           velocity_y = np.sqrt((self.u_box*self.u_box)[:,j] + (self.w_box*self.w_box)[:,j])
           A = np.sum(self.e_box[:,j]*velocity_y)
           # Adds elemetn at the end of the existing list
           ev_sum_now.append(A) 
        ev_sum_now = np.array(ev_sum_now)
        self.dt_array.append(dt)
        self.flux_array.append(ev_sum_now)
        
        # Return the time-step
        return dt 

    def animate_convective_flux(self, save = False): #make an animation of convective flux
        time_steps = np.array(self.dt_array)
        # Each sub-array corresponds to a given time step in time_steps
        # Each element in sub-array corresponds to different depth (y-level)
        convective_flux =  np.array(self.flux_array)  
    
        fig, ax = plt.subplots()
        # Defining 'x-axis', which is the height
        x_axis = np.linspace( 0., self.Y-self.Delta_y, self.ny)/1.e6  
        line, = ax.plot(x_axis, convective_flux[-3,:], color = 'blue')        
        
        # Text to display the current frame
        time_text = ax.text(0.05, 0.95, '', horizontalalignment = 'left',verticalalignment = 'top', transform = ax.transAxes) 

        # Initialize variables 
        def init():
           line.set_ydata(np.ma.array(x_axis, mask = True))
           plt.title("Convective flux as a function of Height", fontsize = 15)
           plt.xlabel("$Y$ $[Mm]$", fontsize = 18) 
           plt.ylabel("$F_C$ $[W/m^2]$", fontsize = 18)
           opt_plot()
           time_text.set_text('')
           
           # Return the variables that will be updated in each frame
           return line ,time_text      
    
        def frame(i):
           line.set_ydata(convective_flux[i])
           time_text.set_text(' Time = %.2f s' % np.sum(time_steps[0:i]))  
           return line , time_text

        anim = animation.FuncAnimation(fig, frame, len(time_steps), init_func = init, interval = 100, blit = True)
        if save == True:
           anim.save('Flux_movie.mp4') 
        plt.show()

if __name__ == '__main__':
    # Create an instance of the visualization class
    vis = FVis.FluidVisualiser() 
    # Here set 'True' if you want to run the sanity check 
    # Set 'False' if you want to run the code with convection
    sanity_check = False 
    if sanity_check == True:
       # ***Perform sanity check***
       
       # Create an instance of the 2D convection simulation
       sanity = twoD_convection() 
       # Turn-off the T-perturbation and show the hydrostatic equilibrium
       sanity.initialize(perturbation = False)
       vis.save_data(60., sanity.hydro_solver, rho = sanity.rho_box, u = sanity.u_box, w = sanity.w_box, e = sanity.e_box, P = sanity.P_box, T = sanity.T_box, sim_fps = 1)
       vis.animate_2D('T', extent = [0.,sanity.X*1e-6,0., sanity.Y*1e-6], matrixLike = False, save = False, units = {'Lx':'Mm', 'Lz':'Mm'}, folder = 'Sanity_check_data', video_fps = 3., anim_time = 60., anim_fps = 1)
    else:
       # Create an instance of the 2D convection simulation
       solver = twoD_convection() 
       solver.initialize(perturbation = True)
       # Save data
       vis.save_data(250., solver.hydro_solver, rho = solver.rho_box, u = solver.u_box, w = solver.w_box, e = solver.e_box, P = solver.P_box, T = solver.T_box, sim_fps = 1) 
       # Animate from 'Results_data' folder
       vis.animate_2D('T', extent = [0.,solver.X*1e-6,0.,solver.Y*1e-6], matrixLike = False, save = False, units = {'Lx':'Mm', 'Lz':'Mm'}, folder = 'Results_data') 
       # Make an animation of evolution of the convective flux (as a function of height, y)
       solver.animate_convective_flux(save = False) 