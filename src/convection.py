# visualiser
import FVis3 as FVis

import smplotlib
import numpy as np
import matplotlib.pyplot as plt

class _2Dconvection:

    def __init__(self,
        do_pert = False,        # True : perturbation, False : no perturbation
        A = 17000,              # [K] Amplitude of perturbation
        sigma = 0.6e6,          # [m] Standard deviation of perturbation
        x0 = 6e6,               # [m] Horizontal position of perturbation
        y0 = 2e6                # [m] Vertical position of perturbation
        ):

        """
        Defining constants and initializing simulation
        """

        ### Parameters:

        self.do_pert = do_pert

        self.A = A
        self.sigma = sigma
        self.x0 = x0
        self.y0 = y0

        ### Constants:

        self.mu = 0.61          # [1] Mean molecular weight
        self.m_u = 1.6605e-27   # [kg] Atomic mass unit
        self.kB = 1.3806e-23    # [J/K] Boltzmann's constant
        self.G = 6.6742e-11     # [Nm^2kg^(-2)] Gravitational constant

        self.gamma = 5/3        # [1] Ratio of specific heats for ideal gas


        ### General solar parameters:

        self.R0 = 6.96e8    # [m] Solar radius
        self.M0 = 1.989e30  # [kg] Solar mass

        ### Photosphere:

        self.P0 = 1.8e4     # [N/m^2] Pressure
        self.T0 = 5778      # [K] Temperature

        self.g = self.G*self.M0/self.R0**2  # [m/s^2] Grav. acceleration
        self.nabla = 2/5 + 1e-3             # [1] Temp. gradient

        ### Simulation

        self.p = 0.1

        ### Box:

        self.nx = 300           # Simulation points in x-direction
        self.ny = 100           # Simulation points in y-direction

        self.x_size = 12e6      # [m] Width of box
        self.y_size = 4e6       # [m] Height of box

        self.dx = self.x_size/self.nx    # [m] Cell size in x-direction
        self.dy = self.y_size/self.ny    # [m] Cell size in y-direction

        ### Initialising:

        self.initialise()

    def initialise(self):

        """
        Initialises simulation box as well as all the variables such as
        velocity, temperature, pressure, density and internal energy
        """

        ### Initialising box:

        x = np.linspace(0,self.x_size,self.nx)
        y = np.linspace(self.y_size,0,self.ny)

        self.X,self.Y = np.meshgrid(x,y)

        # From this we get a radius matrix R so that R[0,:] is the bottom
        # of the box and R[-1,:] is the top
        self.R = self.R0 - self.Y

        # Initialising temperature
        self.T = self.T0 \
        - self.mu*self.m_u*self.g*self.nabla/self.kB*(self.R - self.R0)

        # Initialising pressure
        self.P = self.P0*(self.T/self.T0)**(1/self.nabla)

        # Redefining temperature, with perturbation
        self.perturbation()

        if self.do_pert : self.T += self.G

        # Initialising density
        self.rho = self.mu*self.m_u/(self.kB*self.T)*self.P

        # Initialising energy
        self.e = self.P/(self.gamma - 1)

        # Initialising horizontal and vertical velocities
        self.u = 0*self.R
        self.w = 0*self.R

    def perturbation(self):

        """
        Gaussian temperature perturbation
        """

        A = self.A
        sigma = self.sigma
        x0 = self.x0
        y0 = self.y0

        # Gaussian:
        self.G = A*np.exp(-((self.X-x0)**2+(self.Y-y0)**2)/(2*sigma**2))

    def timestep(self):

        """
        Calculates timestep
        """

        # Function calculating relative change of variable x
        rel = lambda dxdt,x : abs(dxdt/x)

        ## Calculating relative change of primary variables

        rel_rho = rel(self.drho_dt[1:-1],self.rho[1:-1])
        rel_x = rel(self.u[1:-1],self.dx)
        rel_y = rel(self.w[1:-1],self.dy)
        rel_e = rel(self.de_dt[1:-1],self.e[1:-1])

        # Prevents division by zero warning
        np.seterr(divide='ignore', invalid='ignore')

        rel_rhou = np.where(self.rho[1:-1]*self.u[1:-1] < 1e-4, 0,
        rel(self.drhou_dt[1:-1],self.rho[1:-1]*self.u[1:-1]))
        rel_rhow = np.where(self.rho[1:-1]*self.w[1:-1] < 1e-4, 0,
        rel(self.drhow_dt[1:-1],self.rho[1:-1]*self.w[1:-1]))


        delta = np.nanmax([rel_rho,rel_x,rel_y,rel_rhou,rel_rhow,rel_e])

        self.dt = self.p/delta

        # In case the time step is too small or too big:
        if np.isinf(self.dt) or self.dt < 0.001:
            self.dt = 0.001
        elif self.dt > 0.1:
            self.dt = 0.1


    def boundary_conditions(self):

        """
        boundary conditions for energy, density and velocity
        """

        ## Vertical velocity equal to zero
        self.w[0,:] *= 0
        self.w[-1,:] *= 0

        ## Vertical gradient of horizontal velocity equal to zero
        self.u[0,:] = (1/3)*(-self.u[2,:] + 4*self.u[1,:])
        self.u[-1,:] = (1/3)*(-self.u[-3,:] + 4*self.u[-2,:])

        ## Hydrostatic equilibrium at boundary

        # Temperature at boundaries:
        dTdy = -self.mu*self.m_u*self.g*self.nabla/self.kB
        self.T[0,:] = (1/3)*(4*self.T[1,:]-self.T[2,:]-dTdy*2*self.dy)
        self.T[-1,:] = (1/3)*(4*self.T[-2,:]-self.T[-3,:]+dTdy*2*self.dy)

        # Energy at boundaries:
        self.e[0,:] = (4*self.e[1,:] - self.e[2,:])\
        /(3 - self.mu*self.m_u*self.g*2*self.dy/(self.kB*self.T[0,:]))
        self.e[-1,:] = (4*self.e[-2,:] - self.e[-3,:])\
        /(3 + self.mu*self.m_u*self.g*2*self.dy/(self.kB*self.T[-1,:]))

        # Pressure at boundaries:
        self.P[0,:] = (self.gamma - 1)*self.e[0,:]
        self.P[-1,:] = (self.gamma - 1)*self.e[-1,:]

        # Mass density at boundaries:
        self.rho[0,:] = self.mu*self.m_u/(self.kB*self.T[0,:])*self.P[0,:]
        self.rho[-1,:] = self.mu*self.m_u/(self.kB*self.T[-1,:])*self.P[-1,:]


    def central_x(self, var):

        """
        Central difference scheme in x-direction
        """

        var_left = np.roll(var,1,axis=1)
        var_right = np.roll(var,-1,axis=1)

        dvar_dx = (var_right - var_left)/(2*self.dx)

        return dvar_dx

    def central_y(self, var):

        """
        Central difference scheme in y-direction
        """

        var_up = np.roll(var,-1,axis=0)
        var_down = np.roll(var,1,axis=0)

        dvar_dy = (var_up - var_down)/(2*self.dy)

        return dvar_dy

    def upwind_x(self, var, v):

        """
        Upwind difference scheme in x-direction
        """

        var_left = np.roll(var,1,axis=1)
        var_right = np.roll(var,-1,axis=1)

        dvar_dx = np.where(v>=0,
        (var - var_left)/self.dx,
        (var_right - var)/self.dx)

        return dvar_dx


    def upwind_y(self, var, v):

        """
        Upwind difference scheme in y-direction
        """

        var_up = np.roll(var,-1,axis=0)
        var_down = np.roll(var,1,axis=0)

        dvar_dy = np.where(v>=0,
        (var - var_down)/self.dy,
        (var_up - var)/self.dy)

        return dvar_dy

    def hydro_solver(self):

        """
        Hydrodynamic equations solver
        """

        ## Continuity equation:

        self.drho_dt = \
        - self.rho*(self.central_x(self.u) + self.central_y(self.w)) \
        - self.u*self.upwind_x(self.rho,self.u) \
        - self.w*self.upwind_y(self.rho,self.w)

        ## Momentum equation:

        # Horizontal:

        self.drhou_dt = \
        - self.rho*self.u*(self.upwind_x(self.u,self.u) \
        + self.upwind_y(self.w,self.u)) \
        - self.u*self.upwind_x(self.rho*self.u,self.u) \
        - self.w*self.upwind_y(self.rho*self.u,self.w) \
        - self.central_x(self.P)

        # Vertical:

        self.drhow_dt = \
        - self.rho*self.w*(self.upwind_x(self.u,self.w) \
        + self.upwind_y(self.w,self.w)) \
        - self.u*self.upwind_x(self.rho*self.w,self.u) \
        - self.w*self.upwind_y(self.rho*self.w,self.w) \
        - self.central_y(self.P) - self.rho*self.g

        ## Energy equation:

        self.de_dt = \
        - self.u*self.upwind_x(self.e,self.u) \
        - self.w*self.upwind_y(self.e,self.w) \
        - self.gamma*self.e*(self.central_x(self.u) \
        + self.central_y(self.w))

        ## Finding time step size:

        self.timestep()

        dt = self.dt

        ## Calculating new variables

        rho_new = self.rho + self.drho_dt*dt
        u_new = (self.rho*self.u + self.drhou_dt*dt)/rho_new
        w_new = (self.rho*self.w + self.drhow_dt*dt)/rho_new
        e_new = self.e + self.de_dt*dt

        self.rho[1:-1] = rho_new[1:-1]
        self.u[1:-1] = u_new[1:-1]
        self.w[1:-1] = w_new[1:-1]
        self.e[1:-1] = e_new[1:-1]

        self.P[1:-1] = (self.gamma - 1)*self.e[1:-1]
        self.T[1:-1] = self.mu*self.m_u/self.kB*self.P[1:-1]/self.rho[1:-1]

        ## Setting boundary conditions:

        self.boundary_conditions()

        return dt



"""
### SIMULATING:
"""

# Set to True to visualise sanity check (system in equilibrium):
SanityCheck = True

# Set to True to visualise convection (WARNING: takes long to run):
Convection = True

# Set to True to visualise "explosion" (WARNING: takes long to run)
Explosion = True


# Initialising fluid visualiser:
vis = FVis.FluidVisualiser(fontsize = 20)

units = {"Lx":"Mm", "Lz":"Mm", "t":"s"}    # Labeling units
extent = [0, 12, 0, 4]                     # Setting axis extents (in Mm)

if SanityCheck:

    """
    Performs 60 second simulation of gas at hydrostatic equilibrium.
    """

    SC = _2Dconvection(do_pert = False)

    vis.save_data(60, SC.hydro_solver, u = SC.u, w = SC.w, T = SC.T, rho = SC.rho, e = SC.e,
    sim_fps = 1, folder = "SanityCheck")

    vis.animate_2D("T", showQuiver = True, units = units, extent = extent, save = True,
    video_fps = 5, title = "Equilíbrio Hidrostático")

    vis.animate_2D("T", showQuiver = False, units = units, extent = extent, save = True,
    video_fps = 5, title = "Equilíbrio Hidrostático")

if Convection:

    """
    Performs 10 min simulation of gas with temperature perturbation, without
    updating pressure.
    """

    Conv = _2Dconvection(do_pert = True)

    vis.save_data(600, Conv.hydro_solver,rho = Conv.rho, e = Conv.e,
    u = Conv.u, w = Conv.w, T = Conv.T, sim_fps = 1, folder = "Convection")

    # Saving video of temperature and velocity vectors
    vis.animate_2D("T", showQuiver = True, units = units, extent = extent, save = True,
    video_fps = 30, video_name = "Convection_TQuiver", title = "Convecção")

    # Saving video of mass density
    vis.animate_2D("rho", showQuiver = False, units = units,extent = extent, save = True,
    video_fps = 30, video_name = "Convection_rho", title = "Convecção")

    # Saving video of vertical velocity
    vis.animate_2D("w", showQuiver = False, units = units,extent = extent, save = True,
    video_fps = 30, video_name = "Convection_w", title = "Convecção")

if Explosion:

    """
    Performs 10 min simulation of gas with temperature perturbation, updating
    pressure after perturbation.
    """

    Exp = _2Dconvection(do_pert = True)

    # Updating pressure, then mass density and energy density
    Exp.P = Exp.P0*(Exp.T/Exp.T0)**(1/Exp.nabla)
    Exp.rho = Exp.mu*Exp.m_u/(Exp.kB*Exp.T)*Exp.P
    Exp.e = Exp.P/(Exp.gamma - 1)

    vis.save_data(600, Exp.hydro_solver, rho = Exp.rho, e = Exp.e,
    u = Exp.u, w = Exp.w, T = Exp.T, sim_fps = 1, folder = "Explosion")

    # Saving video of temperature and velocity vectors
    vis.animate_2D("T", showQuiver = True, units = units, extent = extent, save = True,
    video_fps = 30, video_name = "Explosion_TQuiver", title = "Explosão")

    # Saving video of mass density
    vis.animate_2D("rho", showQuiver = False, units = units, extent = extent, save = True,
    video_fps = 30, video_name = "Explosion_rho", title = "Explosão")
