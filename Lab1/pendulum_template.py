#!/bin/python3

# Python simulation of a simple planar pendulum with real time animation
# BH, OF, MP, AJ 2020-10-20, latest version 2021-11-02.
import math

from matplotlib import animation
from pylab import *

"""
    This script defines all the classes needed to simulate (and animate) a single pendulum.
    Hierarchy (somehow in order of encapsulation):
    - Oscillator: a struct that stores the parameters of an oscillator (harmonic or pendulum)
    - Observable: a struct that stores the oscillator's coordinates and energy values over time
    - BaseSystem: harmonic oscillators and pendulums are distinguished only by the expression of
                    the return force. This base class defines a virtual force method, which is
                    specified by its child classes
                    -> Harmonic: specifies the return force as -k*t (i.e. spring)
                    -> Pendulum: specifies the return force as -k*sin(t)
    - BaseIntegrator: parent class for all time-marching schemes; function integrate performs
                    a numerical integration steps and updates the quantity of the system provided
                    as input; function timestep wraps the numerical scheme itself and it's not
                    directly implemented by BaseIntegrator, you need to implement it in his child
                    classes (names are self-explanatory)
                    -> EulerIntegrator: ...
                    -> EulerCromerIntegrator: ...
                    -> VerletIntegrator: ...
                    -> RK4Integrator: ...
    - Simulation: this last class encapsulates the whole simulation procedure; functions are 
                    self-explanatory; you can decide whether to just run the simulation or to
                    run while also producing an animation: the latter option is slower
"""

# Global constants
G = 9.8  # gravitational acceleration

taus=[]

class Oscillator:

    """ Class for a general, simple oscillator """

    def __init__(self, m=1, c=9, t0=0, theta0=0, dtheta0=0):
        self.m = m              # mass of the pendulum bob
        self.c = c              # c = g/L
        self.L = G / c          # string length
        self.t = t0             # the time
        self.theta = theta0     # the position/angle
        self.dtheta = dtheta0   # the velocity


class Observables:

    """ Class for storing observables for an oscillator """

    def __init__(self):
        self.time = []          # list to store time
        self.pos = []           # list to store positions
        self.vel = []           # list to store velocities
        self.energy = []        # list to store energy


class BaseSystem:
    
    def force(self, osc):

        """ Virtual method: implemented by the childclasses  """

        pass


class Harmonic(BaseSystem):
    def force(self, osc):
        return -osc.m*osc.c * osc.theta


class DampedHarmonic(BaseSystem):
    def force(self, osc, omega=3, gamma=0.1):
        return osc.m*(-omega**2*osc.theta-gamma*osc.dtheta)


class Pendulum(BaseSystem):
    def force(self, osc):
        return -osc.m*osc.c * np.sin(osc.theta)


class DampedPendulum(BaseSystem):
    def force(self, osc, gamma=1):
        return -osc.m * gamma * osc.dtheta


class BaseIntegrator:

    def __init__(self, _dt=0.1):
        self.dt = _dt   # time step
        self.should_break = False

    def integrate(self, simsystem, osc, obs):

        """ Perform a single integration step """
        
        self.timestep(simsystem, osc, obs)

        # Append observables to their lists
        obs.time.append(osc.t)
        obs.pos.append(osc.theta)
        obs.vel.append(osc.dtheta)


        # Function 'isinstance' is used to check if the instance of the system object is 'Harmonic' or 'Pendulum'
        if isinstance(simsystem, Harmonic) :
            # Harmonic oscillator energy
            obs.energy.append(0.5 * osc.m * osc.L ** 2 * osc.dtheta ** 2 + 0.5 * osc.m * G * osc.L * osc.theta ** 2)
        else:
            obs.energy.append(0.5 * osc.m * osc.L ** 2 * osc.dtheta ** 2 + 0.5 * osc.m * G * osc.L * (1-math.cos(osc.theta)))

        if len(obs.energy)>1 and (obs.energy[-1] >osc.m*G*osc.L*sin(obs.pos[0])/math.e**2-1 and
                                  obs.energy[-1] <osc.m*G*osc.L*sin(obs.pos[0])/math.e**2+1):
            taus.append(osc.t)
            self.should_break = True

    def timestep(self, simsystem, osc, obs):

        """ Virtual method: implemented by the child classes """

        pass


# HERE YOU ARE ASKED TO IMPLEMENT THE NUMERICAL TIME-MARCHING SCHEMES:

class EulerIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        accel = simsystem.force(osc) / osc.m
        osc.t += self.dt
        try:
            osc.theta += obs.vel[-1]*self.dt
            osc.dtheta -= osc.c*obs.pos[-1]*self.dt
        except IndexError:
            theta0 = osc.theta
            osc.theta += osc.dtheta*self.dt
            osc.dtheta -= osc.c*theta0*self.dt


class EulerCromerIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        accel = simsystem.force(osc) / osc.m
        osc.t += self.dt
        osc.dtheta -= osc.c*osc.theta*self.dt
        osc.theta += osc.dtheta*self.dt


class VerletIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        accel = simsystem.force(osc) / osc.m
        osc.t += self.dt
        osc.theta += osc.dtheta*self.dt + accel*self.dt**2/2
        osc.dtheta += (accel + simsystem.force(osc)/osc.m)*self.dt/2


class RK4Integrator(BaseIntegrator):
    def __init__(self, gamma):
        super().__init__()
        self.gamma=gamma

    def timestep(self, simsystem, osc, obs):
        def accel(theta, dtheta, omega=3, gamma=self.gamma):
            if isinstance(simsystem, Pendulum):
                return -osc.c*np.sin(theta)
            if isinstance(simsystem, DampedHarmonic):
                return -omega**2*theta-gamma*dtheta
            if isinstance(simsystem, Harmonic):
                return -osc.c*theta
            if isinstance(simsystem, DampedPendulum):
                return -osc.c*np.sin(theta)-gamma*dtheta
        osc.t += self.dt
        a1 = accel(osc.theta, osc.dtheta)*self.dt
        b1 = osc.dtheta*self.dt
        a2 = accel(osc.theta+b1/2, osc.dtheta+a1/2)*self.dt
        b2 = (osc.dtheta+a1/2)*self.dt
        a3 = accel(osc.theta+b2/2, osc.dtheta+a2/2)*self.dt
        b3 = (osc.dtheta+a2/2)*self.dt
        a4 = accel(osc.theta+b3, osc.dtheta+a3)*self.dt
        b4 = (osc.dtheta+a3)*self.dt
        osc.dtheta += (a1 + 2*a2 + 2*a3 + a4)/6
        osc.theta += (b1 + 2*b2 + 2*b3 + b4)/6


# Animation function which integrates a few steps and return a line for the pendulum
def animate(framenr, simsystem, oscillator, obs, integrator, pendulum_line, stepsperframe):
    
    for it in range(stepsperframe):
        integrator.integrate(simsystem, oscillator, obs)

    x = np.array([0, np.sin(oscillator.theta)])
    y = np.array([0, -np.cos(oscillator.theta)])
    pendulum_line.set_data(x, y)
    return pendulum_line,


class Simulation:

    def reset(self, osc=Oscillator()):
        self.oscillator = osc
        self.obs = Observables()

    def __init__(self, osc=Oscillator()) :
        self.reset(osc)

    def plot_observables(self, title="simulation"):

        plt.clf()
        plt.title(title)
        plt.plot(self.obs.time, self.obs.pos, label="Position")
        plt.plot(self.obs.time, self.obs.vel, label="Velocity")
        plt.plot(self.obs.time, self.obs.energy, label="Energy")
        plt.xlabel('time')
        plt.ylabel('observables')
        plt.legend()
        #plt.savefig(title + ".pdf")
        plt.show()

    def plot_vel_vs_pos(self, title="Velocity vs position"):
        plt.clf()
        plt.title(title)
        plt.plot(self.obs.pos, self.obs.vel)
        plt.xlabel('Position')
        plt.ylabel('Velocity')
        plt.legend()
        # plt.savefig(title +"theta_vs_dtheta"+ ".pdf")
        plt.show()

    # Run without displaying any animation (fast)
    def run(self,
            simsystem,
            integrator,
            tmax=30.,               # final time
            title="simulation",     # Name of output file and title shown at the top
            show_theta_vs_dtheta=False
            ):

        n = int(tmax / integrator.dt)

        for it in range(n):
            integrator.integrate(simsystem, self.oscillator, self.obs)
            if integrator.should_break:
                break

        self.plot_observables(title)
        if show_theta_vs_dtheta:
            self.plot_vel_vs_pos("Velocity vs position for "+title)

    # Run while displaying the animation of a pendolum swinging back and forth (slow-ish)
    def run_animate(self,
            simsystem,
            integrator,
            tmax=30.,               # final time
            stepsperframe=5,        # how many integration steps between visualising frames
            title="simulation",     # Name of output file and title shown at the top
            ):

        numframes = int(tmax / (stepsperframe * integrator.dt))

        plt.clf()
        # fig = plt.figure()
        ax = plt.subplot(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
        plt.axhline(y=0)  # draw a default hline at y=1 that spans the xrange
        plt.axvline(x=0)  # draw a default vline at x=1 that spans the yrange
        pendulum_line, = ax.plot([], [], lw=5)
        plt.title(title)
        # Call the animator, blit=True means only re-draw parts that have changed
        anim = animation.FuncAnimation(plt.gcf(), animate,  # init_func=init,
                                       fargs=[simsystem,self.oscillator,self.obs,integrator,pendulum_line,stepsperframe],
                                       frames=numframes, interval=25, blit=True, repeat=False)
        plt.show()

        # If you experience problems visualizing the animation and/or
        # the following figures comment out the next line 
        plt.waitforbuttonpress(30)
        
        self.plot_observables(title)

# It's good practice to encapsulate the script execution in 
# a main() function (e.g. for profiling reasons)
def main() :

    # Here you can define one or more instances of oscillators, with possibly different parameters, 
    # and pass them to the simulator 

    # Be sure you are passing the correct initial conditions!
    oscillator = Oscillator(m=1, c=9, theta0=pi*0.5, dtheta0=0)

    # Create the simulation object for your oscillator instance:
    simulation = Simulation(oscillator)

    # Examples of calls to 'simulation.run'
    # simulation.run_animate(simsystem=Harmonic(), integrator=EulerIntegrator(), title="Harmonic-Euler")
    #simulation.run(simsystem=Pendulum(), integrator=EulerCromerIntegrator(), title="Harmonic-EulerCromer")
    #simulation.run(simsystem=Harmonic(), integrator=VerletIntegrator(), title="Harmonic-Verlet")
    #simulation.run(simsystem=Pendulum(), integrator=VerletIntegrator(), title="Pendulum-Verlet")

    #simulation.run(simsystem=Pendulum(), integrator=RK4Integrator(), title="Pendulum-RK4")
    #simulation.run(simsystem=Harmonic(), integrator=RK4Integrator(), title="Harmonic-RK4")
    gammas=[]
    for gamma in [0.5,1,2,3,4,5,6,7]:
        gammas.append(gamma)
        simulation.run(simsystem=DampedHarmonic(), integrator=RK4Integrator(gamma=gamma), title="DampedHarmonic-RK4"+str(gamma))
    #simulation.run(simsystem=DampedPendulum(), integrator=RK4Integrator(), title="DampedPendulum-RK4", show_theta_vs_dtheta=True)
    print(taus)
    plt.clf()
    plt.title('tau vs gamma')
    plt.plot(gammas, taus)
    plt.xlabel('gamma')
    plt.ylabel('tau')
    plt.legend()
    # plt.savefig(title +"theta_vs_dtheta"+ ".pdf")
    plt.show()
# Calling 'main()' if the script is executed.
# If the script is instead just imported, main is not called (this can be useful if you want to
# write another script importing and utilizing the functions and classes defined in this one)
if __name__ == "__main__" :
    main()
