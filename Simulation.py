import compute
import numpy as np
from Beam import Beam
from Matboard import Matboard
from pint import UnitRegistry
from config_0 import beam_0, u
import matplotlib.pyplot as plt

class Simulation:
    axle_positions = np.array([-856, -680,-516, -340, -176, 0]) * u.mm
    
    def __init__(self, loadcase, beam, density = 1250):
        self.loadcase = loadcase
        self.beam = beam
        self.density = density
        min,max = self.beam.get_endpoints()
        self.length = max-min
        self.side = (self.length-1200*u.mm)/2
        self.x = np.linspace(
            min, max, density
        )
        
    def loads(self, position):
        loads = np.zeros(len(self.x)) * u.N
        locations = self.axle_positions + position + self.side
        for x_, force in zip(locations, self.loadcase):
            if self.x[0] + self.side <= x_ <= self.x[-1] - self.side:
                i = np.argmin(np.abs(self.x - x_))
                loads[i] = force
        self.loads = loads
        print("loads are defined.")
        
    def reaction(self,position):
        self.loads(position)
        x_L = self.x[0] + self.side
        self.reac_R = np.sum(
            self.loads * (self.x - x_L)
        ) / (1200 * u.mm)
        self.reac_L = np.sum(self.loads) - self.reac_R
        print("Reaction forces are computed.")
        
    def SFD(self,position):
        self.reaction(position)
        forces = -self.loads.copy()
        i_L = np.argmin(np.abs(self.x - (self.x[0] + self.side)))
        i_R = np.argmin(np.abs(self.x - (self.x[-1] - self.side)))
        forces[i_L] += self.reac_L
        forces[i_R] += self.reac_R
        self.sfd = np.cumsum(forces)
        print("Shear forces are computed.")
        return self.sfd

    def BMD(self,position):
        self.SFD(position)
        dx = np.diff(self.x)
        self.bmd = np.zeros(len(self.x)) * u.N * u.mm
        self.bmd[1:] = np.cumsum(self.sfd[:-1] * dx)
        print("Bending Moments are computed.")
        return self.bmd
    
    def I(self):
        self.I = np.zeros(len(self.x)) * u.mm**4
        for i in range(len(self.x)):
            self.I[i]=self.beam.get_2nd_moment(self.x[i])
        print("Moment of Inertia are computed.")
        return self.I
    
    def y_bounds(self):
        self.y_top = np.zeros(len(self.x))* u.mm
        self.y_bot = np.zeros(len(self.x))* u.mm
        for i in range(len(self.x)):
            bounds = self.beam.get_y_bounds(self.x[i])
            self.y_top[i]=bounds[1]
            self.y_bot[i]=bounds[0]
        print("y boundaries are computed.")
        return self.y_top,self.y_bot
    
    def centroids(self):
        self.centroids = np.zeros(len(self.x)) * u.mm
        for i in range(len(self.x)):
            self.centroids[i]=self.beam.get_centroid(self.x[i])
        print("Centroids are computed.")
        return self.centroids
        
    def calculations(self,position):
        
        self.I()
        self.centroids()
        self.y_bounds()
        
        self.BMD(position)
        
        self.flexural_stress_bot = self.bmd*(self.centroids-self.y_bot)/self.I
        self.flexural_stress_top = self.bmd*(self.y_top-self.centroids)/self.I

    def plot(self):
        x = self.x.to(u.mm).magnitude

        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        # SFD
        axes[0].step(x, self.sfd.to(u.N).magnitude,
                    'r-', where='post')
        axes[0].set_title("Shear Force Diagram")
        axes[0].set_ylabel("Shear Force (N)")

        # BMD
        axes[1].plot(x, self.bmd.to(u.N * u.m).magnitude, 'b-')
        axes[1].set_title("Bending Moment Diagram")
        axes[1].set_ylabel("Bending Moment (N·m)")
        axes[1].set_xlabel("Position (mm)")

        for ax in axes:
            ax.axhline(0, color='black', linewidth=0.8)
            ax.grid(True)

        plt.tight_layout()
        plt.show()
        

        
        
    
    
        
        
        
        
