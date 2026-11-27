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
        
    def cal_loads(self, position):
        loads = np.zeros(len(self.x)) * u.N
        locations = self.axle_positions + position + self.side
        for x_, force in zip(locations, self.loadcase):
            if self.x[0] + self.side <= x_ <= self.x[-1] - self.side:
                i = np.argmin(np.abs(self.x - x_))
                loads[i] = force
        self.loads = loads
        #print("loads are defined.")
        
    def reaction(self,position):
        self.cal_loads(position)
        x_L = self.x[0] + self.side
        self.reac_R = np.sum(
            self.loads * (self.x - x_L)
        ) / (1200 * u.mm)
        self.reac_L = np.sum(self.loads) - self.reac_R
        #print("Reaction forces are computed.")
        
    def SFD(self,position):
        self.reaction(position)
        forces = -self.loads.copy()
        i_L = np.argmin(np.abs(self.x - (self.x[0] + self.side)))
        i_R = np.argmin(np.abs(self.x - (self.x[-1] - self.side)))
        forces[i_L] += self.reac_L
        forces[i_R] += self.reac_R
        self.V = np.cumsum(forces)
        #print("Shear forces are computed.")
        return self.V

    def BMD(self,position):
        self.SFD(position)
        dx = np.diff(self.x)
        self.M = np.zeros(len(self.x)) * u.N * u.mm
        self.M[1:] = np.cumsum(self.V[:-1] * dx)
        #print("Bending Moments are computed.")
        return self.M
    
    def cal_I(self):
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
    
    def cal_centroids(self):
        self.centroids = np.zeros(len(self.x)) * u.mm
        for i in range(len(self.x)):
            self.centroids[i]=self.beam.get_centroid(self.x[i])
        print("Centroids are computed.")
        return self.centroids
    
    def cal_first_moments(self,y):
        first_moments = np.zeros(len(self.x)) * u.mm**3
        for i,y_ in zip(range(len(self.x)),y):
            first_moments[i]=self.beam.get_1st_moment(self.x[i],y_)
        print(f"First moments are computed at level y={y}.")
        return first_moments
        
    def cen_widths(self):
        self.widths = np.zeros(len(self.x)) * u.mm
        for i,c in zip(range(len(self.x)),self.centroids):
            self.widths[i]=self.beam.get_width(self.x[i],c)
        print("Central widths are computed.")
        return self.widths
        
    def calculations_location_independent(self):
        self.cal_I()
        self.cal_centroids()
        self.y_bounds()
        self.cen_widths()
        glue = self.beam.glue[0][0]*np.ones(len(self.x))
        self.glue_width = self.beam.glue[0][1]*np.ones(len(self.x))
        self.first_moments_glue = self.cal_first_moments(glue)
        self.first_moments = self.cal_first_moments(self.centroids)
        
    def calculations(self,position):
        self.BMD(position)
        self.applied_stresses(self.V, self.M)

    def applied_stresses(self, shear, moment):
        self.shear_stress_center = shear*self.first_moments/(self.I*self.widths)
        self.shear_stress_glue = shear*self.first_moments_glue/(self.I*self.glue_width)
        self.flexural_stress_bot = moment*(self.centroids-self.y_bot)/self.I
        self.flexural_stress_top = moment*(self.y_top-self.centroids)/self.I

    def buckling_stresses(self):
        (
            middle_flange_width,
            side_flange_width,
            web_top,
            web_height,
            diaphragm_positions
        ) = self.beam.buckling_dimensions.values()
        t = Matboard.t
        plate_factor = np.pi**2 * Matboard.E_m / (12 * (1 - Matboard.v_m**2))
        middle_flange = 4 * plate_factor * (t / middle_flange_width)**2
        side_flange = 0.425 * plate_factor * (t / side_flange_width)**2
        compressed_web_depth = web_top - self.centroids
        web_compression = 6 * plate_factor * (t / compressed_web_depth)**2

        web_shear = np.ones(len(self.x)) * 5 * plate_factor * (t / web_height)**2
        for left, right in zip(diaphragm_positions, diaphragm_positions[1:]):
            in_panel = (left <= self.x) & (self.x < right)
            web_shear[in_panel] += 5 * plate_factor * (t / (right - left))**2
        return middle_flange, side_flange, web_compression, web_shear
        

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
        

        
        
    
    
        
        
        
        
