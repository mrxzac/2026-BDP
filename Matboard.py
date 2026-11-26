from pint import UnitRegistry
import matplotlib.pyplot as plt
import math
from compute import geometric_centroid, polygon_area

u = UnitRegistry()
d = {'x':0,'y':1,'z':2}

class Matboard:
    s_t = 30 * u.MPa 
    """ tensile stress """
    s_c = 6 * u.MPa     
    """ compressive stress """
    s_s = 4 * u.MPa     
    """ shear stress """
    E_m = 4000 * u.MPa  
    """ Young's Modulus """
    v_m = 0.2           
    """ Poisson's Ratio """
    t = 1.27 * u.mm   
    """ thickness """

    def __init__(self, vertices, dirc):
        self.dirc = dirc
        self.vertices = vertices
            
    def get_sec(self, dirc, val):
        
        dim = d[dirc]
        dim_ = d[self.dirc]

        if dirc == self.dirc:
            board_pos = self.vertices[0][dim]
            if board_pos - self.t / 2 <= val <= board_pos + self.t / 2:
                return [tuple(val if j == dim else coord for j, coord in enumerate(p))
                        for p in self.vertices]
            return []

        c1, c2 = None, None
        n = len(self.vertices)
        for i in range(n):
            a = self.vertices[i]
            b = self.vertices[(i + 1) % n]

            if a[dim] == b[dim] or not min(a[dim], b[dim]) <= val <= max(a[dim], b[dim]):
                continue
            t = (val - a[dim]) / (b[dim] - a[dim])
            point = [a[j] + t * (b[j] - a[j]) for j in range(3)]
            if a[dim] > b[dim]:
                c1 = point
            else:
                c2 = point

        if c1 is None or c2 is None or c1 == c2:
            return []

        c1_l = [c - self.t/2 if i == dim_ else c for i, c in enumerate(c1)]
        c2_l = [c - self.t/2 if i == dim_ else c for i, c in enumerate(c2)]
        c1_r = [c + self.t/2 if i == dim_ else c for i, c in enumerate(c1)]
        c2_r = [c + self.t/2 if i == dim_ else c for i, c in enumerate(c2)]

        return [c1_l, c1_r, c2_r, c2_l]
    
    def get_sec_area(self, dirc, val):
        return polygon_area(self.get_sec(dirc,val))
    
    def get_sec_center(self, dirc, val):
        return geometric_centroid(self.get_sec(dirc,val))
    
    
    
