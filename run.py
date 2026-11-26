from config_0 import beam_0, u
import numpy as np
from Simulation import Simulation

x = 0 * u.mm
print(f"Section x = {x:~P}, centroid y = {beam_0.get_centroid(x):~P}")
print(f"I_y = {beam_0.get_2nd_moment(x):,.2f~P}")

print(beam_0.get_endpoints())

loadcase = np.array([67.5, 67.5, 67.5, 67.5,91, 91]) * u.N
position = (1200 + 856) / 2 * u.mm
sim = Simulation(loadcase,beam_0)
sim.loads(position)
sim.reaction()
sim.SFD()
sim.BMD()
sim.I()
sim.y_bounds()
sim.centroids()
sim.calculations()
print(max(sim.bmd))
print(max(sim.y_top))
print(max(sim.flexural_stress_top))

