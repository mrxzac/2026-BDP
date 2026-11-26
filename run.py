from config_0 import beam_0, u
import numpy as np
from Simulation import Simulation

x = 0 * u.mm
print(f"Section x = {x:~P}, centroid y = {beam_0.get_centroid('y', x):~P}")
print(f"I_z = {beam_0.get_2nd_moment('z', x):,.2f~P}")

print(beam_0.get_endpoints())

loadcase = np.array([91, 91, 67.5, 67.5, 67.5, 67.5]) * u.N
position = (1200 - 856) / 2 * u.mm
sim = Simulation(loadcase[::-1],beam_0)
sim.loads(position)
sim.reaction()
print(sim.reac_L,sim.reac_R)
sim.SFD()
sim.BMD()
sim.plot()

