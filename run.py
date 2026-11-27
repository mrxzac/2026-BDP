from config_0 import beam_0, u
import numpy as np
from Simulation import Simulation

x = 0 * u.mm
print(f"Section x = {x:~P}, centroid y = {beam_0.get_centroid(x):~P}")
print(f"I_y = {beam_0.get_2nd_moment(x):,.2f~P}")

print(beam_0.get_endpoints())

loadcase = np.array([67.5, 67.5, 67.5, 67.5, 91, 91]) * u.N
loadcase = np.array([66.7, 66.7, 66.7, 66.7, 66.7, 66.7]) * u.N
position = (1200 + 856) / 2 * u.mm
sim = Simulation(loadcase,beam_0)
max_v = []
max_m = []
max_f_s_t = []
max_f_s_b = []
max_s_s_c = []
max_s_s_g = []
end = beam_0.get_endpoints()[1].magnitude+856
positions = np.linspace(0,end,int(end))*u.mm
sim.calculations_location_independent()

for position in positions:
    print(f"Train head is at x={position}")
    sim.calculations(position)
    max_m.append(max(sim.M))
    max_v.append(max(sim.V))
    max_f_s_t.append(max(sim.flexural_stress_top))
    max_f_s_b.append(max(sim.flexural_stress_bot))
    max_s_s_c.append(max(sim.shear_stress_center))
    max_s_s_g.append(max(sim.shear_stress_glue))
    
print(f"Maximum Shear: {max(max_v):,.2f~P}")
print(f"Maximum Moment: {max(max_m):,.2f~P}")
print(f"Maximum Flexural Stress at top: {max(max_f_s_t).to("MPa"):,.2f~P}")
print(f"Maximum Flexural Stress at bottom: {max(max_f_s_b).to("MPa"):,.2f~P}")
print(f"Maximum Shear Stress at center: {max(max_s_s_c).to("MPa"):,.2f~P}")
print(f"Maximum Shear Stress at glue: {max(max_s_s_g).to("MPa"):,.2f~P}")

