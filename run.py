from config_0 import beam_0, glue_shear_strength, u
from Matboard import Matboard
import numpy as np
from Simulation import Simulation

x = 0 * u.mm
print(f"Section x = {x:~P}, centroid y = {beam_0.get_centroid(x):~P}")
print(f"I_y = {beam_0.get_2nd_moment(x):,.2f~P}")

start, end = beam_0.get_endpoints()
print((start, end))

loadcase = np.full(6, 400 / 6) * u.N  # Load Case 1
#loadcase = np.array([67.5, 67.5, 67.5, 67.5, 91, 91]) * u.N  # Load Case 2
sim = Simulation(loadcase, beam_0, density=round((end - start).m_as(u.mm)) + 1)
travel = sim.length - 2 * sim.side - sim.axle_positions.min()
positions = np.linspace(0, travel.m_as(u.mm), round(travel.m_as(u.mm)) + 1) * u.mm
sim.calculations_location_independent()

middle, side, web_compression, web_shear = sim.buckling_stresses()
critical = {}
max_v = np.zeros(len(sim.x)) * u.N
max_m = np.zeros(len(sim.x)) * u.N * u.mm
for position in positions:
    print(f"Train is at x={position:<20,.2f~P}:",end="")
    perc = int(position/positions[-1]*100)
    for i in range(100):
        if i<perc : print("_",end="")
        if i==perc : print("",end="")
        if i>perc : print("_",end="")
    print("")
    sim.calculations(position)
    max_v = np.maximum(max_v, abs(sim.V))
    max_m = np.maximum(max_m, sim.M)

    # Positive bending puts the top in compression and the bottom in tension.
    top, bottom = sim.flexural_stress_top, sim.flexural_stress_bot
    web_shear_demand = abs(sim.shear_stress_center)
    checks = {
        "Matboard tension": (Matboard.s_t, np.maximum(bottom, -top)),
        "Matboard compression": (Matboard.s_c, np.maximum(top, -bottom)),
        "Matboard shear": (Matboard.s_s, web_shear_demand),
        "Glue shear": (glue_shear_strength, abs(sim.shear_stress_glue)),
        "Middle flange buckling": (middle, top),
        "Side flange buckling": (side, top),
        "Web compression buckling": (web_compression, top),
        "Web shear buckling": (web_shear, web_shear_demand),
    }
    for mode, (capacity, demand) in checks.items():
        capacity_mpa = np.broadcast_to(capacity.m_as(u.MPa), sim.x.shape)
        demand_mpa = demand.m_as(u.MPa)
        fos = np.divide(capacity_mpa, demand_mpa, out=np.full(sim.x.shape, np.inf),
                        where=demand_mpa > 0)
        i = np.argmin(fos)
        if mode not in critical or fos[i] < critical[mode][0]:
            critical[mode] = (fos[i], capacity_mpa[i], demand_mpa[i], sim.x[i], position)

print(f"{'Maximum Shear':<20}: {max_v.max():>7,.2f~P}")
print(f"{'Maximum Moment':<20}: {max_m.max():>7,.2f~P}")

print(f"\n{'Failure mode':<25} {'Capacity (MPa)':>15} {'Applied (MPa)':>15} {'FOS':>8}"
      f" {'Beam x (mm)':>12} {'Head at (mm)':>12}")
for mode, (fos, capacity, demand, section, head) in critical.items():
    print(f"{mode:<25} {capacity:15.3f} {demand:15.3f} {fos:8.3f}"
          f" {section.m_as(u.mm):12.1f} {head.m_as(u.mm):27.1f}")

governing = min(critical, key=lambda mode: critical[mode][0])
print(f"\nOverall FOS = {critical[governing][0]:.3f} ({governing})")
