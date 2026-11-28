from Matboard import Matboard, u
import numpy as np
from Simulation import Simulation
import matplotlib.pyplot as plt

from config_trapezoid_3 import beam_0

x = 0 * u.mm
print(f"Section x = {x:~P}, centroid y = {beam_0.get_centroid(x):~P}")
print(f"I_y = {beam_0.get_2nd_moment(x):,.2f~P}")

start, end = beam_0.get_endpoints()
print((start, end))

loadcase1 = np.full(6, 400 / 6) * u.N  # Load Case 1
loadcase2 = np.array([67.5, 67.5, 67.5, 67.5, 91, 91]) * u.N  # Load Case 2
loadcase3 = np.array([77,77,70,70,107.03,107.03]) *u.N
sim = Simulation(loadcase3, beam_0, density=round((end - start).m_as(u.mm)) + 1)

travel = sim.length - 2 * sim.side - sim.axle_positions.min()
positions = np.linspace(0, travel.m_as(u.mm), round(travel.m_as(u.mm)) + 1) * u.mm
max_v = np.zeros(len(sim.x)) * u.N
max_m = np.zeros(len(sim.x)) * u.N * u.mm

fig,  (ax1, ax2) = plt.subplots(1,2)

sim.calculations_location_independent()

moment_envelope = np.zeros(len(sim.x)) * u.mm*u.N
shear_envelope = np.zeros(len(sim.x)) * u.N
count = 0

for position in positions:
    count+=1
    print(f"Train is at x={position:<20,.2f~P}:",end="")
    perc = int(position/positions[-1]*50)
    for i in range(50):
        if i<perc : print("_",end="")
        if i==perc : print("",end="")
        if i>perc : print("_",end="")
    print("")
    sim.BMD(position)
    shear = abs(sim.V)
    max_v = np.maximum(max_v, shear)
    max_m = np.maximum(max_m, sim.M)
    moment_envelope = np.maximum(moment_envelope,sim.M)
    shear_envelope = np.maximum(shear_envelope,np.abs(sim.V))
    
    if count%100 == 0: ax1.plot(sim.x.m_as("m"),sim.M.m_as("N*m"),color='gray',label='_nolegend_',alpha=0.3)
    if count%100 == 0: ax2.plot(sim.x.m_as("m"),sim.V.m_as("N"),color='gray',label='_nolegend_',alpha=0.3)

ax1.plot(sim.x.m_as("m"),moment_envelope.m_as("N*m"),'b--',label='Envelope')
ax2.plot(sim.x.m_as("m"),shear_envelope.m_as("N"),'b--',label='Envelope')
sim.applied_stresses(max_v, max_m)
middle, side, web_compression, web_shear = sim.buckling_stresses()

top, bottom = sim.flexural_stress_top, sim.flexural_stress_bot
web_shear_demand = abs(sim.shear_stress_center)

checks_M = {
    "Matboard tension": (Matboard.s_t, np.maximum(bottom, -top), ),
    "Matboard compression": (Matboard.s_c, np.maximum(top, -bottom)),
    "Middle flange buckling": (middle, top),
    "Side flange buckling": (side, top),
    "Web compression buckling": (web_compression, top),
}
checks_V = {
    "Matboard shear": (Matboard.s_s, web_shear_demand),
    "Glue shear": (Matboard.s_g, abs(sim.shear_stress_glue)),
    "Web shear buckling": (web_shear, web_shear_demand),
}

for mode, (capacity, demand) in checks_M.items():
    demand_mpa = np.abs(demand.m_as(u.MPa))
    capacity_mpa = np.abs(np.broadcast_to(capacity.m_as(u.MPa), demand_mpa.shape))
    fos = np.divide(
        capacity_mpa,
        demand_mpa,
        out=np.full(demand_mpa.shape, np.inf),
        where=demand_mpa > 0,
    )
    M_fail = moment_envelope * fos
    ax1.plot(sim.x.m_as("m"),M_fail.m_as("N*m"),label=mode)

for mode, (capacity, demand) in checks_V.items():
    demand_mpa = np.abs(demand.m_as(u.MPa))
    capacity_mpa = np.abs(np.broadcast_to(capacity.m_as(u.MPa), demand_mpa.shape))
    fos = np.divide(
        capacity_mpa,
        demand_mpa,
        out=np.full(demand_mpa.shape, np.inf),
        where=demand_mpa > 0,
    )
    V_fail = shear_envelope * fos
    ax2.plot(sim.x.m_as("m"),V_fail.m_as("N"),label=mode)

checks = checks_V | checks_M
critical = {}
for mode, (capacity, demand) in checks.items():
    demand_mpa = demand.m_as(u.MPa)
    capacity_mpa = np.broadcast_to(capacity.m_as(u.MPa), demand_mpa.shape)
    fos = np.divide(
        capacity_mpa,
        demand_mpa,
        out=np.full(demand_mpa.shape, np.inf),
        where=demand_mpa > 0,
    )
    at_x = np.argmin(fos)
    critical[mode] = (fos[at_x],capacity_mpa[at_x],demand_mpa[at_x],sim.x[at_x],)

beam_dimensions = {
    "Centriodal Axis": (np.max(sim.centroids),np.min(sim.centroids),np.mean(sim.centroids)),
    "I Moment of Inertia": (np.max(sim.I),np.min(sim.I),np.mean(sim.I)),
    "First Moment at center": (np.max(sim.first_moments),np.min(sim.first_moments),np.mean(sim.first_moments)),
    "First Moment at glue": (np.max(sim.first_moments_glue),np.min(sim.first_moments_glue),np.mean(sim.first_moments_glue)),
}

print(f"\n{'Beam Dimensions':<25} {'Max':>15} {'Min':>15} {'Average':>15}")
for name, (max_,min_,avg_) in beam_dimensions.items():
    print(f"{name:<25} {max_:15.3f~P} {min_:15.3f~P} {avg_:15.3f~P}")
print("")

print(f"{'Maximum Shear':<20}: {max_v.max():>7,.2f~P}")
print(f"{'Maximum Moment':<20}: {max_m.max():>7,.2f~P}")
print(f"\n{'Failure mode':<25} {'Capacity (MPa)':>15} {'Applied (MPa)':>15} {'FOS':>8} {'Beam x (mm)':>12}")
for mode, (fos, capacity, demand, section) in critical.items():
    print(f"{mode:<25} {capacity:15.3f} {demand:15.3f} {fos:8.3f} {section.m_as(u.mm):12.1f}")


governing = min(critical, key=lambda mode: critical[mode][0])
print(f"\nOverall FOS = {critical[governing][0]:.3f} ({governing})")
print(f"Predicted failure load = {sim.loadcase.sum() * critical[governing][0]:.1f~P}")

ax1.set_xlabel("Length m")
ax1.set_ylabel(f"Torque Nm")
ax1.legend(loc="upper right")

ax2.set_xlabel("Length m")
ax2.set_ylabel(f"Shear N")
ax2.legend(loc="upper right")
plt.show()