from config_0 import beam_0, u

x = 0 * u.mm
print(f"Section x = {beam_0.get_centroid('y', x):~P}")
print(f"I_z = {beam_0.get_2nd_moment('z', x):,.2f~P}")

