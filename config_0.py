
from Beam import Beam
from Matboard import Matboard, u

x0 = 0 * u.mm
length = 1250 * u.mm
t = Matboard.t

deck_y = 75 * u.mm
half_top = 100 * u.mm / 2
half_box = 80 * u.mm / 2
web_z = half_box - t / 2
inner_z = half_box - t
tab_y = deck_y - t / 2
tab_width = 5 * u.mm

sections = [
    ("y", (deck_y + t / 2, -half_top), (deck_y + t / 2, half_top)),
    ("y", (t / 2, -half_box), (t / 2, half_box)),
    ("z", (t, -web_z), (deck_y, -web_z)),
    ("z", (t, web_z), (deck_y, web_z)),
    ("y", (tab_y, -inner_z), (tab_y, -inner_z + tab_width)),
    ("y", (tab_y, inner_z - tab_width), (tab_y, inner_z)),
]

#boards placed
longitudinal_boards = [
    Matboard([
        (x0,     y1, z1),
        (length, y1, z1),
        (length, y2, z2),
        (x0,     y2, z2),
    ], direction)
    for direction, (y1, z1), (y2, z2) in sections
]

span = 1200 * u.mm
support_L = (length - span) / 2
support_R = support_L + span

offset = 15 * u.mm

#diaphragm positions along the beam
diaphragm_positions = [
    support_L - offset,
    support_L + offset,
    support_L + span / 3,
    support_L + span * 2 / 3,
    support_R - offset,
    support_R + offset,
]

#choose the critical buckling dimension for the beam
#i.e. the longest part or other critical axle_positions
#only these types of buckling is considered
buckling_dimensions = {
    "middle_flange_width": 2 * inner_z,
    "side_flange_width": half_top - half_box,
    "web_top": deck_y,
    "web_height": deck_y,
    "diaphragm_positions": diaphragm_positions,
}
glue_shear_strength = 2 * u.MPa 

#diaphragm_section
#will be considered in the moment of inertia calculation 
#if combined into the Beam (optional)
# diaphragm_section = [
#     (t,          -inner_z),
#     (t,           inner_z),
#     (deck_y - t,  inner_z),
#     (deck_y - t,  inner_z - tab_width),
#     (deck_y,      inner_z - tab_width),
#     (deck_y,     -inner_z + tab_width),
#     (deck_y - t, -inner_z + tab_width),
#     (deck_y - t, -inner_z),
# ]

# diaphragms = [
#     Matboard([
#         (x, y, z)
#         for y, z in diaphragm_section
#     ], "x")
#     for x in diaphragm_positions
# ]

glue = [
    (deck_y, 2 * (tab_width+t)),
]

beam_0 = Beam(longitudinal_boards, glue, buckling_dimensions )
