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
    ("y", (deck_y + t / 2, -half_top), (deck_y + t / 2, half_top)),  # top
    ("y", (t / 2, -half_box), (t / 2, half_box)),                    # bottom
    ("z", (t, -web_z), (deck_y, -web_z)),                            # left web
    ("z", (t, web_z), (deck_y, web_z)),                              # right web
    ("y", (tab_y, -inner_z), (tab_y, -inner_z + tab_width)),         # left tab
    ("y", (tab_y, inner_z - tab_width), (tab_y, inner_z)),           # right tab
]

longitudinal_boards = [
    Matboard([
        (x0,     y1, z1),
        (x0 + length, y1, z1),
        (x0 + length, y2, z2),
        (x0,     y2, z2),
    ], direction)
    for direction, (y1, z1), (y2, z2) in sections
]

# Appendix II dimensions the inside stations from the support centrelines.
# The pictured pair over each support is taken as +/-15 mm (undimensioned).
span = 1200 * u.mm
left_support = x0 + (length - span) / 2
right_support = left_support + span
pair_offset = 15 * u.mm
diaphragm_centres = [
    left_support - pair_offset, left_support + pair_offset,
    left_support + span / 3, right_support - span / 3,
    right_support - pair_offset, right_support + pair_offset,
]

# The plate fills the box void. Two 5 mm shoulders fit under the glue tabs.
shoulder_y = deck_y - t
notch_z = inner_z - tab_width
diaphragm_outline = [
    (t, -inner_z), (t, inner_z),
    (shoulder_y, inner_z), (shoulder_y, notch_z),
    (deck_y, notch_z), (deck_y, -notch_z),
    (shoulder_y, -notch_z), (shoulder_y, -inner_z),
]
diaphragms = [
    Matboard([(x, y, z) for y, z in diaphragm_outline], "x")
    for x in diaphragm_centres
]

beam_0 = Beam(longitudinal_boards + diaphragms, [(deck_y, 2 * tab_width)])
