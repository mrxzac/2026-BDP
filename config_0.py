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

beam_0 = Beam([
    Matboard([
        (x0,     y1, z1),
        (length, y1, z1),
        (length, y2, z2),
        (x0,     y2, z2),
    ], direction)
    for direction, (y1, z1), (y2, z2) in sections
], [(deck_y, 2 * tab_width)], diaphragm_spacing=400 * u.mm)
