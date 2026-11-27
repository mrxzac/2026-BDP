from Beam import Beam
from Matboard import Matboard, u

t = Matboard.t.m_as(u.mm)
length, deck, drop, half_box, half_deck = 1260, 120, 25, 31.0, 50.0
web, inner = half_box - t / 2, half_box - t
tabs = 8, 5

profile = [(0, 0), (60, 0), (420, -drop), (840, -drop), (1200, 0), (1260, 0)]

longitudinal_boards = [
    Matboard([(0, deck - t / 2, -half_deck), (310, deck - t / 2, -half_deck), (310, deck - t / 2, half_deck), (0, deck - t / 2, half_deck)], "y"),
    Matboard([(310, deck - t / 2, -half_deck), (length, deck - t / 2, -half_deck), (length, deck - t / 2, half_deck), (310, deck - t / 2, half_deck)], "y"),
    Matboard([(200, deck - 1.5 * t, -inner), (1060, deck - 1.5 * t, -inner), (1060, deck - 1.5 * t, inner), (200, deck - 1.5 * t, inner)], "y"),
    Matboard([(350, deck - 2.5 * t, -inner), (910, deck - 2.5 * t, -inner), (910, deck - 2.5 * t, inner), (350, deck - 2.5 * t, inner)], "y"),
]
longitudinal_boards += [
    Matboard([(0, deck - 1.5 * t, z1), (length, deck - 1.5 * t, z1), (length, deck - 1.5 * t, z2), (0, deck - 1.5 * t, z2)], "y")
    for z1, z2 in [(-half_box - tabs[0], -half_box), (half_box, half_box + tabs[0])]
]
longitudinal_boards += [
    Matboard([(0, deck - t, z), (length, deck - t, z)] + [(x, y + t, z) for x, y in reversed(profile)], "z")
    for z in (-web, web)
]
longitudinal_boards += [
    Matboard([(x, y + offset, z1) for x, y in profile] + [(x, y + offset, z2) for x, y in reversed(profile)], "y")
    for offset, z1, z2 in [(t / 2, -half_box, half_box), (1.5 * t, -inner, -inner + tabs[1]), (1.5 * t, inner - tabs[1], inner)]
]

diaphragm_x = [30, 130, 230, 330, 430, 530, 630, 730, 830, 930, 1030, 1130, 1230]
diaphragm_y = [0, 0, -11.81, -18.75, -25, -25, -25, -25, -25, -18.75, -11.81, 0, 0]
diaphragms = [
    Matboard([(x, y + 2 * t, -inner), (x, deck - 2 * t, -inner), (x, deck - 2 * t, inner), (x, y + 2 * t, inner)], "x")
    for x, y in zip(diaphragm_x, diaphragm_y)
]

buckling_dimensions = {
    "middle_flange_width": 2 * inner * u.mm,
    "side_flange_width": (half_deck - half_box) * u.mm,
    "web_top": (deck - t) * u.mm,
    "web_height": (deck + drop - 2 * t) * u.mm,
    "diaphragm_positions": [x * u.mm for x in diaphragm_x],
}
glue = [((deck - t) * u.mm, 2 * tabs[0] * u.mm)]

beam_0 = Beam(longitudinal_boards, glue, buckling_dimensions)
beam_1 = Beam(longitudinal_boards + diaphragms, glue, buckling_dimensions)