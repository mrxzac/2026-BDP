from Beam import Beam
from Matboard import Matboard, u

t = Matboard.t.m_as(u.mm)
length, deck, drop, half_box, half_deck = 1260, 120, 25, 31.0, 50.0
web, inner = half_box - t / 2, half_box - t
tabs = 8, 5

y_150 = -drop * (150 - 60) / 360
y_270 = -drop * (270 - 60) / 360
y_310 = -drop * (310.5 - 60) / 360
y_390 = -drop * (390 - 60) / 360
y_870 = -drop * (1200 - 870) / 360
y_950 = -drop * (1200 - 950.5) / 360
y_990 = -drop * (1200 - 990) / 360
y_1030 = -drop * (1200 - 1030) / 360
y_1110 = -drop * (1200 - 1110) / 360

web_L = [(0, 0), (60, 0), (420, -drop), (630.5, -drop)]
web_R = [(630.5, -drop), (840, -drop), (1200, 0), (1260, 0)]
bot_L = [(0, 0), (60, 0), (420, -drop), (840, -drop), (950.5, y_950)]
bot_R = [(950.5, y_950), (1200, 0), (1260, 0)]
tab_L = [(0, 0), (60, 0), (310.5, y_310)]
tab_R = [(310.5, y_310), (420, -drop), (840, -drop), (1200, 0), (1260, 0)]
bot_splice = [(870, y_870), (950.5, y_950), (1030, y_1030)]

longitudinal_boards = [
    Matboard([(0, deck - t / 2, -half_deck), (310.5, deck - t / 2, -half_deck), (310.5, deck - t / 2, half_deck), (0, deck - t / 2, half_deck)], "y"),
    Matboard([(310.5, deck - t / 2, -half_deck), (length, deck - t / 2, -half_deck), (length, deck - t / 2, half_deck), (310.5, deck - t / 2, half_deck)], "y"),
    Matboard([(200, deck - 1.5 * t, -inner), (1060, deck - 1.5 * t, -inner), (1060, deck - 1.5 * t, inner), (200, deck - 1.5 * t, inner)], "y"),
    Matboard([(410, deck - 2.5 * t, -inner), (850, deck - 2.5 * t, -inner), (850, deck - 2.5 * t, inner), (410, deck - 2.5 * t, inner)], "y"),
]
longitudinal_boards += [
    Matboard([(x1, deck - 1.5 * t, z1), (x2, deck - 1.5 * t, z1), (x2, deck - 1.5 * t, z2), (x1, deck - 1.5 * t, z2)], "y")
    for x1, x2 in [(0, 950.5), (950.5, length)]
    for z1, z2 in [(-half_box - tabs[0], -half_box), (half_box, half_box + tabs[0])]
]
longitudinal_boards += [
    Matboard([(x0, deck - t, z), (x1, deck - t, z)] + [(x, y + t, z) for x, y in reversed(p)], "z")
    for z in (-web, web)
    for x0, x1, p in [(0, 630.5, web_L), (630.5, length, web_R)]
]
longitudinal_boards += [
    Matboard([(570, deck - 3 * t, z), (690, deck - 3 * t, z), (690, -drop + 2 * t, z), (570, -drop + 2 * t, z)], "z")
    for z in (-inner + t / 2, inner - t / 2)
]
longitudinal_boards += [
    Matboard([(x, y + t / 2, -half_box) for x, y in p] + [(x, y + t / 2, half_box) for x, y in reversed(p)], "y")
    for p in (bot_L, bot_R)
]
longitudinal_boards += [
    Matboard([(x, y + 1.5 * t, -inner + tabs[1]) for x, y in bot_splice] + [(x, y + 1.5 * t, inner - tabs[1]) for x, y in reversed(bot_splice)], "y")
]
longitudinal_boards += [
    Matboard([(x, y + 1.5 * t, z1) for x, y in p] + [(x, y + 1.5 * t, z2) for x, y in reversed(p)], "y")
    for z1, z2 in [(-inner, -inner + tabs[1]), (inner - tabs[1], inner)]
    for p in (tab_L, tab_R)
]

diaphragm_specs = [
    (30, 0, deck - 2 * t, inner),
    (150, y_150, deck - 2 * t, inner),
    (270, y_270, deck - 2 * t, inner),
    (390, y_390, deck - 2 * t, inner),
    (510, -drop, deck - 3 * t, inner),
    (630, -drop, deck - 3 * t, inner - t),
    (750, -drop, deck - 3 * t, inner),
    (870, y_870, deck - 2 * t, inner),
    (990, y_990, deck - 2 * t, inner),
    (1110, y_1110, deck - 2 * t, inner),
    (1230, 0, deck - 2 * t, inner),
]
diaphragms = [
    Matboard([(x, y + 2 * t, -z), (x, yt, -z), (x, yt, z), (x, y + 2 * t, z)], "x")
    for x, y, yt, z in diaphragm_specs
]

buckling_dimensions = {
    "middle_flange_width": 2 * inner * u.mm,
    "side_flange_width": (half_deck - half_box) * u.mm,
    "web_top": (deck - t) * u.mm,
    "web_height": (deck + drop - 2 * t) * u.mm,
    "diaphragm_positions": [x * u.mm for x, _, _, _ in diaphragm_specs],
}
glue = [((deck - t) * u.mm, 2 * tabs[0] * u.mm)]

beam_0 = Beam(longitudinal_boards, glue, buckling_dimensions)
beam_1 = Beam(longitudinal_boards + diaphragms, glue, buckling_dimensions)