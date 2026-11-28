from Beam import Beam
from Matboard import Matboard, u

t = Matboard.t.m_as(u.mm)
length, deck, half_box, half_deck = 1260, 100, 31.5, 50.5
web, inner = half_box - t / 2, half_box - t
tabs = 12, 7

# Underside profile (x, y) in mm.
profile = [(0, 0), (60, 0), (442.5, -17), (817.5, -17), (1200, 0), (1260, 0)]

longitudinal_boards = [
    Matboard([(0, y, z1), (length, y, z1), (length, y, z2), (0, y, z2)], "y")
    for y, z1, z2 in [
        (deck - t / 2, -half_deck, half_deck),
        (deck - 1.5 * t, -half_box - tabs[0], -half_box),
        (deck - 1.5 * t, half_box, half_box + tabs[0]),
    ]
]
longitudinal_boards += [
    Matboard(
        [(0, deck - t, z), (length, deck - t, z)]
        + [(x, y + t, z) for x, y in reversed(profile)],
        "z",
    )
    for z in (-web, web)
]
longitudinal_boards += [
    Matboard(
        [(x, y + offset, z1) for x, y in profile]
        + [(x, y + offset, z2) for x, y in reversed(profile)],
        "y",
    )
    for offset, z1, z2 in [
        (t / 2, -half_box, half_box),
        (1.5 * t, -inner, -inner + tabs[1]),
        (1.5 * t, inner - tabs[1], inner),
    ]
]
longitudinal_boards += [
    Matboard([
        (start, y, 0), (start + 80, y, -inner + 0.25),
        (length - start - 80, y, -inner + 0.25), (length - start, y, 0),
        (length - start - 80, y, inner - 0.25), (start + 80, y, inner - 0.25),
    ], "y")
    for start, y in [(200, deck - 1.5 * t), (320, deck - 2.5 * t)]
]

diaphragm_positions = [x * u.mm for x in range(30, 1231, 120)]
buckling_dimensions = {
    "middle_flange_width": 2 * inner * u.mm,
    "side_flange_width": (half_deck - half_box) * u.mm,
    "web_top": (deck - t) * u.mm,
    "web_height": (deck - min(y for x, y in profile) - 2 * t) * u.mm,
    "diaphragm_positions": diaphragm_positions,
}
glue = [((deck - t) * u.mm, 2 * tabs[0] * u.mm)]
beam_0 = Beam(longitudinal_boards, glue, buckling_dimensions)
