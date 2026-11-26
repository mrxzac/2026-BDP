"""Tapered box girder compatible with the supplied Beam/Matboard/strength API.
Coordinates in mm. Nominal strength model: exactly two continuous webs.
Splice plates are recorded separately; their contribution and strength are
not credited by the current two-web strength solver. Not a fabrication approval.
"""
import numpy as np
from Beam import Beam
from Matboard import Matboard, u

x0 = 0 * u.mm
length = 1250 * u.mm
span = 1200 * u.mm
t = Matboard.t
T = t.to(u.mm).magnitude
L = 1250.0
H = 80.0                      # exposed deck elevation above supports
B = 60.0                      # outside box width
W = 100.0                     # deck width
TAB = 6.0
web_top = H - 2*T
web_z = B/2 - T/2
inner_z = B/2 - T
notch_z = inner_z - TAB

# Horizontal bearing regions, tapered haunches, deep central section.
profile = [(0, 0), (60, 0), (525, -50), (725, -50),
           (1190, 0), (1250, 0)]


def bottom(x):
    for (a, ya), (b, yb) in zip(profile, profile[1:]):
        if a <= x <= b:
            return ya + (yb-ya)*(x-a)/(b-a)
    raise ValueError('x outside bridge')


def knots(a, b):
    return [a] + [x for x, _ in profile if a < x < b] + [b]


class _LongitudinalPanel(Matboard):
    """A piece of a continuous plate with half-open x ownership.

    At an internal profile break, the segment on the right owns the section.
    This prevents duplicate webs/flanges at x=60, 525, 725 and 1190 mm.
    The final segment includes the bridge's rightmost endpoint.
    """
    def get_sec_mm(self, dirc, val):
        if dirc == 'x':
            x = u.Quantity(val, u.mm).magnitude
            a = self.vertices[:, 0].min()
            b = self.vertices[:, 0].max()
            if x < a or x > b or (x == b and b != L):
                return np.empty((0, 3))
        return super().get_sec_mm(dirc, val)


def board(vertices, direction):
    # Numeric vertices are already millimetres, as accepted by compute.py.
    cls = Matboard if direction == 'x' else _LongitudinalPanel
    return cls(np.asarray(vertices, dtype=float), direction)


def horizontal(a, b, y, z1, z2):
    return board([(a,y,z1), (b,y,z1), (b,y,z2), (a,y,z2)], 'y')


def web(a, b, z):
    # Split at each slope break: every individual web is a convex
    # rectangle or trapezoid with exactly four vertices.
    p = knots(a, b)
    return [board([(x, bottom(x)+T, z), (y, bottom(y)+T, z),
                   (y, web_top, z), (x, web_top, z)], 'z')
            for x, y in zip(p, p[1:])]


def bottom_strips(a, b, z1, z2, offset):
    # Each face is planar but some are sloped. The y-thickness convention
    # approximates normal thickness on the ramps (slope ~0.108).
    p = knots(a,b)
    return [board([(x,bottom(x)+offset,z1), (y,bottom(y)+offset,z1),
                   (y,bottom(y)+offset,z2), (x,bottom(x)+offset,z2)], 'y')
            for x,y in zip(p,p[1:])]


def web_with_tabs(a,b,z):
    z1,z2 = ((-inner_z,-notch_z) if z < 0 else (notch_z,inner_z))
    return web(a,b,z) + [horizontal(a,b,web_top-T/2,z1,z2)] + \
           bottom_strips(a,b,z1,z2,1.5*T)


# Analysis geometry split at profile breaks. Adjacent pieces of a physical
# web are subdivisions, not extra cuts. Physical blanks have folded tabs.
longitudinal_boards = [horizontal(0,L,H-T/2,-W/2,W/2),
                       horizontal(0,L,H-1.5*T,-W/2,W/2)]
longitudinal_boards += bottom_strips(0,L,-B/2,B/2,T/2)
for z in [-web_z,web_z]:
    longitudinal_boards += web_with_tabs(0,L,z)

# Staggered cut seams. Bottom strips at profile bends are faces of the
# same folded physical blank, not extra independent cut pieces.
cut_boards = []
for y,seam in [(H-T/2,725),(H-1.5*T,525)]:
    for a,b in [(0,seam),(seam,L)]:
        cut_boards.append(horizontal(a,b,y,-W/2,W/2))
for z,seam in [(-web_z,525),(web_z,725)]:
    for a,b in [(0,seam),(seam,L)]:
        cut_boards += web_with_tabs(a,b,z)
for a,b in [(0,400),(400,850),(850,L)]:
    cut_boards += bottom_strips(a,b,-B/2,B/2,T/2)

left_support = 25*u.mm
right_support = 1225*u.mm
diaphragm_stations = [10,40,175,310,445,575,675,805,940,1075,1210,1240]
diaphragm_centres = [x*u.mm for x in diaphragm_stations]


def diaphragm(x):
    # The physical notched diaphragm is one piece. Represent it for the
    # section solver as three non-overlapping rectangular regions.
    lo = max(bottom(x-T/2), bottom(x+T/2))+T
    hi = web_top
    def rect(y1, y2, z1, z2):
        return board([(x,y1,z1), (x,y2,z1), (x,y2,z2), (x,y1,z2)], 'x')
    return [rect(lo, hi, -notch_z, notch_z),
            rect(lo+T, hi-T, -inner_z, -notch_z),
            rect(lo+T, hi-T, notch_z, inner_z)]


diaphragms = [p for x in diaphragm_stations for p in diaphragm(x)]
web_splices = [p for x,z in [(525,-B/2-T/2),(725,B/2+T/2)]
               for p in web(x-50,x+50,z)]
deck_splices = [horizontal(x-50,x+50,web_top-T/2,-notch_z,notch_z)
                for x in [525,725]]
bottom_splices = [p for x in [400,850] for p in
                  bottom_strips(x-50,x+50,-notch_z,notch_z,1.5*T)]
reinforcements = web_splices+deck_splices+bottom_splices
physical_boards = cut_boards+diaphragms+reinforcements

# strength._length_at already accepts callables receiving a Pint length.
def lower_glue_y(x):
    return (bottom(u.Quantity(x, u.mm).magnitude) + T) * u.mm


glue_joints = [
    (web_top*u.mm, 2*TAB*u.mm),
    ((H-T)*u.mm, W*u.mm),
    (lower_glue_y, 2*TAB*u.mm),
]


def glue_joints_at(x):
    return [(h(u.Quantity(x, u.mm)) if callable(h) else h, b) for h,b in glue_joints]


# Keep splice doublers in physical_boards, not the nominal strength model:
# strength._section explicitly requires exactly two vertical web strips.
# The calculation assumes fully effective splices, without credit for doublers.
beam_0 = Beam(longitudinal_boards + diaphragms, glue_joints)
beam_1 = beam_0

# Optional visualization only; do NOT pass this object to strength.analyze.
physical_beam = Beam(physical_boards, glue_joints)
