"""Elastic beam and CIV102 thin-plate checks for a two-web matboard box girder.

Internal units are mm, N and MPa. Positive moment compresses the upper flange.
The box has two vertical webs and continuous horizontal top and bottom flanges.
Longitudinal plates carry bending and shear. Transverse diaphragms set the web
panel length through beam.diaphragm_spacing (a length or a function of x).
"""

from dataclasses import dataclass
from math import inf, pi

import numpy as np

from compute import u
from Matboard import Matboard


MODE_LABELS = {
    'tension': 'Matboard tension',
    'compression': 'Matboard compression',
    'shear': 'Matboard shear',
    'glue': 'Glue shear',
    'middle_flange': 'Middle flange buckling',
    'outstand': 'Side flange buckling',
    'web_compression': 'Web compression buckling',
    'web_shear': 'Web shear buckling',
}
SHEAR_MODES = {'shear', 'glue', 'web_shear'}
GLUE_STRENGTH = 2.0  # MPa for fully cured contact cement
PLATE_FACTOR = pi**2 * Matboard.E_m.m_as(u.MPa) / (12 * (1 - Matboard.v_m**2))


@dataclass(frozen=True)
class Strip:
    normal: str
    y0: float
    y1: float
    z0: float
    z1: float

    @property
    def width(self):
        return self.z1 - self.z0

    @property
    def height(self):
        return self.y1 - self.y0

    @property
    def y(self):
        return (self.y0 + self.y1) / 2

    @property
    def area(self):
        return self.width * self.height

    def first_moment_above(self, elevation, neutral_axis):
        cut = min(max(elevation, self.y0), self.y1)
        return self.width * (self.y1 - cut) * ((self.y1 + cut) / 2 - neutral_axis)


@dataclass(frozen=True)
class Limit:
    strength: float          # MPa
    stress_per_unit: float   # MPa per N or per N mm

    @property
    def capacity(self):
        return self.strength / self.stress_per_unit if self.stress_per_unit else inf


ABSENT = Limit(inf, 0.)


def _weakest(limits):
    return min(limits, key=lambda check: check.capacity, default=ABSENT)


def _bending_checks(strips, webs, ybar, inertia, sign):
    """Plate widths run to web centrelines. Only full-panel plies add buckling thickness."""
    top, bottom = max(p.y1 for p in strips), min(p.y0 for p in strips)
    compression_arm, tension_arm = (top - ybar, ybar - bottom)
    if sign < 0:
        compression_arm, tension_arm = tension_arm, compression_arm
    supports = [(p.z0 + p.z1) / 2 for p in webs]
    left, right = supports
    flange = [p for p in strips if p.normal == 'y' and sign * (p.y - ybar) > 0]
    continuous = [p for p in flange if p.z0 <= left and right <= p.z1]
    if not continuous:
        raise ValueError("Each flange must span both web centrelines.")
    mid_t = sum(p.height for p in continuous)
    middle = Limit(4 * PLATE_FACTOR * (mid_t / (right - left))**2, compression_arm / inertia)

    outstands = []
    for side, support, face in [(-1, left, webs[0].z0), (1, right, webs[1].z1)]:
        edge = min(p.z0 for p in flange) if side < 0 else max(p.z1 for p in flange)
        if side * (edge - face) > 1e-9:
            t = sum(p.height for p in flange if p.z0 <= min(edge, support)
                    and p.z1 >= max(edge, support))
            if not t:
                raise ValueError("An outstand needs a plate continuous to its supporting web.")
            outstands.append(Limit(0.425 * PLATE_FACTOR * (t / (edge - support))**2,
                                   compression_arm / inertia))
    web_checks = []
    for web in webs:
        arm = web.y1 - ybar if sign > 0 else ybar - web.y0
        depth = min(arm, web.height)
        if depth > 0:
            web_checks.append(Limit(6 * PLATE_FACTOR * (web.width / depth)**2, arm / inertia))
    return {
        'tension': Limit(Matboard.s_t.m_as(u.MPa), tension_arm / inertia),
        'compression': Limit(Matboard.s_c.m_as(u.MPa), compression_arm / inertia),
        'middle_flange': middle,
        'outstand': _weakest(outstands),
        'web_compression': _weakest(web_checks),
    }


def _section(strips, joints, panel_length):
    webs = sorted((p for p in strips if p.normal == 'z'), key=lambda p: p.z0)
    if len(webs) != 2 or webs[0].z1 >= webs[1].z0:
        raise ValueError("Box-strength checks require two distinct vertical webs.")
    area = sum(p.area for p in strips)
    ybar = sum(p.area * p.y for p in strips) / area
    inertia = sum(p.area * (p.height**2 / 12 + (p.y - ybar)**2) for p in strips)
    if not all(web.y0 < ybar < web.y1 for web in webs):
        raise ValueError("The neutral axis must cross both webs for the web shear model.")
    width = sum(p.width for p in strips if p.y0 < ybar < p.y1)
    q = sum(p.first_moment_above(ybar, ybar) for p in strips)
    shear_per_n = q / (inertia * width)
    properties = {'area_mm2': area, 'ybar_mm': ybar, 'Iz_mm4': inertia,
                  'y_top_mm': max(p.y1 for p in strips) - ybar,
                  'y_bottom_mm': ybar - min(p.y0 for p in strips),
                  'Q_NA_mm3': q, 'b_NA_mm': width, 'panel_length_mm': panel_length}
    glue = []
    for i, (height, bonded_width) in enumerate(joints, 1):
        if not min(p.y0 for p in strips) < height < max(p.y1 for p in strips) or bonded_width <= 0:
            raise ValueError("Glue joints need an interior elevation and a positive total bonded width.")
        q_joint = sum(p.first_moment_above(height, ybar) for p in strips)
        properties.update({f'Q_glue_{i}_mm3': q_joint, f'b_glue_{i}_mm': bonded_width})
        glue.append(Limit(GLUE_STRENGTH, q_joint / (inertia * bonded_width)))
    if not np.isfinite(panel_length) or panel_length <= 0:
        raise ValueError("Use a positive, finite diaphragm spacing.")
    shear_buckling = min(5 * PLATE_FACTOR * ((w.width / w.height)**2
                                             + (w.width / panel_length)**2) for w in webs)
    shear = {'shear': Limit(Matboard.s_s.m_as(u.MPa), shear_per_n),
             'glue': _weakest(glue), 'web_shear': Limit(shear_buckling, shear_per_n)}
    positive = _bending_checks(strips, webs, ybar, inertia, 1)
    negative = _bending_checks(strips, webs, ybar, inertia, -1)
    return properties, positive | shear, negative | shear


def _length_at(value, x):
    if callable(value):
        value = value(x * u.mm)
    return float(u.Quantity(value, u.mm).magnitude)


def _section_key(beam, x):
    strips = []
    for plate in beam.matboards:
        if plate.dirc == 'x':
            continue
        points = plate.get_sec_mm('x', x)
        if len(points):
            y0, z0 = points[:, 1:].min(axis=0)
            y1, z1 = points[:, 1:].max(axis=0)
            strips.append(Strip(plate.dirc, y0, y1, z0, z1))
    joints = tuple((_length_at(h, x), _length_at(b, x)) for h, b in beam.glue)
    if beam.diaphragm_spacing is None:
        raise ValueError("Set beam.diaphragm_spacing to calculate web shear buckling.")
    return tuple(strips), joints, _length_at(beam.diaphragm_spacing, x)


def check_section(beam, x, cache=None):
    """Get geometry and (sagging, hogging) mode limits at x in mm."""
    key = _section_key(beam, x)
    if cache is None:
        return _section(*key)
    if key not in cache:
        cache[key] = _section(*key)
    return cache[key]


def analyze(beam, stations, envelope):
    """Return numeric columns for geometry, demands, capacities and all local FOS."""
    sections = {}
    rows = []
    for i, x in enumerate(stations.m_as(u.mm)):
        properties, positive, negative = check_section(beam, x, sections)
        row = {'x_mm': x, **properties}
        row.update({name: record.values[i] for name, record in envelope.items()})
        for mode in MODE_LABELS:
            prefix = 'V' if mode in SHEAR_MODES else 'M'
            plus, minus = positive[mode], negative[mode]
            dp, dm = row[prefix + '_max'], -row[prefix + '_min']
            util_p, util_m = dp / plus.capacity, dm / minus.capacity
            limit, demand = (minus, dm) if util_m > util_p else (plus, dp)
            utilization = max(util_p, util_m)
            row[mode + '_fos'] = 1 / utilization if utilization else inf
            row[mode + '_stress_MPa'] = demand * limit.stress_per_unit
            row[mode + '_limit_MPa'] = limit.strength
            row[mode + '_capacity_pos'] = plus.capacity
            row[mode + '_capacity_neg'] = -minus.capacity
        rows.append(row)
    return {key: np.array([row[key] for row in rows]) for key in rows[0]}
