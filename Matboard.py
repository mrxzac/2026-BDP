import numpy as np

import compute

u = compute.u
d = compute.d


class Matboard:
    s_t = 30 * u.MPa 
    """ tensile stress """
    s_c = 6 * u.MPa     
    """ compressive stress """
    s_s = 4 * u.MPa     
    """ shear stress """
    E_m = 4000 * u.MPa  
    """ Young's Modulus """
    v_m = 0.2           
    """ Poisson's Ratio """
    t = 1.27 * u.mm   
    """ thickness """

    def __init__(self, vertices, dirc):
        self.dirc = dirc
        self.vertices = compute.points_in_mm(vertices)

    def get_sec_mm(self, dirc, val):
        """Cut a flat plate with a single connected intersection at the given plane."""
        axis = d[dirc]
        normal = d[self.dirc]
        position = u.Quantity(val, u.mm).magnitude
        vertices = self.vertices
        half_t = self.t.m_as(u.mm) / 2
        empty = np.empty((0, 3))

        if axis == normal:
            centre = vertices[0, axis]
            if not centre - half_t <= position <= centre + half_t:
                return empty
            section = vertices.copy()
            section[:, axis] = position
            return section

        if position < vertices[:, axis].min() or position > vertices[:, axis].max():
            return empty
        ends = np.roll(vertices, -1, axis=0)
        start_coord, end_coord = vertices[:, axis], ends[:, axis]
        hits = ((start_coord != end_coord)
                & (position >= np.minimum(start_coord, end_coord))
                & (position <= np.maximum(start_coord, end_coord)))
        starts, ends = vertices[hits], ends[hits]
        if not len(starts):
            return empty
        fraction = (position - starts[:, axis]) / (ends[:, axis] - starts[:, axis])
        points = starts + fraction[:, None] * (ends - starts)
        tangent = 3 - axis - normal
        crossings = np.unique(points[:, tangent])
        if len(crossings) < 2:
            return empty
        if len(crossings) > 2:
            # A cut along notch shoulders can touch several vertices but still
            # be one strip. Check each interval for polygon interior or boundary.
            mid = ((crossings[:-1] + crossings[1:]) / 2)[:, None]
            straddles = (starts[:, axis] > position) != (ends[:, axis] > position)
            interior = np.count_nonzero(points[straddles, tangent] > mid, axis=1) % 2
            flat = (start_coord == position) & (end_coord == position)
            t0 = vertices[:, tangent][flat]
            t1 = np.roll(vertices[:, tangent], -1)[flat]
            boundary = ((np.minimum(t0, t1) <= mid) & (mid <= np.maximum(t0, t1))).any(axis=1)
            if not np.all(interior | boundary):
                raise ValueError("Split a plate with disjoint cross-sections into separate matboards.")
        low, high = points[:, tangent].argmin(), points[:, tangent].argmax()
        if starts[low, axis] < ends[low, axis]:
            low, high = high, low
        descending, ascending = points[low], points[high]
        offset = np.zeros(3)
        offset[normal] = half_t
        return np.stack([descending - offset, descending + offset,
                         ascending + offset, ascending - offset])

    def _section(self, dirc, val):
        return self.get_sec_mm(dirc, val) * u.mm

    def get_sec(self, dirc, val):
        points = self._section(dirc, val)
        if dirc == self.dirc:
            return [tuple(point) for point in points]
        return [list(point) for point in points]

    def get_sec_area(self, dirc, val):
        return compute.polygon_area(self.get_sec_mm(dirc, val))

    def get_sec_center(self, dirc, val):
        return compute.geometric_centroid(self.get_sec_mm(dirc, val))
