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

    def _section(self, dirc, val):
        axis = d[dirc]
        normal_axis = d[self.dirc]
        vertices = u.Quantity(self.vertices, u.mm)
        position = u.Quantity(val, u.mm)
        half_thickness = self.t / 2

        if dirc == self.dirc:
            board_pos = vertices[0, axis]
            if board_pos - half_thickness <= position <= board_pos + half_thickness:
                section = vertices.copy()
                section[:, axis] = position
                return section
            return np.empty((0, 3)) * u.mm

        descending = ascending = None
        for start, end in zip(vertices, np.roll(vertices, -1, axis=0)):
            low, high = sorted((start[axis], end[axis]))
            if low == high or not low <= position <= high:
                continue

            fraction = (position - start[axis]) / (end[axis] - start[axis])
            point = (1 - fraction) * start + fraction * end
            point[axis] = position
            if start[axis] > end[axis]:
                descending = point
            else:
                ascending = point

        if descending is None or ascending is None or (descending == ascending).all():
            return np.empty((0, 3)) * u.mm

        offset = np.zeros(3) * u.mm
        offset[normal_axis] = half_thickness
        return np.stack([
            descending - offset, descending + offset,
            ascending + offset, ascending - offset,
        ])

    def get_sec_mm(self, dirc, val):
        return self._section(dirc, val).magnitude

    def get_sec(self, dirc, val):
        points = self._section(dirc, val)
        if dirc == self.dirc:
            return [tuple(point) for point in points]
        return [list(point) for point in points]

    def get_sec_area(self, dirc, val):
        return compute.polygon_area(self.get_sec_mm(dirc, val))

    def get_sec_center(self, dirc, val):
        return compute.geometric_centroid(self.get_sec_mm(dirc, val))
