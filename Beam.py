import numpy as np

import compute

class Beam:
    def __init__(self, matboards, glue):
        self.matboards = matboards
        self.glue = glue

    def get_endpoints(self):
        bounds = []
        for board in self.matboards:
            x = board.vertices[:, 0]
            start = x.min() * compute.u.mm
            end = x.max() * compute.u.mm
            if board.dirc == "x":
                start -= board.t / 2
                end += board.t / 2
            bounds.append(start)
            bounds.append(end)
        return min(bounds), max(bounds)
    
    def get_y_bounds(self, x):
        bounds = []
        for board in self.matboards:
            points = board.get_sec_mm("x", x)
            if len(points) < 3:
                continue
            y = points[:, 1]
            bounds.extend([y.min(), y.max()])
        return min(bounds) * compute.u.mm, max(bounds) * compute.u.mm

    def _sections_at(self, x):
        sections = []
        for board in self.matboards:
            points = board.get_sec_mm("x", x)
            if len(points) < 3:
                continue
            area, center = compute.polygon_properties_mm(points)
            if area:
                sections.append((points, area, center))
        if not sections:
            raise ValueError(f"No cross-sectional area at x={x}.")
        return sections

    def get_centroid(self, x, dirc = 'y'):
        center = compute._centroid(self._sections_at(x))
        return center[compute.d[dirc]] * compute.u.mm

    def get_2nd_moment(self, x, dirc = 'y'):
        sections = self._sections_at(x)
        _, y_bar, z_bar = compute._centroid(sections)
        iy = iz = 0.0
        for points, area, center in sections:
            local_iy, local_iz = compute._local_moments(points, center)
            iy += local_iy + area * (center[2] - z_bar)**2
            iz += local_iz + area * (center[1] - y_bar)**2
        moment = {"x": iy + iz, "z": iy, "y": iz}[dirc]
        return moment * compute.u.mm**4
