import numpy as np

import compute

class Beam:
    def __init__(self, matboards, glue, buckling_dimensions ):
        self.matboards = matboards
        self.glue = glue
        self.buckling_dimensions = buckling_dimensions 

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
    
    def get_width(self, x, y):
        level = compute.u.Quantity(y, compute.u.mm).magnitude
        intervals = []
        for board in self.matboards:
            points = board.get_sec_mm("x", x)
            if len(points) < 3:
                continue
            intersections = []
            for start, end in zip(points, np.roll(points, -1, axis=0)):
                y1, z1 = start[1:]
                y2, z2 = end[1:]
                if y1 == y2:
                    if level == y1:
                        intervals.append(sorted((z1, z2)))
                    continue
                if min(y1, y2) <= level < max(y1, y2):
                    fraction = (level - y1) / (y2 - y1)
                    z = z1 + fraction * (z2 - z1)
                    intersections.append(z)
            intersections.sort()
            for i in range(0, len(intersections) - 1, 2):
                intervals.append(intersections[i:i + 2])
                
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])

        width = sum(end - start for start, end in merged)

        return width * compute.u.mm

    def get_1st_moment(self, x, y, dirc='y'):
        sections = self._sections_at(x)
        axis = compute.d[dirc]
        centroid = compute._centroid(sections)[axis]
        Q = 0.0
        for board in self.matboards:
            points = board.get_sec_below_mm(x, y, dirc)
            if len(points) < 3:
                continue
            area, center = compute.polygon_properties_mm(points)
            Q += area * (centroid - center[axis])
        return abs(Q) * compute.u.mm**3
            
