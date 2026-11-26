from Matboard import u
from compute import geometric_centroid, polygon_area


def _centroid_coordinate(sections, dirc):
    index = {"x": 0, "y": 1, "z": 2}[dirc]
    total_area = sum(area for _, area, _ in sections)
    return sum(area * center[index] for _, area, center in sections) / total_area


def _local_moments(points, center):
    yz = [((y - center[1]).m_as("mm"), (z - center[2]).m_as("mm"))
          for _, y, z in points]
    iy = iz = 0.0
    for (y0, z0), (y1, z1) in zip(yz, yz[1:] + yz[:1]):
        cross = y0 * z1 - y1 * z0 / 2
        iy += cross * (z0**2 + z0 * z1 + z1**2)
        iz += cross * (y0**2 + y0 * y1 + y1**2)
    return abs(iy) / 6 * u.mm**4, abs(iz) / 6 * u.mm**4


class Beam:
    def __init__(self, matboards):
        self.matboards = matboards

    def _sections_at(self, x):
        sections = []
        for board in self.matboards:
            points = board.get_sec("x", x)
            area = polygon_area(points)
            if area:
                sections.append((points, area, geometric_centroid(points)))
        if not sections:
            raise ValueError(f"No cross-sectional area at x={x}.")
        return sections

    def get_centroid(self, dirc, x):
        """Coordinate of the area centroid of the cross-section at x."""
        return _centroid_coordinate(self._sections_at(x), dirc)

    def get_2nd_moment(self, dirc, x):
        """Area moment about centroidal axis dirc (x gives the polar area moment)."""
        sections = self._sections_at(x)
        y_bar = _centroid_coordinate(sections, "y")
        z_bar = _centroid_coordinate(sections, "z")
        iy = iz = 0 * u.mm**4
        for points, area, center in sections:
            local_iy, local_iz = _local_moments(points, center)
            iy += local_iy + area * (center[2] - z_bar)**2
            iz += local_iz + area * (center[1] - y_bar)**2
        return {"x": iy + iz, "y": iy, "z": iz}[dirc]
