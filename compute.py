import numpy as np
from pint import UnitRegistry

u = UnitRegistry()
d = {"x": 0, "y": 1, "z": 2}


def points_in_mm(points):
    if isinstance(points, np.ndarray) and points.dtype != object:
        return np.asarray(points, dtype=float)
    return np.array([
        [u.Quantity(coord, u.mm).magnitude for coord in point]
        for point in points
    ], dtype=float)


def polygon_properties_mm(points):
    offsets = points[1:] - points[0]
    crosses = np.cross(offsets[:-1], offsets[1:])
    normal = crosses.sum(axis=0)
    area = np.linalg.norm(normal) / 2
    weights = crosses[:, np.argmax(np.abs(normal))]
    triangle_centers = (offsets[:-1] + offsets[1:]) / 3
    center = points[0] + np.average(triangle_centers, axis=0, weights=weights)
    return area, center


def geometric_centroid(points):
    points = points_in_mm(points)
    _, center = polygon_properties_mm(points)
    return tuple(center * u.mm)


def polygon_area(points):
    area, _ = polygon_properties_mm(points_in_mm(points))
    return area * u.mm**2

def _centroid(sections):
    _, areas, centers = zip(*sections)
    return np.average(centers, axis=0, weights=areas)


def _local_moments(points, center):
    y = points[:, 1] - center[1]
    z = points[:, 2] - center[2]
    next_y = np.roll(y, -1)
    next_z = np.roll(z, -1)

    cross = (y * next_z - next_y * z )/2
    iy = np.sum(cross * (z**2 + z * next_z + next_z**2)) / 6
    iz = np.sum(cross * (y**2 + y * next_y + next_y**2)) / 6
    return abs(iy), abs(iz)


