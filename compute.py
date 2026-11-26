import math

from pint import Quantity


def _points_in_mm(points):
    # Matboard imports compute, so use its registry after module initialization.
    from Matboard import u

    coords = [
        tuple(v.m_as("mm") if isinstance(v, Quantity) else v for v in point)
        for point in points
    ]
    return coords, u.mm


def _triangle_cross_products(points):
    """Cross products for the triangle fan about the first vertex, in mm²."""
    x0, y0, z0 = points[0]
    offsets = [(x - x0, y - y0, z - z0) for x, y, z in points[1:]]
    return [
        (ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx)
        for (ax, ay, az), (bx, by, bz) in zip(offsets, offsets[1:])
    ]


def geometric_centroid(points):
    """Return the (x, y, z) area centroid of an ordered planar polygon in mm."""
    points, mm = _points_in_mm(points)
    if len(points) < 3:
        raise ValueError("A polygon requires at least 3 points.")

    crosses = _triangle_cross_products(points)
    normal = tuple(math.fsum(components) for components in zip(*crosses))
    axis = max(range(3), key=lambda i: abs(normal[i]))
    if normal[axis] == 0:
        raise ValueError("Polygon has zero area.")

    # Signed projected areas give each triangle its share of the polygon.
    origin = points[0]
    moments = [0.0, 0.0, 0.0]
    for p1, p2, cross in zip(points[1:], points[2:], crosses):
        for i in range(3):
            tri_offset = ((p1[i] - origin[i]) + (p2[i] - origin[i])) / 3
            moments[i] += cross[axis] * tri_offset

    return tuple((p + moment / normal[axis]) * mm for p, moment in zip(origin, moments))


def polygon_area(points):
    """Return the area of an ordered planar polygon in mm²."""
    points, mm = _points_in_mm(points)
    if len(points) < 3:
        return 0.0 * mm**2

    crosses = _triangle_cross_products(points)
    normal = tuple(math.fsum(components) for components in zip(*crosses))
    return (math.hypot(*normal) / 2) * mm**2