import compute

class Beam:
    def __init__(self, matboards, glue, diaphragm_spacing=None):
        self.matboards = matboards
        self.glue = glue  # (interface elevation, total bonded width) per joint
        self.diaphragm_spacing = diaphragm_spacing

    def web_panel_length(self, x):
        """Centre-to-centre brace spacing in mm; at a shared brace use the longer bay."""
        position = compute.u.Quantity(x, compute.u.mm)
        centres = sorted({b.vertices[0, 0] for b in self.matboards if b.dirc == "x"})
        if centres:
            bays = [b - a for a, b in zip(centres, centres[1:])
                    if a <= position.magnitude <= b]
            # At a free end there is no second diaphragm to brace a finite bay.
            return max(bays, default=float("inf")) * compute.u.mm
        if self.diaphragm_spacing is None:
            return float("inf") * compute.u.mm
        spacing = self.diaphragm_spacing
        return compute.u.Quantity(spacing(position) if callable(spacing) else spacing, compute.u.mm)

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

    def _sections_at(self, x):
        sections = []
        for board in self.matboards:
            points = board.get_sec_mm("x", x)
            area, center = compute.polygon_properties_mm(points)
            if area:
                sections.append((points, area, center))
        if not sections:
            raise ValueError(f"No cross-sectional area at x={x}.")
        return sections

    def get_centroid(self, dirc, x):
        center = compute._centroid(self._sections_at(x))
        return center[compute.d[dirc]] * compute.u.mm

    def get_2nd_moment(self, dirc, x):
        sections = self._sections_at(x)
        _, y_bar, z_bar = compute._centroid(sections)
        iy = iz = 0.0
        for points, area, center in sections:
            local_iy, local_iz = compute._local_moments(points, center)
            iy += local_iy + area * (center[2] - z_bar)**2
            iz += local_iz + area * (center[1] - y_bar)**2
        moment = {"x": iy + iz, "y": iy, "z": iz}[dirc]
        return moment * compute.u.mm**4
