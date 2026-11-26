"""Statics of a six-axle train on a simply supported bridge with end overhangs."""

from dataclasses import dataclass

import numpy as np

from compute import u


@dataclass
class Extreme:
    values: np.ndarray
    first_axle_mm: np.ndarray  # Relative to the left support
    reversed: np.ndarray      # True when the supplied axle weights were reversed
    approach: np.ndarray      # -1 / +1 denote the left / right train-position limits


class Simulation:
    axle_positions = np.array([0, 176, 340, 516, 680, 856]) * u.mm
    span = 1200 * u.mm

    def __init__(self, loadcase, beam, density=1251):
        self.loadcase = u.Quantity(loadcase, u.N)
        weights = self.loadcase.magnitude
        if (weights.shape != (6,) or not np.isfinite(weights).all()
                or (weights < 0).any() or not weights.sum() > 0):
            raise ValueError("Supply six finite, nonnegative axle forces with positive total weight.")
        if not isinstance(density, (int, np.integer)) or density < 2:
            raise ValueError("density must be an integer of at least two stations.")
        self.beam = beam
        self.density = density
        start, end = beam.get_endpoints()
        self.length = end - start
        self.side = (self.length - self.span) / 2
        if self.side.magnitude < 0:
            raise ValueError("The bridge must cover the 1200 mm support span.")
        self.start, self.end = start.m_as(u.mm), end.m_as(u.mm)
        self.left, self.right = (start + self.side).m_as(u.mm), (end - self.side).m_as(u.mm)
        self.x = np.unique(np.r_[np.linspace(self.start, self.end, density),
                                 self.left, self.right]) * u.mm

    def _active(self, axles, weights, approach):
        # Round coordinate differences only, so coincident metric inputs compare equally.
        at_start = np.round(axles - self.start, 9)
        at_end = np.round(axles - self.end, 9)
        entered = at_start >= 0 if approach >= 0 else at_start > 0
        remaining = at_end <= 0 if approach <= 0 else at_end < 0
        return np.where(entered & remaining, weights, 0.)

    def _reactions(self, axles, active):
        right = (active * (axles - self.left)).sum(axis=-1) / (self.right - self.left)
        return active.sum(axis=-1) - right, right

    def _response(self, x, first, weights, approach=0):
        """V on each face of a cut and sagging-positive M, in N and N mm."""
        axles = np.asarray(first)[..., None] + self.axle_positions.magnitude
        active = self._active(axles, weights, approach)
        left, right = self._reactions(axles, active)
        arm = np.round(x[:, None] - axles, 9)
        before = arm > 0 if approach >= 0 else arm >= 0
        through = arm >= 0 if approach <= 0 else arm > 0
        v_left = left * (x > self.left) + right * (x > self.right) - (active * before).sum(axis=-1)
        v_right = left * (x >= self.left) + right * (x >= self.right) - (active * through).sum(axis=-1)
        moment = (left * np.maximum(x - self.left, 0) + right * np.maximum(x - self.right, 0)
                  - (active * np.maximum(arm, 0)).sum(axis=-1))
        moment[np.abs(moment) < 1e-8] = 0
        return v_left, v_right, moment

    def loads(self, position):
        """Place the leftmost axle this distance from the left support (positive right)."""
        self.first = self.left + u.Quantity(position, u.mm).magnitude
        axles = self.first + self.axle_positions.magnitude
        active = self._active(axles, self.loadcase.magnitude, 0)
        self.reac_L, self.reac_R = (r * u.N for r in self._reactions(axles, active))
        # The grid vector is for inspection. Statics use the true axle coordinates.
        self.load_vector = np.zeros(len(self.x)) * u.N
        nearest = np.abs(self.x.magnitude[:, None] - axles).argmin(axis=0)
        np.add.at(self.load_vector.magnitude, nearest, active)
        _, shear, moment = self._response(self.x.magnitude, self.first, self.loadcase.magnitude)
        self.sfd, self.bmd = shear * u.N, moment * u.N * u.mm
        return self.load_vector

    def reaction(self):
        return self.reac_L, self.reac_R

    def SFD(self):
        return self.sfd

    def BMD(self):
        return self.bmd

    def envelopes(self):
        """Extrema over continuous train travel in both orientations at each station.

        With the cut fixed, the influence functions are piecewise linear in train
        position. Breaks occur when an axle enters, exits, or crosses that cut.
        Evaluate the event and both one-sided limits to retain shear jumps.
        """
        x = self.x.magnitude
        unloaded = self.start - self.left - self.axle_positions.magnitude[-1] - 1
        envelope = {name: Extreme(np.zeros_like(x), np.full_like(x, unloaded),
                                  np.zeros(x.size, dtype=bool), np.zeros(x.size, dtype=int))
                    for name in ('V_min', 'V_max', 'M_min', 'M_max')}
        weights = self.loadcase.magnitude
        orders = [weights] if np.array_equal(weights, weights[::-1]) else [weights, weights[::-1]]
        for order, forces in enumerate(orders):
            for offset in self.axle_positions.magnitude:
                for event in (self.start, self.end, x):
                    first = event - offset
                    for approach in (0, -1, 1):
                        vl, vr, m = self._response(x, first, forces, approach)
                        candidates = (np.minimum(vl, vr), np.maximum(vl, vr), m, m)
                        for (name, record), value in zip(envelope.items(), candidates):
                            better = value < record.values if name.endswith('min') else value > record.values
                            np.copyto(record.values, value, where=better)
                            np.copyto(record.first_axle_mm, first - self.left, where=better)
                            record.reversed[better] = bool(order)
                            record.approach[better] = approach
        return envelope

    def diagram(self, position, reverse=False, approach=0):
        """Return exact event coordinates, shear on both faces, and bending moment."""
        first = self.left + u.Quantity(position, u.mm).magnitude
        if approach:
            # Separate coincident wheel and support jumps for a limiting-position plot.
            first += approach * 1e-6
        axles = first + self.axle_positions.magnitude
        events = np.unique(np.r_[self.start, self.left, self.right, self.end,
                                 axles[(axles >= self.start) & (axles <= self.end)]])
        weights = self.loadcase.magnitude
        return (events, *self._response(events, first, weights[::-1] if reverse else weights))

    def plot(self, show=True):
        import matplotlib.pyplot as plt

        x, vl, vr, m = self.diagram(self.first - self.left)
        fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 6), layout='constrained')
        axes[0].plot(np.repeat(x, 2), np.column_stack((vl, vr)).ravel())
        axes[0].set_ylabel('Shear (N)')
        axes[1].plot(x, m / 1000)
        axes[1].set_ylabel('Moment (N·m)')
        axes[1].set_xlabel('Bridge x (mm)')
        for ax in axes:
            ax.axhline(0, color='0.3', linewidth=0.7)
            ax.grid(alpha=0.3)
        fig.suptitle(f'Train with first axle {self.first - self.left:g} mm from left support')
        if show:
            plt.show()
        return fig
