"""Regression checks against point-load equilibrium and Design 0 hand calculations."""

import unittest

import numpy as np

from Beam import Beam
from config_0 import beam_0, u
from Matboard import Matboard
from Simulation import Simulation
from strength import analyze, check_section


class BridgeTests(unittest.TestCase):
    def test_repeated_off_grid_loads_on_translated_bridge(self):
        shift = 137.5
        plates = [Matboard(p.vertices + [shift, 0, 0], p.dirc) for p in beam_0.matboards]
        beam = Beam(plates, beam_0.glue, beam_0.diaphragm_spacing)
        sim = Simulation([120, 0, 0, 0, 0, 0] * u.N, beam, density=4)
        for distance in (-25., 237.25, 1225., 1225.01):
            sim.loads(distance * u.mm)
            total = 120. if distance <= 1225. else 0.
            right = total * distance / 1200
            np.testing.assert_allclose([r.m_as(u.N) for r in sim.reaction()], [total - right, right])
            self.assertAlmostEqual(sim.load_vector.magnitude.sum(), total)
            x, _, _, m = sim.diagram(distance)
            self.assertAlmostEqual(m[0], 0)
            self.assertAlmostEqual(m[-1], 0)
            if distance == 237.25:
                self.assertAlmostEqual(m[x == sim.left + distance].item(), (total - right) * distance)

    def test_centered_base_train(self):
        sim = Simulation([67.5, 67.5, 67.5, 67.5, 91., 91.] * u.N, beam_0)
        sim.loads(172 * u.mm)
        np.testing.assert_allclose([r.m_as(u.N) for r in sim.reaction()], [212.683333333, 239.316666667])
        # The shear changes sign at the fourth axle, 688 mm from the left support.
        peak = 212.683333333333 * 688 - 67.5 * (516 + 340 + 176)
        self.assertAlmostEqual(sim.BMD().magnitude.max(), peak)
        self.assertAlmostEqual(sim.SFD().magnitude[-1], 0)

    def test_single_axle_signed_envelope_including_overhangs(self):
        sim = Simulation([120, 0, 0, 0, 0, 0] * u.N, beam_0, density=51)
        env = sim.envelopes()
        s = sim.x.magnitude - sim.left
        inside = (s >= 0) & (s <= 1200)
        np.testing.assert_allclose(env['M_max'].values[inside],
                                   120 * s[inside] * (1200 - s[inside]) / 1200, atol=1e-8)
        # Load at a free end: support moment and negative-moment influence over span.
        np.testing.assert_allclose(env['M_min'].values[inside],
                                   -120 * 25 * np.maximum(s[inside], 1200 - s[inside]) / 1200)
        self.assertAlmostEqual(env['V_max'].values.max(), 120.)
        self.assertAlmostEqual(env['V_min'].values.min(), -120.)

    def test_local_reinforcement_only_contributes_where_present(self):
        patch = Matboard(np.array([[.5, .08, -.05], [.7, .08, -.05],
                                   [.7, .08, .05], [.5, .08, .05]]) * u.m, 'y')
        beam = Beam([*beam_0.matboards, patch], beam_0.glue, beam_0.diaphragm_spacing)
        self.assertEqual(patch.get_sec_area('x', 400 * u.mm).magnitude, 0)
        self.assertAlmostEqual(patch.get_sec_area('x', 600 * u.mm).magnitude, 127)
        reference = beam_0.get_2nd_moment('z', 0 * u.mm).magnitude
        self.assertAlmostEqual(beam.get_2nd_moment('z', 400 * u.mm).magnitude, reference)
        self.assertGreater(beam.get_2nd_moment('z', 600 * u.mm).magnitude, reference)

    def test_diaphragms_fill_the_box_and_clear_the_tabs(self):
        plates = [p for p in beam_0.matboards if p.dirc == 'x']
        self.assertEqual([p.vertices[0, 0] for p in plates], [10, 40, 425, 825, 1210, 1240])
        t = 1.27
        clear_width, clear_height = 80 - 2*t, 75 - t
        plate_area = clear_width * clear_height - 2*5*t
        for plate in plates:
            x = plate.vertices[0, 0]
            for cut in [x - t/2, x, x + t/2]:
                self.assertAlmostEqual(plate.get_sec_area('x', cut).magnitude, plate_area)
            self.assertEqual(plate.get_sec_area('x', x + t/2 + .001).magnitude, 0)
            total = sum(p.get_sec_area('x', x).magnitude for p in beam_0.matboards)
            self.assertAlmostEqual(total, 80*75 + 100*t)
        plate = plates[0]
        for y, width in [(75-t, clear_width), (75-t/2, clear_width-10)]:
            self.assertAlmostEqual(plate.get_sec_area('y', y).magnitude, t*width)
        self.assertAlmostEqual(plate.get_sec_area('z', clear_width/2-5).magnitude, t*clear_height)

    def test_web_buckling_uses_the_local_diaphragm_bay(self):
        for x, spacing in [(0, np.inf), (25, 30), (41, 385), (425, 400),
                           (625, 400), (900, 385), (1225, 30), (1250, np.inf)]:
            self.assertEqual(beam_0.web_panel_length(x).m_as(u.mm), spacing)
        # Removing the interior brace at 425 merges the 385 and 400 mm bays.
        boards = [p for p in beam_0.matboards if p.dirc != 'x' or p.vertices[0, 0] != 425]
        longer_bay = Beam(boards, beam_0.glue)
        original, positive, _ = check_section(beam_0, 425)
        modified, new_positive, _ = check_section(longer_bay, 425)
        self.assertEqual(modified['panel_length_mm'], 785)
        self.assertLess(new_positive['web_shear'].capacity, positive['web_shear'].capacity)
        self.assertEqual(original['Iz_mm4'], modified['Iz_mm4'])

    def test_design_zero_governing_mode(self):
        sim = Simulation(np.full(6, 400 / 6) * u.N, beam_0)
        env = sim.envelopes()
        data = analyze(beam_0, sim.x, env)
        # Independent rectangle decomposition (width, height, centroid elevation).
        t = 1.27
        layers = [(100, t, 75 + t/2), (80, t, t/2), (2*t, 75-t, (75+t)/2), (10, t, 75-t/2)]
        y = sum(b*h*c for b, h, c in layers) / sum(b*h for b, h, c in layers)
        inertia = sum(b*h**3/12 + b*h*(c-y)**2 for b, h, c in layers)
        q = 100*t*(75+t/2-y) + 10*t*(75-t/2-y) + t*(75-y)**2
        np.testing.assert_allclose(data['ybar_mm'], y)
        np.testing.assert_allclose(data['Iz_mm4'], inertia)
        self.assertAlmostEqual(beam_0.get_2nd_moment('z', 0).magnitude, inertia)
        np.testing.assert_allclose(data['Q_NA_mm3'], q)
        np.testing.assert_allclose(data['Q_glue_1_mm3'], 100*t*(75+t/2-y))
        # Axle 3 and the train resultant straddle midspan equally (556 and 644 mm).
        p = 400 / 6
        m = (400 * (1200-644)/1200) * 556 - p * (340 + 164)
        self.assertAlmostEqual(env['M_max'].values.max(), m)
        self.assertAlmostEqual(env['V_max'].values.max(), 400 * (1200-428)/1200)
        sigma_cr = 4*np.pi**2*4000/(12*(1-0.2**2)) * (t/(80-t))**2
        expected = sigma_cr * inertia / ((75+t-y) * m)
        minimum = min(np.min(v) for k, v in data.items() if k.endswith('_fos'))
        self.assertAlmostEqual(minimum, expected)
        self.assertAlmostEqual(data['middle_flange_fos'].min(), expected)


if __name__ == '__main__':
    unittest.main()
