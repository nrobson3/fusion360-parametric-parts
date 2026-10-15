"""Tests for the parts of each script that do not need Fusion.

Run from the repository root:  python -m unittest discover tests
"""
import importlib.util
import math
import os
import random
import struct
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(rel):
    spec = importlib.util.spec_from_file_location(os.path.basename(rel)[:-3], os.path.join(ROOT, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gear = load('planetary_gearset/planetary_gearset.py')
grip = load('compliant_gripper/compliant_gripper.py')
snap = load('snap_fit_housing/snap_calc.py')
draft = load('verify/draft_check/draft_check.py')
step = load('verify/step_roundtrip/step_roundtrip.py')
stl = load('verify/stl_check/stl_check.py')


class GearSet(unittest.TestCase):
    def test_default_set_is_valid(self):
        errors, info = gear.check_set(18, 18, 54, 3, 1.25, 20.0)
        self.assertEqual(errors, [])
        self.assertAlmostEqual(info['ratio_per_stage'], 4.0)
        self.assertAlmostEqual(info['centre_distance'], 22.5)

    def test_rejections(self):
        self.assertTrue(any('undercut' in e for e in gear.check_set(15, 15, 45, 3, 1.25, 20.0)[0]))
        self.assertTrue(any('assembly' in e for e in gear.check_set(18, 17, 52, 3, 1.25, 20.0)[0]))
        self.assertTrue(any('coaxial' in e for e in gear.check_set(18, 18, 55, 3, 1.25, 20.0)[0]))
        self.assertTrue(any('too close' in e for e in gear.check_set(18, 30, 78, 6, 1.25, 20.0)[0]))

    def test_min_teeth(self):
        self.assertEqual(gear.min_teeth_no_undercut(20.0), 18)   # 2/sin^2(20) = 17.1
        self.assertEqual(gear.min_teeth_no_undercut(25.0), 12)

    def test_reduces_to_hand_derived_form_for_18_18_54(self):
        # P = 10 + 2T - S and P = 3R - 2T - 10 (mod 20), derived by hand for this set
        rnd = random.Random(1)
        for _ in range(200):
            th, s, r = rnd.uniform(0, 360), rnd.uniform(0, 360), rnd.uniform(0, 360)
            self.assertAlmostEqual(gear.planet_phase_external(th, s, 18, 18),
                                   gear.wrap(10 + 2 * th - s, 20), places=9)
            self.assertAlmostEqual(gear.planet_phase_internal(th, r, 18, 54),
                                   gear.wrap(3 * r - 2 * th - 10, 20), places=9)

    def test_solved_stages_mesh_exactly(self):
        rnd = random.Random(2)
        combos = [(ns, np_, ns + 2 * np_, n) for ns in range(17, 31) for np_ in range(17, 31)
                  for n in (3, 4, 5) if (2 * ns + 2 * np_) % n == 0
                  and not gear.check_set(ns, np_, ns + 2 * np_, n, 1.0, 20.0)[0]]
        self.assertGreater(len(combos), 20)
        for ns, np_, nr, n in combos:
            phi_r = rnd.uniform(0, 360)
            for kw in ({'theta': rnd.uniform(0, 360)}, {'phi_s': rnd.uniform(0, 360)}):
                phi_s, planets = gear.solve_stage(ns, np_, nr, n, phi_r, **kw)
                self.assertEqual(len(planets), n)
                for e, i in gear.mesh_errors(ns, np_, nr, phi_s, phi_r, planets):
                    self.assertLess(e, 1e-9)
                    self.assertLess(i, 1e-9)

    def test_helix_rates_keep_mesh_along_the_face(self):
        # if both mesh conditions hold at z=0 they must hold at every z
        ns, np_, nr, m, b = 18, 18, 54, 1.25, 15.0
        ws, wp, wr = (gear.helix_rate(ns, m, b, +1), gear.helix_rate(np_, m, b, -1),
                      gear.helix_rate(nr, m, b, -1))
        phi_s, planets = gear.solve_stage(ns, np_, nr, 3, 0.0, theta=0.0)
        for z in (0.0, 3.3, 8.0, 21.7):
            moved = [(th, pp + wp * z) for th, pp in planets]
            for e, i in gear.mesh_errors(ns, np_, nr, phi_s + ws * z, wr * z, moved):
                self.assertLess(e, 1e-9)
                self.assertLess(i, 1e-9)

    def test_tooth_thickness_at_pitch_circle(self):
        m, n, thin = 1.25, 18, 0.2
        pts = gear.external_profile(m, n, thin, 0.0, 20.0, npts=400)
        rp = m * n / 2
        near = min(pts, key=lambda p: abs(math.hypot(*p) - rp) + (0 if p[1] > 0 else 1e9))
        half = math.atan2(near[1], near[0])
        self.assertAlmostEqual(2 * half * rp, math.pi * m / 2 - thin, delta=0.005)


class Gripper(unittest.TestCase):
    def test_strain_formula(self):
        self.assertAlmostEqual(grip.flexure_strain(0.8, 12.0, 8.0), 0.8 * math.radians(12) / 16)

    def test_default_sizing(self):
        s = grip.size()
        self.assertEqual(s['problems'], [])
        self.assertAlmostEqual(s['strain'], 0.010472, places=5)
        self.assertAlmostEqual(s['materials']['PLA']['safety'], 1.91, places=2)
        self.assertAlmostEqual(s['materials']['PETG']['safety'], 3.82, places=2)
        self.assertAlmostEqual(s['kinematics']['travel_per_jaw'], 5.198, places=3)

    def test_max_thickness_inverts_strain(self):
        t = grip.max_thickness(12.0, 8.0, 0.02)
        self.assertAlmostEqual(grip.flexure_strain(t, 12.0, 8.0), 0.02)

    def test_layout_matches_the_built_part(self):
        lay = grip.layout()
        self.assertEqual(lay['base'], (-42.0, 0.0, 42.0, 12.0))
        self.assertEqual(lay['half']['outer_link'], (28.0, 12.0, 36.0, 45.0))
        self.assertEqual(lay['flexure_y'], (16.0, 41.0))
        v = grip.expected_volumes(lay, 0.8)
        # the as-built Fusion part measured 24.3200 / 20.6336 / 20.6040 cm3
        self.assertAlmostEqual(v['outline'], 24320.0, places=3)
        self.assertAlmostEqual(v['notched'], 20633.6, places=3)
        self.assertAlmostEqual(v['filleted'], 20642.4, delta=0.05)
        self.assertAlmostEqual(v['final'], 20604.0, delta=0.05)


class SnapFit(unittest.TestCase):
    def test_published_numbers(self):
        out = snap.calc()
        self.assertAlmostEqual(out['allowable']['repeated'], 0.01092, places=5)
        self.assertAlmostEqual(out['assembly_nominal']['strain'], 0.007377, places=6)
        self.assertAlmostEqual(out['push_on']['two_arms_N'], 13.92, places=2)
        self.assertAlmostEqual(out['retention']['barb_shear_N_per_arm'], 500.76, places=2)

    def test_stop_caps_strain(self):
        self.assertLess(snap.strain(snap.stop_gap), snap.eps_once)


class Draft(unittest.TestCase):
    def test_bands(self):
        self.assertEqual(draft.classify([1.5, 2.0]), 'plus')
        self.assertEqual(draft.classify([-1.5, -0.9]), 'minus')
        self.assertEqual(draft.classify([90.0, 89.9]), 'flat')
        self.assertEqual(draft.classify([0.0, 0.2]), 'low')
        self.assertEqual(draft.classify([-30.0, 30.0]), 'split')

    def test_angle(self):
        self.assertAlmostEqual(draft.draft_deg((math.cos(math.radians(1.5)), 0, math.sin(math.radians(1.5))),
                                               (0, 0, 1)), 1.5)


STEP_MM = """ISO-10303-21;
HEADER;
ENDSEC;
DATA;
#1=PRODUCT('PB-100','PB-100','',(#2));
#3=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );
#4=( CONVERSION_BASED_UNIT('DEGREE',#5) NAMED_UNIT(*) PLANE_ANGLE_UNIT() );
#6=MANIFOLD_SOLID_BREP('',#7);
#8=BREP_WITH_VOIDS('',#9,(#10));
ENDSEC;
END-ISO-10303-21;
"""
STEP_INCH = STEP_MM.replace(
    "#4=", "#11=( CONVERSION_BASED_UNIT('inch',#12) LENGTH_UNIT() NAMED_UNIT(#13) );\n#4=").replace(
    "#6=MANIFOLD_SOLID_BREP('',#7);\n", "#6=OPEN_SHELL('',(#7));\n#14=PRODUCT('2026-01-01-10-00-00-000','x','',(#2));\n")


class Step(unittest.TestCase):
    def test_mm_file(self):
        info = step.inspect_step(STEP_MM)
        self.assertTrue(info['header_ok'] and info['terminator_ok'])
        self.assertEqual(info['length_unit'], 'millimetre')   # not fooled by DEGREE
        self.assertEqual(info['solids'], 2)
        self.assertEqual(step.file_problems(info, 2), [])
        self.assertTrue(step.file_problems(info, 3))

    def test_inch_and_surfaces(self):
        info = step.inspect_step(STEP_INCH)
        self.assertEqual(info['length_unit'], 'inch')
        self.assertTrue(any('surface' in p for p in step.file_problems(info, 1)))

    def test_truncated(self):
        self.assertTrue(step.file_problems(step.inspect_step(STEP_MM[:-20]), 2))

    def test_compare(self):
        a = {'bodies': 2, 'volume_mm3': 1000.0, 'bbox_mm': (0, 0, 0, 10, 10, 10)}
        self.assertEqual(step.compare(a, dict(a, volume_mm3=1000.001)), [])
        self.assertTrue(step.compare(a, dict(a, bodies=1)))
        self.assertTrue(step.compare(a, dict(a, bbox_mm=(0, 0, 0, 10, 10, 10.5))))


def cube(s=10.0, o=(0.0, 0.0, 0.0)):
    x, y, z = o
    v = [(x + i * s, y + j * s, z + k * s) for i in (0, 1) for j in (0, 1) for k in (0, 1)]
    q = [(0, 1, 3, 2), (4, 6, 7, 5), (0, 4, 5, 1), (2, 3, 7, 6), (0, 2, 6, 4), (1, 5, 7, 3)]
    tris = []
    for a, b, c, d in q:
        tris += [(v[a], v[b], v[c]), (v[a], v[c], v[d])]
    return tris


class Stl(unittest.TestCase):
    def test_closed_cube(self):
        r = stl.check_mesh(cube())
        self.assertTrue(r['watertight'])
        self.assertEqual((r['shells'], r['euler']), (1, 2))
        self.assertAlmostEqual(r['volume'], 1000.0)
        self.assertEqual(stl.problems(r, 1000.0), [])

    def test_hole(self):
        r = stl.check_mesh(cube()[:-1])
        self.assertEqual(r['boundary_edges'], 3)
        self.assertFalse(r['watertight'])

    def test_flipped_triangle(self):
        t = cube()
        a, b, c = t[0]
        t[0] = (a, c, b)
        r = stl.check_mesh(t)
        self.assertGreater(r['inconsistent_winding'], 0)

    def test_inside_out_and_loose_shell(self):
        inv = [(a, c, b) for a, b, c in cube()]
        self.assertTrue(any('inward' in p for p in stl.problems(stl.check_mesh(inv))))
        self.assertEqual(stl.check_mesh(cube() + cube(o=(20, 0, 0)))['shells'], 2)

    def test_binary_roundtrip(self):
        path = os.path.join(tempfile.mkdtemp(), 'c.stl')
        stl.write_binary_stl(path, cube())
        self.assertEqual(len(stl.read_stl(path)), 12)
        self.assertEqual(os.path.getsize(path), 84 + 12 * 50)
        self.assertAlmostEqual(stl.check_mesh(stl.read_stl(path))['volume'], 1000.0, places=3)


if __name__ == '__main__':
    unittest.main()
