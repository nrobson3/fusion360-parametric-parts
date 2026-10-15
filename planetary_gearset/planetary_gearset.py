"""Helical planetary gear set generator for Autodesk Fusion.

Builds the sun, planets and ring of a one- or two-stage planetary set from tooth
counts, checks that the set can actually be built (assembly condition, undercut,
planet clearance), and places every planet with its teeth phased into mesh.
With two stages both share one ring, which is the case that needs the coupled
phasing solution in solve_stage().

Run it from Utilities > Add-Ins > Scripts and Add-Ins. It opens a new design;
results go to the Text Commands window and a summary dialog.
"""
import math
import traceback

try:
    import adsk.core
    import adsk.fusion
except ImportError:  # lets the gear math be imported and tested outside Fusion
    adsk = None

# ------------------------------------------------------------------ parameters
MODULE = 1.25            # mm, transverse module
PRESSURE_ANGLE = 20.0    # deg, transverse
HELIX_ANGLE = 15.0       # deg at the pitch cylinder; 0 gives spur gears
N_SUN = 18
N_PLANET = 18
N_RING = 54
N_PLANETS = 3            # planets per stage
STAGES = 2               # 1 or 2; stage 2 sits above stage 1 inside the same ring
FACE_WIDTH = 8.0         # mm, sun and planet face width
STAGE_GAP = 4.0          # mm, axial gap between stages
TOOTH_THINNING = 0.2     # mm of arc thickness taken off every tooth at the pitch circle;
                         # two thinned gears in mesh give 2x this as circumferential backlash
RING_OD = 80.0           # mm
SUN_BORE = 5.0           # mm
PLANET_BORE = 5.4        # mm, 5 mm pin + 0.2 radial clearance
LOFT_SECTIONS = 5        # profiles per helical tooth; 4-6 is plenty
ALLOW_UNDERCUT = False   # True builds below the minimum tooth count anyway
MIN_PLANET_TIP_GAP = 1.0  # mm, planet tip circle to neighbouring planet tip circle


# ------------------------------------------------------------------ gear math
def involute(a):
    return math.tan(a) - a


def wrap(v, p):
    """v modulo p, in [0, p)."""
    r = math.fmod(v, p)
    if r < 0:
        r += p
    return 0.0 if r >= p else r   # -1e-17 + p rounds to p


def min_teeth_no_undercut(pressure_angle_deg):
    """Smallest external tooth count a standard rack cutter does not undercut."""
    # round first: 2/sin^2(30 deg) evaluates to 8.000000000000002
    return math.ceil(round(2.0 / math.sin(math.radians(pressure_angle_deg)) ** 2, 9))


def check_set(n_sun, n_planet, n_ring, n_planets, module, pressure_angle_deg,
              allow_undercut=False, min_tip_gap=1.0):
    """Return (errors, info). Errors mean the set cannot be built as specified."""
    errors, info = [], {}
    if n_planets < 2:
        return ['N_PLANETS must be at least 2'], info
    # coaxial: the planets must bridge sun and ring exactly
    if n_ring != n_sun + 2 * n_planet:
        errors.append('coaxial condition fails: N_RING must equal N_SUN + 2*N_PLANET '
                      '(%d != %d)' % (n_ring, n_sun + 2 * n_planet))
    # assembly: equally spaced planets only fit if (Ns + Nr) / n is a whole number
    if (n_sun + n_ring) % n_planets:
        errors.append('assembly condition fails: (N_SUN + N_RING) / N_PLANETS = %.3f, '
                      'must be a whole number' % ((n_sun + n_ring) / n_planets))
    # undercut: below 2/sin^2(alpha) teeth the cutter removes the involute root
    zmin = min_teeth_no_undercut(pressure_angle_deg)
    info['min_teeth_no_undercut'] = zmin
    for name, n in (('N_SUN', n_sun), ('N_PLANET', n_planet)):
        if n < zmin and not allow_undercut:
            errors.append('%s = %d undercuts at %.1f deg pressure angle (minimum %d)'
                          % (name, n, pressure_angle_deg, zmin))
    # internal gear: the ring tip circle must stay outside its base circle
    rp_ring = module * n_ring / 2.0
    if rp_ring - module <= rp_ring * math.cos(math.radians(pressure_angle_deg)):
        errors.append('ring tip circle falls inside its base circle; use more ring teeth')
    if n_ring - n_planet < 12:
        errors.append('N_RING - N_PLANET = %d; below about 12 the planet tips foul '
                      'the ring teeth (tip interference)' % (n_ring - n_planet))
    # neighbouring planets must not touch
    centre = module * (n_sun + n_planet) / 2.0
    pitch_chord = 2.0 * centre * math.sin(math.pi / n_planets)
    tip_dia = module * (n_planet + 2)
    info['planet_tip_gap'] = pitch_chord - tip_dia
    if pitch_chord - tip_dia < min_tip_gap:
        errors.append('planets too close: tip-to-tip gap %.2f mm (need %.2f)'
                      % (pitch_chord - tip_dia, min_tip_gap))
    info['centre_distance'] = centre
    info['ratio_per_stage'] = 1.0 + n_ring / n_sun
    return errors, info


def helix_rate(n, module, helix_deg, hand):
    """Tooth phase change in deg per mm of face width. hand: +1 or -1."""
    return hand * math.degrees(math.tan(math.radians(helix_deg)) / (module * n / 2.0))


# Phase conventions (all angles in degrees, global frame, at one common height z):
#   phi_s  angle of a sun tooth centreline
#   phi_r  angle of a ring tooth centreline
#   phi_p  angle of a planet tooth centreline (about the planet's own centre)
#   theta  angular position of the planet centre
# External mesh (sun-planet): a sun tooth on the line of centres meets a planet gap.
#   Ns(theta - phi_s)/360 + Np(theta + 180 - phi_p)/360 = 1/2  (mod 1)
# Internal mesh (planet-ring): a ring tooth on the line of centres meets a planet gap.
#   Np(theta - phi_p)/360 - Nr(theta - phi_r)/360 = 1/2  (mod 1)
# Eliminating phi_p gives the stage consistency condition
#   180 Np + (Ns + Nr) theta - Ns phi_s - Nr phi_r = 0  (mod 360)
# which holds for every planet at theta + 360k/n exactly when (Ns + Nr)/n is whole:
# the assembly condition falls out of the same equation.
# For 18/18/54 these reduce to P = 10 + 2T - S and P = 3R - 2T - 10 (mod 20).

def planet_phase_external(theta, phi_s, ns, np_):
    return wrap(theta + 180.0 + (ns / np_) * (theta - phi_s) - 180.0 / np_, 360.0 / np_)


def planet_phase_internal(theta, phi_r, np_, nr):
    return wrap(theta - (nr / np_) * (theta - phi_r) - 180.0 / np_, 360.0 / np_)


def sun_phase_for(theta, phi_r, ns, np_, nr):
    """Sun phase that puts a planet at theta into mesh with both sun and ring."""
    return wrap((180.0 * np_ + (ns + nr) * theta - nr * phi_r) / ns, 360.0 / ns)


def carrier_angle_for(phi_s, phi_r, ns, np_, nr):
    """Planet position that meshes with a sun whose phase is already fixed."""
    return wrap(ns * phi_s + nr * phi_r - 180.0 * np_, 360.0) / (ns + nr)


def solve_stage(ns, np_, nr, n_planets, phi_r, theta=None, phi_s=None):
    """Solve one stage. Give either the carrier angle theta (sun is free) or the
    sun phase phi_s (carrier is free). Returns (phi_s, [(theta_k, phi_p_k), ...])."""
    if phi_s is None:
        theta = 0.0 if theta is None else theta
        phi_s = sun_phase_for(theta, phi_r, ns, np_, nr)
    else:
        theta = carrier_angle_for(phi_s, phi_r, ns, np_, nr)
    planets = []
    for k in range(n_planets):
        th = theta + 360.0 * k / n_planets
        planets.append((th, planet_phase_external(th, phi_s, ns, np_)))
    return phi_s, planets


def mesh_errors(ns, np_, nr, phi_s, phi_r, planets):
    """Residual of both mesh conditions for each planet, in planet teeth (0 = exact)."""
    out = []
    for th, pp in planets:
        e = (ns * (th - phi_s) + np_ * (th + 180.0 - pp)) / 360.0 - 0.5
        i = (np_ * (th - pp) - nr * (th - phi_r)) / 360.0 - 0.5
        out.append((abs(e - round(e)), abs(i - round(i))))
    return out


def external_profile(m, n, thinning, phi, pa_deg, npts=16, ntip=5, overlap=1.0):
    """Closed outline of one external involute tooth, centred on angle phi (rad).
    Starts 'overlap' mm below the root circle so it fuses with the root disk."""
    pa = math.radians(pa_deg)
    rp = m * n / 2.0
    rb = rp * math.cos(pa)
    ra = rp + m
    rf = rp - 1.25 * m
    psi = math.pi / (2.0 * n) - (thinning / 2.0) / rp   # half tooth angle at the pitch circle
    ipa = involute(pa)

    def th(r):
        ar = 0.0 if r <= rb else math.acos(max(-1.0, min(1.0, rb / r)))
        return psi + ipa - involute(ar)

    r_lo = max(rf, rb)
    flank = [(rf - overlap, th(r_lo))]
    for i in range(npts + 1):
        r = r_lo + (ra - r_lo) * (i / float(npts)) ** 1.15
        flank.append((r, th(r)))
    tt = th(ra)
    tip = [(ra, tt - 2.0 * tt * i / float(ntip)) for i in range(1, ntip)]
    pol = flank + tip + [(r, -t) for (r, t) in reversed(flank)]
    return [(r * math.cos(t + phi), r * math.sin(t + phi)) for (r, t) in pol]


def internal_profile(m, n, thinning, phi, pa_deg, npts=14, ntip=5, overlap=1.0):
    """Closed outline of one ring tooth, centred on angle phi (rad). Runs 'overlap'
    mm beyond the root circle so it fuses with the ring annulus."""
    pa = math.radians(pa_deg)
    rp = m * n / 2.0
    rb = rp * math.cos(pa)
    ra = rp - m
    rf = rp + 1.25 * m
    psi = math.pi / (2.0 * n) - (thinning / 2.0) / rp
    ipa = involute(pa)

    def th(r):
        ar = 0.0 if r <= rb else math.acos(max(-1.0, min(1.0, rb / r)))
        return psi - ipa + involute(ar)

    r_lo = max(ra, rb)
    flank = [(r_lo + (rf - r_lo) * i / float(npts), th(r_lo + (rf - r_lo) * i / float(npts)))
             for i in range(npts + 1)]
    flank.append((rf + overlap, th(rf)))
    tt = th(r_lo)
    tip = [(r_lo, -tt + 2.0 * tt * i / float(ntip)) for i in range(1, ntip)]
    pol = [(r, -t) for (r, t) in reversed(flank)] + tip + flank
    return [(r * math.cos(t + phi), r * math.sin(t + phi)) for (r, t) in pol]


# ------------------------------------------------------------------ Fusion build
def _log(msg):
    app = adsk.core.Application.get()
    try:
        app.log(msg)
    except Exception:
        print(msg)


def _cm(v):
    return v / 10.0


def _plane_at(comp, z, name):
    pi = comp.constructionPlanes.createInput()
    pi.setByOffset(comp.xYConstructionPlane, adsk.core.ValueInput.createByReal(_cm(z)))
    pl = comp.constructionPlanes.add(pi)
    pl.name = name
    pl.isLightBulbOn = False
    return pl


def _circle_extrude(comp, radii, height, op, body=None):
    """Extrude the region between concentric circles (or a full disk for one radius)."""
    sk = comp.sketches.add(comp.xYConstructionPlane)
    for r in radii:
        sk.sketchCurves.sketchCircles.addByCenterRadius(adsk.core.Point3D.create(0, 0, 0), _cm(r))
    want = len(radii)  # disk: 1 loop; annulus: 2 loops
    prof = [p for p in sk.profiles if p.profileLoops.count == want][0]
    ei = comp.features.extrudeFeatures.createInput(prof, op)
    ei.setDistanceExtent(False, adsk.core.ValueInput.createByReal(_cm(height)))
    if body is not None:
        ei.participantBodies = [body]
    return comp.features.extrudeFeatures.add(ei)


def _tooth_loft(comp, body, profile_fn, height, twist_rad):
    """Helical tooth as a loft through rotated copies of one spline outline."""
    li = comp.features.loftFeatures.createInput(adsk.fusion.FeatureOperations.JoinFeatureOperation)
    li.participantBodies = [body]
    for i in range(LOFT_SECTIONS):
        f = i / float(LOFT_SECTIONS - 1)
        sk = comp.sketches.add(_plane_at(comp, height * f, 'Tooth section %d' % i))
        sk.isVisible = False
        pts = adsk.core.ObjectCollection.create()
        for x, y in profile_fn(twist_rad * f):
            pts.add(adsk.core.Point3D.create(_cm(x), _cm(y), 0))
        # one closed fitted spline gives a single profile and a smooth flank;
        # a polyline renders as visible facets
        spl = sk.sketchCurves.sketchFittedSplines.add(pts)
        spl.isClosed = True
        li.loftSections.add(max(sk.profiles, key=lambda p: p.areaProperties().area))
    return comp.features.loftFeatures.add(li)


def _volume(body):
    # spline-faced bodies need the high-accuracy integration; the default
    # .volume is off by ~3e-4, enough to hide a real tooth overlap
    acc = adsk.fusion.CalculationAccuracy.VeryHighCalculationAccuracy
    return body.getPhysicalProperties(acc).volume


def _pattern(comp, feature, n):
    coll = adsk.core.ObjectCollection.create()
    coll.add(feature)
    ci = comp.features.circularPatternFeatures.createInput(coll, comp.zConstructionAxis)
    ci.quantity = adsk.core.ValueInput.createByReal(n)
    ci.totalAngle = adsk.core.ValueInput.createByString('360 deg')
    # identical copies: the default 'adjust' mode re-solves each loft instance
    ci.patternComputeOption = adsk.fusion.PatternComputeOptions.IdenticalPatternCompute
    return comp.features.circularPatternFeatures.add(ci)


def _external_gear(root, name, n, hand, bore):
    occ = root.occurrences.addNewComponent(adsk.core.Matrix3D.create())
    comp = occ.component
    comp.name = name
    rf = MODULE * n / 2.0 - 1.25 * MODULE
    disk = _circle_extrude(comp, [rf], FACE_WIDTH, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    body = disk.bodies.item(0)
    body.name = name
    v0 = _volume(body)
    twist = hand * FACE_WIDTH * math.tan(math.radians(HELIX_ANGLE)) / (MODULE * n / 2.0)
    tooth = _tooth_loft(comp, body,
                        lambda phi: external_profile(MODULE, n, TOOTH_THINNING, phi, PRESSURE_ANGLE),
                        FACE_WIDTH, twist)
    v1 = _volume(body)
    _pattern(comp, tooth, n)
    v2 = _volume(body)
    _circle_extrude(comp, [bore / 2.0], FACE_WIDTH, adsk.fusion.FeatureOperations.CutFeatureOperation, body)
    return occ, _pattern_check(name, n, v1 - v0, v2 - v0)


def _ring_gear(root, height):
    occ = root.occurrences.addNewComponent(adsk.core.Matrix3D.create())
    comp = occ.component
    comp.name = 'Ring'
    rf = MODULE * N_RING / 2.0 + 1.25 * MODULE
    ann = _circle_extrude(comp, [RING_OD / 2.0, rf], height,
                          adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    body = ann.bodies.item(0)
    body.name = 'Ring'
    v0 = _volume(body)
    # internal teeth take the planet's hand
    twist = -height * math.tan(math.radians(HELIX_ANGLE)) / (MODULE * N_RING / 2.0)
    tooth = _tooth_loft(comp, body,
                        lambda phi: internal_profile(MODULE, N_RING, TOOTH_THINNING, phi, PRESSURE_ANGLE),
                        height, twist)
    v1 = _volume(body)
    _pattern(comp, tooth, N_RING)
    v2 = _volume(body)
    return occ, _pattern_check('Ring', N_RING, v1 - v0, v2 - v0)


def _pattern_check(name, n, one, total):
    # n teeth that do not overlap add n times one tooth's volume. Integration noise on
    # the spline faces is ~1e-4 of the total; any real overlap costs several % of a tooth.
    return {'part': name, 'teeth': n, 'tooth_mm3': round(one * 1000, 3),
            'pattern_mm3': round(total * 1000, 3), 'expected_mm3': round(n * one * 1000, 3),
            'ok': abs(total - n * one) <= 1e-4 * n * one}


def _matrix(rot_deg, x, y, z):
    m = adsk.core.Matrix3D.create()
    m.setToRotation(math.radians(rot_deg), adsk.core.Vector3D.create(0, 0, 1),
                    adsk.core.Point3D.create(0, 0, 0))
    m.translation = adsk.core.Vector3D.create(_cm(x), _cm(y), _cm(z))
    return m


def build(design):
    """Build the gear set into design.rootComponent. Returns a report dict."""
    errors, info = check_set(N_SUN, N_PLANET, N_RING, N_PLANETS, MODULE, PRESSURE_ANGLE,
                             ALLOW_UNDERCUT, MIN_PLANET_TIP_GAP)
    if errors:
        raise ValueError('Gear set rejected:\n- ' + '\n- '.join(errors))
    if STAGES not in (1, 2):
        raise ValueError('STAGES must be 1 or 2')
    if not 0 <= TOOTH_THINNING < 0.4 * MODULE:
        raise ValueError('TOOTH_THINNING must be between 0 and 0.4 x MODULE')
    root = design.rootComponent
    info['backlash_per_mesh'] = 2.0 * TOOTH_THINNING
    report = {'info': info, 'patterns': []}

    stage_z = [0.0, FACE_WIDTH + STAGE_GAP][:STAGES]
    ring_h = stage_z[-1] + FACE_WIDTH

    # geometry first, placement last: occurrence transforms are not timeline
    # features, so a later recompute would drop them back to the origin
    sun_occ, chk = _external_gear(root, 'Sun', N_SUN, +1, SUN_BORE)
    report['patterns'].append(chk)
    planet_occ, chk = _external_gear(root, 'Planet', N_PLANET, -1, PLANET_BORE)
    report['patterns'].append(chk)
    ring_occ, chk = _ring_gear(root, ring_h)
    report['patterns'].append(chk)

    w_r = helix_rate(N_RING, MODULE, HELIX_ANGLE, -1)
    a = info['centre_distance']
    stages = []
    placements = []   # (occurrence, matrix), applied together at the end
    for k, z in enumerate(stage_z):
        # the shared ring presents a different tooth phase at each stage height
        phi_r = w_r * z
        if k == 0:
            phi_s, planets = solve_stage(N_SUN, N_PLANET, N_RING, N_PLANETS, phi_r, theta=0.0)
        else:
            # in a gearbox the stage-2 sun is part of the stage-1 carrier, so its phase
            # is not free. Fix it (at 0) and solve the stage-2 carrier angle instead.
            phi_s, planets = solve_stage(N_SUN, N_PLANET, N_RING, N_PLANETS, phi_r, phi_s=0.0)
        resid = mesh_errors(N_SUN, N_PLANET, N_RING, phi_s, phi_r, planets)
        stages.append({'z': z, 'ring_phase': round(phi_r, 4), 'sun_phase': round(phi_s, 4),
                       'planets': [(round(t, 4), round(p, 4)) for t, p in planets],
                       'max_residual_teeth': max(max(r) for r in resid)})
        sun = sun_occ if k == 0 else root.occurrences.addExistingComponent(
            sun_occ.component, adsk.core.Matrix3D.create())
        placements.append((sun, _matrix(phi_s, 0, 0, z)))
        for j, (th, pp) in enumerate(planets):
            occ = planet_occ if (k == 0 and j == 0) else root.occurrences.addExistingComponent(
                planet_occ.component, adsk.core.Matrix3D.create())
            t = math.radians(th)
            placements.append((occ, _matrix(pp, a * math.cos(t), a * math.sin(t), z)))
    report['stages'] = stages

    # adding an occurrence can drop the transforms set on earlier ones, so place
    # everything in one pass after the last add, then read the poses back
    for occ, m in placements:
        occ.transform2 = m
    def off(occ, m):
        v = occ.transform2.translation
        v.subtract(m.translation)
        return v.length * 10.0
    report['placement_error_mm'] = max(off(occ, m) for occ, m in placements)
    if report['placement_error_mm'] > 1e-6:
        raise RuntimeError('occurrence placement did not hold (%.4f mm off)' % report['placement_error_mm'])

    report['clearances'] = _measure(root)
    return report


def _measure(root):
    """Minimum distance for every mating pair, measured on the placed solids."""
    app = adsk.core.Application.get()
    occs = [root.occurrences.item(i) for i in range(root.occurrences.count)]
    suns = [o for o in occs if o.component.name == 'Sun']
    planets = [o for o in occs if o.component.name == 'Planet']
    ring = [o for o in occs if o.component.name == 'Ring'][0]

    def gap(a, b):
        return round(app.measureManager.measureMinimumDistance(
            a.bRepBodies.item(0), b.bRepBodies.item(0)).value * 10.0, 4)

    def z_of(o):
        return round(o.transform2.translation.z * 10.0, 3)

    out = []
    for p in planets:
        sun = [s for s in suns if abs(z_of(s) - z_of(p)) < 0.01][0]
        out.append({'planet': p.name, 'z': z_of(p),
                    'to_sun': gap(p, sun), 'to_ring': gap(p, ring)})
    return out


def _summary(rep):
    lines = ['Ratio %.3f:1 per stage, %.3f:1 overall; backlash %.2f mm per mesh' % (
        rep['info']['ratio_per_stage'], rep['info']['ratio_per_stage'] ** STAGES,
        rep['info']['backlash_per_mesh'])]
    for c in rep['patterns']:
        lines.append('%s: %d teeth, pattern volume %s' % (c['part'], c['teeth'],
                                                         'OK' if c['ok'] else 'MISMATCH (teeth overlap)'))
    for s in rep['stages']:
        lines.append('Stage at z=%.1f: sun phase %.3f deg, phasing residual %.1e teeth'
                     % (s['z'], s['sun_phase'], s['max_residual_teeth']))
    for c in rep['clearances']:
        lines.append('%s (z=%.1f): to sun %.4f mm, to ring %.4f mm'
                     % (c['planet'], c['z'], c['to_sun'], c['to_ring']))
    lines.append('Planets in one stage should show identical gaps: that is the check '
                 'that the phasing is exact.')
    return '\n'.join(lines)


def run(context):
    app = adsk.core.Application.get()
    ui = app.userInterface
    try:
        doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
        design = adsk.fusion.Design.cast(doc.products.itemByProductType('DesignProductType'))
        design.designType = adsk.fusion.DesignTypes.ParametricDesignType
        rep = build(design)
        text = _summary(rep)
        _log(text)
        ui.messageBox(text, 'Planetary gear set')
    except Exception:
        ui.messageBox('Failed:\n{}'.format(traceback.format_exc()), 'Planetary gear set')
