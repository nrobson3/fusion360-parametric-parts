"""Single-piece compliant parallel-jaw gripper for Autodesk Fusion.

Prints as one part with no pins, screws or springs. Each jaw hangs from a
four-bar parallelogram of two links with a flexure hinge at each end, so the
jaws stay parallel through the stroke. Squeezing the handles closes the jaws.

The flexure is sized from the material's strain limit before anything is
modelled; the script refuses a geometry that would yield. Then it builds the
part and measures it back against the numbers it was built from.

Run it from Utilities > Add-Ins > Scripts and Add-Ins. It opens a new design;
results go to the Text Commands window and a summary dialog.
"""
import math
import traceback

try:
    import adsk.core
    import adsk.fusion
except ImportError:  # lets the sizing math be imported and tested outside Fusion
    adsk = None

# ------------------------------------------------------------------ parameters
NOZZLE = 0.4             # mm, printer nozzle / extrusion width
FLEXURE_LINES = 2        # perimeters across the flexure; t = 2 x 0.4 = 0.8 mm
FLEXURE_LENGTH = 8.0     # mm, necked length of each hinge
LINK_LENGTH = 25.0       # mm, flexure centre to flexure centre
DESIGN_ANGLE = 12.0      # deg, hinge rotation at full squeeze
MATERIAL = 'PLA'         # sized against this one; the report covers all of MATERIALS
MIN_SAFETY = 1.5         # required yield strain / working strain
MATERIALS = {            # tensile strain at yield (typical datasheet values for FDM grades)
    'PLA': 0.020,
    'PETG': 0.040,
}
LINK_WIDTH = 8.0         # mm, rigid part of each link
LINK_GAP = 6.0           # mm, between the inner and outer link
JAW_GAP = 14.0           # mm, jaw opening at rest
JAW_DEPTH = 7.0          # mm, jaw from its face to the inner link
BASE_HEIGHT = 12.0       # mm
BASE_MARGIN = 6.0        # mm, base beyond the outer link
CARRIER_HEIGHT = 8.0     # mm, jaw carrier joining the two links
HANDLE_LENGTH = 32.0     # mm
THICKNESS = 8.0          # mm, extruded in Z; prints flat in this orientation
ROOT_FILLET = 0.4        # mm, at both ends of every flexure
TOOTH_PITCH = 2.0        # mm, V-grooves across each jaw face
TOOTH_DEPTH = 0.6        # mm


# ------------------------------------------------------------------ sizing math
def flexure_thickness(nozzle=NOZZLE, lines=FLEXURE_LINES):
    # a whole number of perimeters prints as solid walls with no infill; at an
    # in-between width the slicer leaves a void where the stress peaks
    return nozzle * lines


def flexure_strain(t, theta_deg, length):
    """Peak bending strain of a uniform flexure bent through theta.
    Curvature theta/L, outer fibre at t/2: eps = t * theta / (2 L)."""
    return t * math.radians(theta_deg) / (2.0 * length)


def max_thickness(theta_deg, length, eps_allow):
    """Thickest flexure that stays inside eps_allow at theta."""
    return 2.0 * length * eps_allow / math.radians(theta_deg)


def kinematics(link_len, theta_deg, jaw_gap):
    """Parallelogram jaw motion: each jaw moves link_len*sin(theta) inward and
    drops link_len*(1 - cos(theta))."""
    th = math.radians(theta_deg)
    travel = link_len * math.sin(th)
    return {'travel_per_jaw': travel, 'drop': link_len * (1.0 - math.cos(th)),
            'gap_closed': jaw_gap - 2.0 * travel}


def size(material=MATERIAL):
    """Check the hinge against every material. Raises if the design material fails."""
    t = flexure_thickness()
    eps = flexure_strain(t, DESIGN_ANGLE, FLEXURE_LENGTH)
    table = {}
    for name, eps_y in MATERIALS.items():
        table[name] = {'yield_strain': eps_y, 'safety': eps_y / eps,
                       't_max_for_min_safety': max_thickness(DESIGN_ANGLE, FLEXURE_LENGTH,
                                                             eps_y / MIN_SAFETY)}
    if DESIGN_ANGLE <= 0 or not 0 < t < LINK_WIDTH:
        raise ValueError('need DESIGN_ANGLE > 0 and 0 < flexure thickness < LINK_WIDTH')
    kin = kinematics(LINK_LENGTH, DESIGN_ANGLE, JAW_GAP)
    problems = []
    if table[material]['safety'] < MIN_SAFETY:
        problems.append('%s safety %.2f < %.2f: thin the flexure to %.2f mm or lengthen it'
                        % (material, table[material]['safety'], MIN_SAFETY,
                           table[material]['t_max_for_min_safety']))
    if kin['gap_closed'] <= 0:
        problems.append('jaws collide before %.1f deg (closed gap %.2f mm)'
                        % (DESIGN_ANGLE, kin['gap_closed']))
    if 2 * ROOT_FILLET >= FLEXURE_LENGTH:
        problems.append('root fillets use up the whole flexure length')
    return {'t': t, 'strain': eps, 'materials': table, 'kinematics': kin, 'problems': problems}


def layout():
    """Rectangles (x0, y0, x1, y1) for the +X half, plus derived stations. mm."""
    x_jaw = JAW_GAP / 2.0
    x_in0 = x_jaw + JAW_DEPTH
    x_in1 = x_in0 + LINK_WIDTH
    x_out0 = x_in1 + LINK_GAP
    x_out1 = x_out0 + LINK_WIDTH
    y_f1 = BASE_HEIGHT + FLEXURE_LENGTH / 2.0       # lower flexure centre
    y_f2 = y_f1 + LINK_LENGTH                       # upper flexure centre
    y_top = y_f2 + FLEXURE_LENGTH / 2.0             # top of the links
    y_car = y_top + CARRIER_HEIGHT
    half = x_out1 + BASE_MARGIN
    parts = {
        'inner_link': (x_in0, BASE_HEIGHT, x_in1, y_top),
        'outer_link': (x_out0, BASE_HEIGHT, x_out1, y_top),
        'carrier': (x_in0, y_top, x_out1, y_car),
        'jaw': (x_jaw, y_top, x_in0, y_car),
        'handle': (x_out0, y_car, x_out1, y_car + HANDLE_LENGTH),
    }
    return {'base': (-half, 0.0, half, BASE_HEIGHT), 'half': parts,
            'flexure_x': ((x_in0 + x_in1) / 2.0, (x_out0 + x_out1) / 2.0),
            'flexure_y': (y_f1, y_f2), 'window': (x_in1, x_out0, BASE_HEIGHT, y_top),
            'y_top': y_top, 'y_car': y_car}


def expected_volumes(lay, t):
    """Volumes in mm^3 at each build step, from plain geometry."""
    def area(r):
        return (r[2] - r[0]) * (r[3] - r[1])
    outline = area(lay['base']) + 2 * sum(area(r) for r in lay['half'].values())
    notches = 8 * (LINK_WIDTH - t) * FLEXURE_LENGTH          # 8 hinges
    fillets = 32 * ROOT_FILLET ** 2 * (1 - math.pi / 4)      # 4 concave corners per hinge
    n_teeth = int(round(CARRIER_HEIGHT / TOOTH_PITCH, 9))
    teeth = 2 * n_teeth * TOOTH_DEPTH * TOOTH_PITCH / 2.0
    v0 = outline * THICKNESS
    v1 = v0 - notches * THICKNESS
    v2 = v1 + fillets * THICKNESS
    v3 = v2 - teeth * THICKNESS
    return {'outline': v0, 'notched': v1, 'filleted': v2, 'final': v3, 'n_teeth': n_teeth}


# ------------------------------------------------------------------ Fusion build
def _log(msg):
    app = adsk.core.Application.get()
    try:
        app.log(msg)
    except Exception:
        print(msg)


def _cm(v):
    return v / 10.0


def _rect(lines, x0, y0, x1, y1):
    lines.addTwoPointRectangle(adsk.core.Point3D.create(_cm(x0), _cm(y0), 0),
                               adsk.core.Point3D.create(_cm(x1), _cm(y1), 0))


def _cut(comp, body, sketch):
    coll = adsk.core.ObjectCollection.create()
    for p in sketch.profiles:
        coll.add(p)
    ext = comp.features.extrudeFeatures
    ci = ext.createInput(coll, adsk.fusion.FeatureOperations.CutFeatureOperation)
    ci.participantBodies = [body]
    ci.setDistanceExtent(False, adsk.core.ValueInput.createByReal(_cm(THICKNESS)))
    return ext.add(ci)


def build(design):
    """Build the gripper into design.rootComponent. Returns a report dict."""
    sizing = size()
    if sizing['problems']:
        raise ValueError('Sizing rejected:\n- ' + '\n- '.join(sizing['problems']))
    t = sizing['t']
    lay = layout()
    want = expected_volumes(lay, t)
    root = design.rootComponent
    rep = {'sizing': sizing, 'volumes': {}}

    def vol(body):
        return body.volume * 1000.0

    # 1. outline: base plus mirrored links, carriers, jaws and handles
    sk = root.sketches.add(root.xYConstructionPlane)
    sk.name = 'Gripper outline'
    sk.isComputeDeferred = True
    lines = sk.sketchCurves.sketchLines
    _rect(lines, *lay['base'])
    for (x0, y0, x1, y1) in lay['half'].values():
        _rect(lines, x0, y0, x1, y1)
        _rect(lines, -x1, y0, -x0, y1)
    sk.isComputeDeferred = False

    # the region between two links is enclosed by base, links and carrier, so it
    # is a profile too. Extruding it would weld the links solid: skip it by centroid.
    wx0, wx1, wy0, wy1 = lay['window']
    profs = adsk.core.ObjectCollection.create()
    skipped = 0
    for p in sk.profiles:
        c = p.areaProperties().centroid
        if wx0 < abs(c.x * 10) < wx1 and wy0 < c.y * 10 < wy1:
            skipped += 1
            continue
        profs.add(p)
    ext = root.features.extrudeFeatures
    ei = ext.createInput(profs, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    ei.setDistanceExtent(False, adsk.core.ValueInput.createByReal(_cm(THICKNESS)))
    body = ext.add(ei).bodies.item(0)
    body.name = 'Gripper'
    rep['windows_skipped'] = skipped
    rep['volumes']['outline'] = vol(body)

    # 2. flexures: neck each link to t over FLEXURE_LENGTH at both ends
    sk = root.sketches.add(root.xYConstructionPlane)
    sk.name = 'Flexure notches'
    sk.isComputeDeferred = True
    lines = sk.sketchCurves.sketchLines
    for cx in lay['flexure_x']:
        for cy in lay['flexure_y']:
            ya, yb = cy - FLEXURE_LENGTH / 2.0, cy + FLEXURE_LENGTH / 2.0
            for s in (1.0, -1.0):
                c = cx * s
                _rect(lines, c - LINK_WIDTH / 2.0, ya, c - t / 2.0, yb)
                _rect(lines, c + t / 2.0, ya, c + LINK_WIDTH / 2.0, yb)
    sk.isComputeDeferred = False
    _cut(root, body, sk)
    rep['volumes']['notched'] = vol(body)

    # 3. root fillets: eps = t*theta/2L assumes a uniform beam; real flexures fail
    #    at the sharp internal corners, and a fillet there roughly halves the peak
    corners = []
    for cx in lay['flexure_x']:
        for cy in lay['flexure_y']:
            for s in (1.0, -1.0):
                for x in (cx * s - t / 2.0, cx * s + t / 2.0):
                    for y in (cy - FLEXURE_LENGTH / 2.0, cy + FLEXURE_LENGTH / 2.0):
                        corners.append((x, y))
    edges = adsk.core.ObjectCollection.create()
    for e in body.edges:
        a, b = e.startVertex.geometry, e.endVertex.geometry
        if abs(a.x - b.x) > 1e-6 or abs(a.y - b.y) > 1e-6:
            continue  # only the vertical edges, which run the full thickness
        if any(abs(a.x * 10 - x) < 0.01 and abs(a.y * 10 - y) < 0.01 for x, y in corners):
            edges.add(e)
    rep['fillet_edges'] = edges.count
    fi = root.features.filletFeatures.createInput()
    fi.edgeSetInputs.addConstantRadiusEdgeSet(
        edges, adsk.core.ValueInput.createByReal(_cm(ROOT_FILLET)), False)
    root.features.filletFeatures.add(fi)
    rep['volumes']['filleted'] = vol(body)

    # 4. V-grooves across both jaw faces
    sk = root.sketches.add(root.xYConstructionPlane)
    sk.name = 'Jaw teeth'
    sk.isComputeDeferred = True
    lines = sk.sketchCurves.sketchLines
    x_face = JAW_GAP / 2.0
    for k in range(want['n_teeth']):
        yc = lay['y_top'] + TOOTH_PITCH * (k + 0.5)
        for s in (1.0, -1.0):
            p1 = adsk.core.Point3D.create(_cm(x_face * s), _cm(yc - TOOTH_PITCH / 2.0), 0)
            p2 = adsk.core.Point3D.create(_cm((x_face + TOOTH_DEPTH) * s), _cm(yc), 0)
            p3 = adsk.core.Point3D.create(_cm(x_face * s), _cm(yc + TOOTH_PITCH / 2.0), 0)
            lines.addByTwoPoints(p1, p2)
            lines.addByTwoPoints(p2, p3)
            lines.addByTwoPoints(p3, p1)
    sk.isComputeDeferred = False
    _cut(root, body, sk)
    rep['volumes']['final'] = vol(body)

    rep['expected'] = want
    rep['bodies'] = root.bRepBodies.count
    rep['lumps'] = body.lumps.count
    rep['measured_flexure_t'] = _measure_flexure(body, lay)
    return rep


def _measure_flexure(body, lay):
    """Flexure thickness from the planar faces, not from vertices: fillet tangent
    points sit inside the hinge band and make a vertex scan read 2t."""
    cx = lay['flexure_x'][0]
    cy = lay['flexure_y'][0]
    xs = set()
    for f in body.faces:
        if f.geometry.surfaceType != adsk.core.SurfaceTypes.PlaneSurfaceType:
            continue
        n = adsk.core.Plane.cast(f.geometry).normal
        if abs(abs(n.x) - 1.0) > 1e-6:
            continue
        bb = f.boundingBox
        x = bb.minPoint.x * 10
        if abs(x - cx) < LINK_WIDTH / 2.0 - 0.01 and abs(bb.minPoint.y * 10 - cy) < FLEXURE_LENGTH / 2.0 \
                and abs(bb.maxPoint.y * 10 - cy) < FLEXURE_LENGTH / 2.0:
            xs.add(round(x, 4))
    return round(max(xs) - min(xs), 4) if len(xs) >= 2 else None


def _summary(rep):
    s = rep['sizing']
    lines = ['Flexure t = %.2f mm (%d x %.1f mm lines), strain at %.0f deg = %.3f %%'
             % (s['t'], FLEXURE_LINES, NOZZLE, DESIGN_ANGLE, s['strain'] * 100)]
    for name, m in s['materials'].items():
        lines.append('  %s: yield %.1f %%, safety %.2fx' % (name, m['yield_strain'] * 100, m['safety']))
    k = s['kinematics']
    lines.append('Jaw travel %.3f mm each, gap %.2f -> %.3f mm, carrier drop %.3f mm'
                 % (k['travel_per_jaw'], JAW_GAP, k['gap_closed'], k['drop']))
    for step in ('outline', 'notched', 'filleted', 'final'):
        got, want = rep['volumes'][step], rep['expected'][step]
        lines.append('Volume %-8s %10.3f mm3  expected %10.3f  %s'
                     % (step, got, want, 'OK' if abs(got - want) <= 1e-5 * want else 'CHECK'))
    lines.append('Bodies %d, lumps %d (must both be 1), windows skipped %d (expect 2), '
                 'fillet edges %d (expect 32)' % (rep['bodies'], rep['lumps'],
                                                  rep['windows_skipped'], rep['fillet_edges']))
    lines.append('Measured flexure thickness %s mm' % rep['measured_flexure_t'])
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
        ui.messageBox(text, 'Compliant gripper')
    except Exception:
        ui.messageBox('Failed:\n{}'.format(traceback.format_exc()), 'Compliant gripper')
