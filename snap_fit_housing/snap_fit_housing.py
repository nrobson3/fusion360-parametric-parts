"""Two-part injection-moulded wall sensor housing for Autodesk Fusion.

SH-100 Body and SH-200 Cover in ABS, held together by two cantilever snap arms.
Every feature is laid out for a straight-pull, two-plate tool: no side actions.
The two places that would normally need one are the snap catch faces, because
a catch always faces its own part's closed end:
  - body: the catch pocket is a pass-through core. It continues as a channel out
    through a slot in the hidden back face, so core steel forms it from behind.
  - cover: a notch in the plate edge over each arm lets cavity steel reach down
    and form the barb's catch face. The notch doubles as the service slot.

Snap strain, forces and the section rules are worked in snap_calc.py.

Run it from Utilities > Add-Ins > Scripts and Add-Ins. It opens a new design.
The main dimensions are written to the design as user parameters, but the
sketches are laid out from their values when the script runs: to change the
part, edit PARAMS here and run it again rather than editing parameters in Fusion.
"""
import math
import traceback

try:
    import adsk.core
    import adsk.fusion
except ImportError:
    adsk = None

# ------------------------------------------------------------------ parameters
# name, expression, unit, comment. Fusion rejects some short names (t, H) because
# they collide with unit symbols, hence the descriptive ones.
PARAMS = [
    ('W_out', '62 mm', 'mm', 'Body footprint at the wall (back face), square'),
    ('R_out', '10 mm', 'mm', 'Body plan corner radius at the back face'),
    ('H_body', '20 mm', 'mm', 'Overall depth, back face to front face'),
    ('wall', '2 mm', 'mm', 'Nominal wall (ABS)'),
    ('draft', '1.5 deg', 'deg', 'Body walls, inside and outside (fine texture on outside)'),
    ('R_back', '0.5 mm', 'mm', 'Back edge break: kept small so the catch slots cannot break out through it'),
    ('R_in', 'wall / 2', 'mm', 'Inside floor/wall radius (0.5t)'),
    ('R_front', '0.8 mm', 'mm', 'Cover front edge round is half this (the body rim is a 0.6 chamfer)'),
    ('clr', '0.1 mm', 'mm', 'Cover-to-body radial clearance at the cover underside'),
    ('lid_t', 'wall', 'mm', 'Cover plate thickness'),
    ('arm_w', '12 mm', 'mm', 'Snap arm width at root'),
    ('arm_t', '1.2 mm', 'mm', 'Snap arm thickness (0.6t)'),
    ('R_root', 'arm_t / 2', 'mm', 'Snap arm root radius (R/t = 0.5, Kt ~1.5)'),
    ('z_cf', '5.9 mm', 'mm', 'Height of the barb catch face above the back face'),
    ('barb_h', '0.7 mm', 'mm', 'Barb height (overhang off the arm face)'),
    ('crest', '0.3 mm', 'mm', 'Barb crest land'),
    ('lead', '30 deg', 'deg', 'Lead-in angle from the insertion axis'),
    ('groove_w', '12.6 mm', 'mm', 'Catch channel width at the back face'),
    ('groove_d', '0.75 mm', 'mm', 'Catch pocket depth into the wall at the pocket top'),
    ('groove_draft', '0.5 deg', 'deg', 'Draft on the pass-through core faces'),
    ('rib_t', '0.4 * wall', 'mm', 'Stop rib thickness at the tip (root ~0.54t after draft)'),
    ('rib_draft', '0.5 deg', 'deg', 'Rib side draft'),
    ('cradle_t', '0.45 * wall', 'mm', 'Cradle rib thickness at the top (root 0.5t)'),
    ('cradle_h', '5 mm', 'mm', 'Cradle rib height above the floor (2.5t)'),
    ('stop_gap', '0.85 mm', 'mm', 'Over-deflection stop behind the arm tip (caps strain ~1.05%)'),
    ('boss_od', '4.4 mm', 'mm', 'PCB boss OD (wall 0.63t)'),
    ('boss_id', '1.9 mm', 'mm', 'Pilot for a PT 2.2 thread-forming screw'),
    ('boss_h', '2.5 mm', 'mm', 'PCB boss height below the cover plate'),
    ('boss_draft', '0.5 deg', 'deg', 'Boss/pin draft'),
]


# ------------------------------------------------------------------ helpers
def _log(msg):
    app = adsk.core.Application.get()
    try:
        app.log(msg)
    except Exception:
        print(msg)


def _cm(v):
    return v / 10.0


def _pt(sk, x, y, z):
    return sk.modelToSketchSpace(adsk.core.Point3D.create(_cm(x), _cm(y), _cm(z)))


def _poly(sk, pts3):
    """Closed polyline through model-space points (mm), with shared end points."""
    L = sk.sketchCurves.sketchLines
    sps = [_pt(sk, *p) for p in pts3]
    lines = []
    for i in range(len(sps)):
        a = lines[-1].endSketchPoint if lines else sps[0]
        b = sps[i + 1] if i < len(sps) - 1 else lines[0].startSketchPoint
        lines.append(L.addByTwoPoints(a, b))
    return lines


def _rect(sk, x0, y0, x1, y1, z):
    return _poly(sk, [(x0, y0, z), (x1, y0, z), (x1, y1, z), (x0, y1, z)])


def _rrect(sk, hx, hy, r, z):
    """Rounded rectangle centred on the Z axis, with explicit tangent arcs."""
    L = sk.sketchCurves.sketchLines
    A = sk.sketchCurves.sketchArcs
    segs = [((-hx + r, -hy), (hx - r, -hy)), ((hx, -hy + r), (hx, hy - r)),
            ((hx - r, hy), (-hx + r, hy)), ((-hx, hy - r), (-hx, -hy + r))]
    cents = [(hx - r, -hy + r), (hx - r, hy - r), (-hx + r, hy - r), (-hx + r, -hy + r)]
    lines = [L.addByTwoPoints(_pt(sk, a[0], a[1], z), _pt(sk, b[0], b[1], z)) for a, b in segs]
    for i in range(4):
        c = cents[i]
        nxt = lines[(i + 1) % 4].startSketchPoint
        arc = A.addByCenterStartSweep(_pt(sk, c[0], c[1], z), lines[i].endSketchPoint, math.pi / 2)
        # sweep direction depends on the sketch orientation: fix it by checking the end point
        if arc.endSketchPoint.geometry.distanceTo(nxt.geometry) > 1e-4:
            arc.deleteMe()
            arc = A.addByCenterStartSweep(_pt(sk, c[0], c[1], z), lines[i].endSketchPoint, -math.pi / 2)
        e = arc.endSketchPoint if arc.endSketchPoint.geometry.distanceTo(nxt.geometry) < 1e-4 \
            else arc.startSketchPoint
        e.merge(nxt)
    return lines


def _plane_z(comp, z, name):
    pi = comp.constructionPlanes.createInput()
    pi.setByOffset(comp.xYConstructionPlane, adsk.core.ValueInput.createByReal(_cm(z)))
    pl = comp.constructionPlanes.add(pi)
    pl.name = name
    return pl


def _profs(sk):
    oc = adsk.core.ObjectCollection.create()
    for p in sk.profiles:
        oc.add(p)
    return oc


def _sketch(comp, plane, name):
    sk = comp.sketches.add(plane)
    sk.name = name
    return sk


def _extrude(comp, prof, dist, direction, taper, op, name, bodies=None):
    # taper > 0 grows the profile along the extrude direction
    VI = adsk.core.ValueInput
    ei = comp.features.extrudeFeatures.createInput(prof, op)
    ei.setOneSideExtent(adsk.fusion.DistanceExtentDefinition.create(VI.createByString(dist)),
                        direction, VI.createByString(taper))
    if bodies is not None:
        ei.participantBodies = bodies
    f = comp.features.extrudeFeatures.add(ei)
    f.name = name
    return f


def _edge_set(edges):
    ec = adsk.core.ObjectCollection.create()
    for e in edges:
        ec.add(e)
    return ec


# ------------------------------------------------------------------ build
def build(design):
    """Build both parts into design.rootComponent. Returns a report dict."""
    VI = adsk.core.ValueInput
    NEW = adsk.fusion.FeatureOperations.NewBodyFeatureOperation
    JOIN = adsk.fusion.FeatureOperations.JoinFeatureOperation
    CUT = adsk.fusion.FeatureOperations.CutFeatureOperation
    POS = adsk.fusion.ExtentDirections.PositiveExtentDirection
    NEG = adsk.fusion.ExtentDirections.NegativeExtentDirection
    root = design.rootComponent
    rep = {}

    up = design.userParameters
    for n, e, u, c in PARAMS:
        p = up.itemByName(n)
        if p is None:
            up.add(n, VI.createByString(e), u, c)
        else:
            p.expression = e
            p.comment = c

    def P(name):  # parameter value in mm or deg
        p = design.allParameters.itemByName(name)
        return p.value * (10.0 if p.unit == 'mm' else 180.0 / math.pi)

    W_out, R_out, H, t = P('W_out'), P('R_out'), P('H_body'), P('wall')
    d = math.radians(P('draft'))
    k = math.tan(d)
    clr, lid_t = P('clr'), P('lid_t')
    arm_w, arm_t, z_cf, barb_h, crest = P('arm_w'), P('arm_t'), P('z_cf'), P('barb_h'), P('crest')
    lead = math.radians(P('lead'))
    groove_w, groove_d = P('groove_w'), P('groove_d')
    kg = math.tan(math.radians(P('groove_draft')))
    z_pk = z_cf + 0.1                    # pocket top (0.1 axial play)
    z_lid = H - lid_t                    # cover underside

    def W(z):   # inner wall half-width of the body at height z (plan)
        return W_out / 2 + z * k - t / math.cos(d)

    def O(z):   # snap arm outer face: clr off the wall at the catch, drafted with the cavity
        return W(z_cf) - clr - (z - z_cf) * k

    def T(z):   # snap arm inner face
        return O(z) - arm_t

    def planar_at(bdy, z, nz):
        res = []
        for f in bdy.faces:
            if isinstance(f.geometry, adsk.core.Plane):
                bb = f.boundingBox
                if abs(bb.minPoint.z - _cm(z)) < 1e-5 and abs(bb.maxPoint.z - _cm(z)) < 1e-5:
                    ok, n = f.evaluator.getNormalAtPoint(f.pointOnFace)
                    if n.z * nz > 0.9:
                        res.append(f)
        return res

    # ================================================================ SH-100 BODY
    ob = root.occurrences.addNewComponent(adsk.core.Matrix3D.create())
    cb = ob.component
    cb.name = 'SH-100 Body'
    sk = _sketch(cb, cb.xYConstructionPlane, 'B01 Footprint at wall')
    _rrect(sk, W_out / 2, W_out / 2, R_out, 0)
    env = _extrude(cb, sk.profiles.item(0), 'H_body', POS, 'draft', NEW, 'B02 Envelope (draft out to front)')
    body = env.bodies.item(0)
    body.name = 'SH-100 Body'

    # back perimeter R_back (the shell turns it into 0.5t inside)
    back = [f for f in body.faces if isinstance(f.geometry, adsk.core.Plane)
            and abs(f.boundingBox.maxPoint.z) < 1e-6][0]
    fi = cb.features.filletFeatures.createInput()
    fi.addConstantRadiusEdgeSet(_edge_set(back.edges), VI.createByString('R_back'), True)
    cb.features.filletFeatures.add(fi).name = 'B03 Back edge break'

    front = [f for f in body.faces if isinstance(f.geometry, adsk.core.Plane)
             and abs(f.boundingBox.minPoint.z - _cm(H)) < 1e-6][0]
    rm = adsk.core.ObjectCollection.create()
    rm.add(front)
    si = cb.features.shellFeatures.createInput(rm, False)
    si.insideThickness = VI.createByString('wall')
    cb.features.shellFeatures.add(si).name = 'B04 Shell (uniform wall)'

    # inside floor/wall radius 0.5t
    floor = planar_at(body, t, 1)[0]
    fi = cb.features.filletFeatures.createInput()
    fi.addConstantRadiusEdgeSet(_edge_set(floor.edges), VI.createByString('R_in'), True)
    cb.features.filletFeatures.add(fi).name = 'B04a Inside floor/wall R0.5t'

    # front rim: a chamfer, not a round. The parting line needs a real edge here;
    # a round on a 1.5 deg wall straddles the PL and back-drafts the core side.
    rim = planar_at(body, H, 1)[0]
    ec = _edge_set([lp for lp in rim.loops if lp.isOuter][0].edges)
    chi = cb.features.chamferFeatures.createInput2()
    chi.chamferEdgeSets.addEqualDistanceChamferEdgeSet(ec, VI.createByString('0.6 mm'), True)
    cb.features.chamferFeatures.add(chi).name = 'B04b Front rim chamfer 0.6 (PL edge)'

    # catch channels: pass-through core from the back face up to the pocket top
    sk = _sketch(cb, cb.xYConstructionPlane, 'B05 Catch channel (pass-through core)')
    G0 = W(z_pk) + groove_d + z_pk * kg
    for s in (1, -1):
        y0, y1 = s * (G0 - 1.76), s * G0
        _rect(sk, -groove_w / 2, min(y0, y1), groove_w / 2, max(y0, y1), 0)
    _extrude(cb, _profs(sk), '%.4f mm' % z_pk, POS, '-groove_draft', CUT,
             'B06 Catch pockets + back slots', [body])
    rep['wall_behind_channel_mm'] = {'back': W_out / 2 - G0,
                                     'pocket_top': W_out / 2 + z_pk * k - (G0 - z_pk * kg)}

    # cover stop ribs on the +/-Y walls: the cover plate seats on their tops
    pl = _plane_z(cb, z_lid, 'B07 Cover seat plane')
    sk = _sketch(cb, pl, 'B08 Stop ribs')
    rt = P('rib_t')
    for sx in (1, -1):
        for sy in (1, -1):
            yin, yout = W(z_lid) - 2.0, W(z_lid) + 0.6
            _rect(sk, sx * 18 - rt / 2, sy * yin, sx * 18 + rt / 2, sy * yout, z_lid)
    _extrude(cb, _profs(sk), '%.4f mm' % (z_lid - t), NEG, 'rib_draft', JOIN, 'B09 Stop ribs (cover seat)', [body])

    # battery cradle ribs (2x AAA along X)
    ct, ch = P('cradle_t'), P('cradle_h')
    pl = _plane_z(cb, t + ch, 'B10 Cradle top plane')
    sk = _sketch(cb, pl, 'B11 Cradle ribs')
    for x in (-12, 0, 12):
        _rect(sk, x - ct / 2, -13, x + ct / 2, 13, t + ch)
    _extrude(cb, _profs(sk), 'cradle_h', NEG, 'rib_draft', JOIN, 'B12 Cradle ribs', [body])
    sk = _sketch(cb, cb.yZConstructionPlane, 'B13 AAA saddles')
    for y in (6, -6):
        sk.sketchCurves.sketchCircles.addByCenterRadius(_pt(sk, 0, y, 8.05), _cm(5.4))
    ei = cb.features.extrudeFeatures.createInput(_profs(sk), CUT)
    ei.setSymmetricExtent(VI.createByString('13 mm'), False)
    ei.participantBodies = [body]
    cb.features.extrudeFeatures.add(ei).name = 'B14 AAA saddles'

    # over-deflection stops: floor ribs just inboard of each arm tip, so a service
    # blade levered hard cannot take the arm past ~1.05 % strain
    z_arm_tip = z_cf - crest - barb_h / math.tan(lead) - 0.2
    y_stop = T(z_arm_tip) - P('stop_gap')
    pl = _plane_z(cb, z_cf + 0.3, 'B15 Stop top plane')
    sk = _sketch(cb, pl, 'B16 Over-deflection stops')
    for s in (1, -1):
        y0, y1 = s * (y_stop - 0.5 * t), s * y_stop
        _rect(sk, -5, min(y0, y1), 5, max(y0, y1), z_cf + 0.3)
    _extrude(cb, _profs(sk), '%.4f mm' % (z_cf + 0.3 - t), NEG, 'rib_draft', JOIN,
             'B17 Over-deflection stops', [body])

    # wall-mount screw holes, on a diagonal so their weld lines miss the catch strips
    sk = _sketch(cb, cb.xYConstructionPlane, 'B18 Mount holes')
    for x, y in ((-14, 18), (14, -18)):
        sk.sketchCurves.sketchCircles.addByCenterRadius(_pt(sk, x, y, 0), _cm(2.1))
    _extrude(cb, _profs(sk), '2.5 mm', POS, 'boss_draft', CUT, 'B19 Mount holes 2x d4.2', [body])

    # sprue gate well on the back face: the gate vestige stays below the mounting plane
    sk = _sketch(cb, cb.xYConstructionPlane, 'B20 Gate well')
    sk.sketchCurves.sketchCircles.addByCenterRadius(_pt(sk, 0, 0, 0), _cm(3.0))
    _extrude(cb, _profs(sk), '0.5 mm', POS, '-5 deg', CUT, 'B21 Sprue gate well d6 x0.5', [body])

    # ================================================================ SH-200 COVER
    oc = root.occurrences.addNewComponent(adsk.core.Matrix3D.create())
    cc = oc.component
    cc.name = 'SH-200 Cover'
    pl = _plane_z(cc, z_lid, 'C01 Cover underside plane')
    sk = _sketch(cc, pl, 'C02 Cover outline')
    hl = W(z_lid) - clr
    rl = R_out + z_lid * k - t / math.cos(d) - clr
    _rrect(sk, hl, hl, rl, z_lid)
    plate = _extrude(cc, sk.profiles.item(0), 'lid_t', POS, '-draft', NEW, 'C03 Cover plate (draft in to front)')
    lid = plate.bodies.item(0)
    lid.name = 'SH-200 Cover'
    top = [f for f in lid.faces if isinstance(f.geometry, adsk.core.Plane)
           and abs(f.boundingBox.minPoint.z - _cm(H)) < 1e-5 and abs(f.boundingBox.maxPoint.z - _cm(H)) < 1e-5][0]
    fi = cc.features.filletFeatures.createInput()
    fi.addConstantRadiusEdgeSet(_edge_set([lp for lp in top.loops if lp.isOuter][0].edges),
                                VI.createByString('0.5 * R_front'), True)
    cc.features.filletFeatures.add(fi).name = 'C03a Front edge round'
    rep['seam_at_front_mm'] = W(H) - (hl - lid_t * k)

    # snap arm (+Y): profile in YZ. Its catch face is formed by cavity steel through a notch.
    sk = _sketch(cc, cc.yZConstructionPlane, 'C04 Snap arm profile')
    zr = z_lid + 0.5
    zc = z_cf - crest
    zf = zc - barb_h / math.tan(lead)
    ztip = zf - 0.2
    pts_yz = [(T(zr), zr), (O(zr), zr), (O(z_cf), z_cf), (O(z_cf) + barb_h, z_cf),
              (O(zc) + barb_h, zc), (O(zf), zf), (O(ztip), ztip), (T(ztip), ztip)]
    _poly(sk, [(0, y, z) for y, z in pts_yz])
    ei = cc.features.extrudeFeatures.createInput(sk.profiles.item(0), JOIN)
    ei.setSymmetricExtent(VI.createByString('arm_w / 2'), False)
    ei.participantBodies = [lid]
    arm = cc.features.extrudeFeatures.add(ei)
    arm.name = 'C05 Snap arm +Y'
    rep['arm'] = {'length_root_to_catch': z_lid - z_cf, 'tip_z': ztip}

    mi = cc.features.mirrorFeatures.createInput(
        adsk.core.ObjectCollection.createWithArray([arm]), cc.xZConstructionPlane)
    cc.features.mirrorFeatures.add(mi).name = 'C06 Snap arm -Y (mirror)'

    # arm side draft: core-formed faces (core pulls -Z) need n_z <= 0, so the arm is
    # full width at the root and narrows 0.5 deg/side toward the tip
    sk = _sketch(cc, cc.xZConstructionPlane, 'C07 Arm side draft')
    kb = math.tan(math.radians(0.5))
    for s in (1, -1):
        _poly(sk, [(s * arm_w / 2, 0, z_lid), (s * (arm_w / 2 + 2.0), 0, z_lid),
                   (s * (arm_w / 2 + 2.0), 0, ztip - 0.5),
                   (s * (arm_w / 2 - (z_lid - ztip + 0.5) * kb), 0, ztip - 0.5)])
    ei = cc.features.extrudeFeatures.createInput(_profs(sk), CUT)
    ei.setSymmetricExtent(VI.createByString('40 mm'), False)
    ei.participantBodies = [lid]
    cc.features.extrudeFeatures.add(ei).name = 'C08 Arm sides 0.5 deg'

    # arm root radius on the inner side, which is in tension during release
    ec = adsk.core.ObjectCollection.create()
    for e in lid.edges:
        a, b = e.startVertex.geometry, e.endVertex.geometry
        if (abs(a.z - _cm(z_lid)) < 1e-5 and abs(b.z - _cm(z_lid)) < 1e-5
                and abs(abs(a.y) - _cm(T(z_lid))) < 2e-4 and abs(abs(b.y) - _cm(T(z_lid))) < 2e-4):
            ec.add(e)
    rep['arm_root_edges'] = ec.count
    fi = cc.features.filletFeatures.createInput()
    fi.addConstantRadiusEdgeSet(ec, VI.createByString('R_root'), False)
    cc.features.filletFeatures.add(fi).name = 'C08a Arm root R0.5t'

    # shut-off notches over the arms: cavity steel passes through to form the catch
    # face. A 2 mm blade pushed in here also releases the arm for service.
    sk = _sketch(cc, pl, 'C09 Shut-off / service notches')
    nw = groove_w
    for s in (1, -1):
        y0, y1 = s * O(z_lid), s * (W_out / 2 + 2)
        _rect(sk, -nw / 2, min(y0, y1), nw / 2, max(y0, y1), z_lid)
    _extrude(cc, _profs(sk), '3 mm', POS, 'draft', CUT, 'C10 Shut-off notches (service slot)', [lid])

    # PCB bosses on the underside, outside the AAA footprint so screw heads clear the cells
    sk = _sketch(cc, pl, 'C11 PCB bosses')
    for x in (19, -19):
        for y in (13, -13):
            sk.sketchCurves.sketchCircles.addByCenterRadius(_pt(sk, x, y, z_lid), _cm(P('boss_od') / 2))
    _extrude(cc, _profs(sk), 'boss_h', NEG, '-boss_draft', JOIN, 'C12 PCB bosses 4x', [lid])
    pb = _plane_z(cc, z_lid - P('boss_h'), 'C13 Boss tip plane')
    sk = _sketch(cc, pb, 'C14 Pilot holes')
    for x in (19, -19):
        for y in (13, -13):
            sk.sketchCurves.sketchCircles.addByCenterRadius(_pt(sk, x, y, z_lid - P('boss_h')),
                                                            _cm(P('boss_id') / 2))
    _extrude(cc, _profs(sk), 'boss_h', POS, '-boss_draft', CUT, 'C15 Pilots d1.9 (PT 2.2)', [lid])

    # vents: two banks (inlet low, outlet high); slots run along X, the flow direction
    sk = _sketch(cc, pl, 'C16 Vent slots')
    for yc in (17.5, 20.5, 23.5, -17.5, -20.5, -23.5):
        _rect(sk, -7, yc - 0.6, 7, yc + 0.6, z_lid)
    _extrude(cc, _profs(sk), '3 mm', POS, '1 deg', CUT, 'C17 Vent slots 6x 1.2 x 14', [lid])

    # ---------------------------------------------------------------- checks
    tm = adsk.fusion.TemporaryBRepManager.get()
    b, lb = ob.bRepBodies.item(0), oc.bRepBodies.item(0)
    rep['body'] = {'volume_mm3': round(b.volume * 1000, 3), 'lumps': b.lumps.count, 'faces': b.faces.count}
    rep['cover'] = {'volume_mm3': round(lb.volume * 1000, 3), 'lumps': lb.lumps.count, 'faces': lb.faces.count}
    A, B = tm.copy(b), tm.copy(lb)
    try:
        ok = tm.booleanOperation(A, B, adsk.fusion.BooleanTypes.IntersectionBooleanType)
        rep['interference_mm3'] = round(A.volume * 1000, 4) if ok else 'unknown (boolean failed)'
    except Exception as e:
        # two solids that do not touch can make the intersection throw; say so
        # rather than reporting a clean zero that was never measured
        rep['interference_mm3'] = 'unknown (%s)' % e
    return rep


def _summary(rep):
    return '\n'.join([
        'SH-100 Body: %.1f mm3, %d lump(s), %d faces' % (rep['body']['volume_mm3'], rep['body']['lumps'], rep['body']['faces']),
        'SH-200 Cover: %.1f mm3, %d lump(s), %d faces' % (rep['cover']['volume_mm3'], rep['cover']['lumps'], rep['cover']['faces']),
        'Body/cover interference: %s mm3 (must be 0)' % rep['interference_mm3'],
        'Wall behind catch channel: %.2f mm at the back face, %.2f mm at the pocket top'
        % (rep['wall_behind_channel_mm']['back'], rep['wall_behind_channel_mm']['pocket_top']),
        'Snap arm root edges filleted: %d (expect 2)' % rep['arm_root_edges'],
        'Run verify/draft_check to map draft per face for the +Z pull.',
    ])


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
        ui.messageBox(text, 'Snap-fit housing')
    except Exception:
        ui.messageBox('Failed:\n{}'.format(traceback.format_exc()), 'Snap-fit housing')
