"""Draft check for moulded parts, from the B-rep face normals.

For every face of every visible solid in the active design, samples the surface
normal on a grid and works out its draft against the pull direction. Each face
is then sorted:
  plus   released by the half that moves along +PULL (draft >= MIN_DRAFT)
  minus  released by the half that moves along -PULL
  flat   perpendicular to the pull (parting / shut-off planes); needs no draft
  low    under MIN_DRAFT somewhere: will drag or need a side action
  split  normals on both sides of the pull: the face crosses the parting line
         and has to be split there, or it is an undercut

Fusion's own Draft Analysis has no API, which is why this exists. Optionally
colours the faces by band so the result can be captured as an image (this
changes face appearances in the design; parametric bodies only).

Sampling is a SAMPLES x SAMPLES grid per face, so a split confined to a very
thin strip of a large face can be missed.

Run it from Utilities > Add-Ins > Scripts and Add-Ins with the part open.
"""
import math
import traceback

try:
    import adsk.core
    import adsk.fusion
except ImportError:
    adsk = None

# ------------------------------------------------------------------ parameters
PULL = (0.0, 0.0, 1.0)   # mould opening direction, model axes
MIN_DRAFT = 0.5          # deg, below this a wall is flagged
FLAT_TOL = 89.5          # deg, faces this close to perpendicular count as flat
SAMPLES = 7              # grid samples per direction on each face
COLOUR_FACES = False     # True paints faces by band (parametric bodies only)


# ------------------------------------------------------------------ classification
def draft_deg(normal, pull):
    """Draft of a surface with unit normal n against pull p: asin(n . p).
    0 is a vertical wall, +90 a face looking along the pull."""
    d = sum(a * b for a, b in zip(normal, pull))
    return math.degrees(math.asin(max(-1.0, min(1.0, d))))


def classify(drafts, min_draft=MIN_DRAFT, flat_tol=FLAT_TOL):
    """Band for one face from its sampled drafts (deg)."""
    # a face modelled at exactly the minimum evaluates to 0.4999...; do not flag it
    lo, hi = min(drafts) + 1e-3, max(drafts) - 1e-3
    if lo >= flat_tol or hi <= -flat_tol:
        return 'flat'
    if lo >= min_draft:
        return 'plus'
    if hi <= -min_draft:
        return 'minus'
    if lo < -min_draft and hi > min_draft:
        return 'split'
    return 'low'


def _unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v)


# ------------------------------------------------------------------ Fusion side
def _face_drafts(face, pull):
    ev = face.evaluator
    rng = ev.parametricRange()
    u0, v0, u1, v1 = rng.minPoint.x, rng.minPoint.y, rng.maxPoint.x, rng.maxPoint.y
    out = []
    for i in range(SAMPLES):
        for j in range(SAMPLES):
            p = adsk.core.Point2D.create(u0 + (u1 - u0) * (i + 0.5) / SAMPLES,
                                         v0 + (v1 - v0) * (j + 0.5) / SAMPLES)
            if not ev.isParameterOnFace(p):
                continue  # trimmed-away part of the underlying surface
            ok, n = ev.getNormalAtParameter(p)
            if ok:
                out.append(draft_deg((n.x, n.y, n.z), pull))
    if not out:
        ok, n = ev.getNormalAtPoint(face.pointOnFace)
        if ok:
            out.append(draft_deg((n.x, n.y, n.z), pull))
    return out or [0.0]


COLOURS = {'plus': (60, 150, 190), 'minus': (235, 160, 60), 'flat': (200, 203, 208),
           'low': (245, 200, 50), 'split': (220, 40, 40)}


def _appearance(design, band):
    name = 'Draft check - ' + band
    a = design.appearances.itemByName(name)
    if a is None:
        app = adsk.core.Application.get()
        lib = app.materialLibraries.itemByName('Fusion Appearance Library')
        base = lib.appearances.itemByName('Plastic - Glossy (Grey)') if lib else None
        if base is None:
            return None
        a = design.appearances.addByCopy(base, name)
        rgb = COLOURS[band]
        for i in range(a.appearanceProperties.count):
            p = a.appearanceProperties.item(i)
            if p.id == 'opaque_albedo' and isinstance(p, adsk.core.ColorProperty):
                p.value = adsk.core.Color.create(rgb[0], rgb[1], rgb[2], 255)
    return a


def check(design):
    """Return {body name: {'counts': {...}, 'flagged': [...], 'min_wall_draft': deg}}."""
    pull = _unit(PULL)
    result = {}
    for name, body in _visible_bodies(design):
        counts, flagged, walls = {}, [], []
        for i, f in enumerate(body.faces):
            ds = _face_drafts(f, pull)
            band = classify(ds)
            counts[band] = counts.get(band, 0) + 1
            if band in ('plus', 'minus'):
                walls.append(min(abs(d) for d in ds))
            if band in ('low', 'split'):
                c = f.pointOnFace
                flagged.append({'face': i, 'band': band, 'draft_range': (round(min(ds), 2), round(max(ds), 2)),
                                'at_mm': (round(c.x * 10, 2), round(c.y * 10, 2), round(c.z * 10, 2)),
                                'area_mm2': round(f.area * 100, 2)})
            if COLOUR_FACES:
                a = _appearance(design, band)
                if a is not None:
                    f.appearance = a
        result[name] = {'counts': counts, 'flagged': flagged,
                             'min_wall_draft': round(min(walls), 3) if walls else None}
    return result


def _visible_bodies(design):
    """(name, body) for every visible solid, as occurrence proxies so positions are
    world. Named by occurrence: instances of one component share a body name."""
    root = design.rootComponent
    out = [(b.name, b) for b in root.bRepBodies if b.isSolid and b.isVisible]
    for occ in root.allOccurrences:
        if occ.isVisible:
            out += [(occ.name + '/' + b.name, b) for b in occ.bRepBodies if b.isSolid and b.isVisible]
    return out


def _summary(res):
    lines = ['Pull %s, minimum draft %.2f deg' % (PULL, MIN_DRAFT)]
    for name, r in res.items():
        c = r['counts']
        lines.append('%s: %s | smallest wall draft %s deg' % (
            name, ', '.join('%s %d' % (k, c[k]) for k in sorted(c)), r['min_wall_draft']))
        for f in r['flagged'][:12]:
            lines.append('   face %d %s, draft %s deg, at %s, %.1f mm2'
                         % (f['face'], f['band'], f['draft_range'], f['at_mm'], f['area_mm2']))
        if len(r['flagged']) > 12:
            lines.append('   ... %d more in the Text Commands window' % (len(r['flagged']) - 12))
    return '\n'.join(lines)


def run(context):
    app = adsk.core.Application.get()
    ui = app.userInterface
    try:
        design = adsk.fusion.Design.cast(app.activeProduct)
        if design is None:
            ui.messageBox('Open a design first (the Design workspace must be active).')
            return
        res = check(design)
        text = _summary(res)
        try:
            app.log(text)
            for name, r in res.items():
                for f in r['flagged']:
                    app.log('%s face %d %s %s %s' % (name, f['face'], f['band'], f['draft_range'], f['at_mm']))
        except Exception:
            print(text)
        ui.messageBox(text, 'Draft check')
    except Exception:
        ui.messageBox('Failed:\n{}'.format(traceback.format_exc()), 'Draft check')
