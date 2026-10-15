"""STEP export with a round-trip proof, for Autodesk Fusion.

Exporting a STEP proves nothing on its own: a hidden occurrence exports as an
empty shell, an inch design can look like millimetres if you only read SI_UNIT,
and surfaces can go out instead of solids. This script exports the active
design, inspects the file text, then imports it into a new document and
compares body count, total volume and bounding box against the original.
The imported document is closed without saving.

Run it from Utilities > Add-Ins > Scripts and Add-Ins with the design open.
"""
import os
import re
import tempfile
import traceback

try:
    import adsk.core
    import adsk.fusion
except ImportError:  # the file checks below run without Fusion
    adsk = None

# ------------------------------------------------------------------ parameters
OUTPUT_DIR = ''          # '' writes to the system temp folder
VOLUME_TOL = 1e-4        # relative; B-rep round-trip noise is ~1e-6
BBOX_TOL_MM = 0.01


# ------------------------------------------------------------------ file checks
def inspect_step(text):
    """Structural facts about a STEP file's text."""
    text = text.lstrip('\ufeff')  # a byte-order mark would hide the header
    def count(entity):
        return len(re.findall(r'\b' + entity + r'\s*\(', text))
    # the length unit is the instance that carries LENGTH_UNIT(); an inch file has
    # SI millimetres too, wrapped in a CONVERSION_BASED_UNIT, so read that first.
    # (Angles use CONVERSION_BASED_UNIT('DEGREE') in mm files, hence the filter.)
    length_unit = None
    for inst in text.split(';'):
        if 'LENGTH_UNIT' not in inst:
            continue
        conv = re.search(r"CONVERSION_BASED_UNIT\s*\(\s*'([^']+)'", inst)
        si = re.search(r'SI_UNIT\s*\(\s*([^,]*?)\s*,\s*\.METRE\.\s*\)', inst)
        if conv:
            length_unit = conv.group(1).lower()
            break
        if si and length_unit is None:
            prefix = si.group(1).strip().strip('.').lower()
            length_unit = (prefix if prefix not in ('', '$') else '') + 'metre'
    return {
        'header_ok': text.lstrip().startswith('ISO-10303-21;'),
        'terminator_ok': text.rstrip().endswith('END-ISO-10303-21;'),
        'solids': count('MANIFOLD_SOLID_BREP') + count('BREP_WITH_VOIDS'),
        'open_shells': count('OPEN_SHELL'),
        'surface_models': count('SHELL_BASED_SURFACE_MODEL'),
        'length_unit': length_unit,
        'products': re.findall(r"PRODUCT\s*\(\s*'([^']*)'", text),
    }


UNIT_NAMES = {'mm': 'millimetre', 'cm': 'centimetre', 'm': 'metre', 'in': 'inch', 'ft': 'foot'}


def file_problems(info, expected_bodies, expected_unit=None):
    """expected_bodies: distinct body definitions. A STEP assembly stores an instanced
    component once and places it several times, so visible instances do not count."""
    p = []
    if expected_unit and info['length_unit'] != expected_unit:
        p.append('file length unit is %s, design unit is %s' % (info['length_unit'], expected_unit))
    if not info['header_ok'] or not info['terminator_ok']:
        p.append('file is truncated or not ISO 10303-21')
    if info['open_shells'] or info['surface_models']:
        p.append('%d open shells / %d surface models: some bodies went out as surfaces'
                 % (info['open_shells'], info['surface_models']))
    # a shelled body with a fully enclosed cavity exports as BREP_WITH_VOIDS, so the
    # count is solids-like entities, and it may exceed the body count but not trail it
    if info['solids'] < expected_bodies:
        p.append('%d solids in the file for %d distinct bodies (hidden occurrences export empty)'
                 % (info['solids'], expected_bodies))
    return p


def compare(orig, back, vol_tol=VOLUME_TOL, bbox_tol=BBOX_TOL_MM):
    """orig/back: {'bodies': n, 'volume_mm3': v, 'bbox_mm': (x0,y0,z0,x1,y1,z1)}"""
    p = []
    if orig['bodies'] != back['bodies']:
        p.append('body count %d -> %d' % (orig['bodies'], back['bodies']))
    if abs(back['volume_mm3'] - orig['volume_mm3']) > vol_tol * max(orig['volume_mm3'], 1e-9):
        p.append('volume %.4f -> %.4f mm3' % (orig['volume_mm3'], back['volume_mm3']))
    if max(abs(a - b) for a, b in zip(orig['bbox_mm'], back['bbox_mm'])) > bbox_tol:
        p.append('bounding box %s -> %s' % (orig['bbox_mm'], back['bbox_mm']))
    return p


# ------------------------------------------------------------------ Fusion side
def _solids(design):
    root = design.rootComponent
    out = [b for b in root.bRepBodies if b.isSolid and b.isVisible]
    for occ in root.allOccurrences:
        if occ.isVisible:
            out += [b for b in occ.bRepBodies if b.isSolid and b.isVisible]
    return out


def _stats(design):
    bodies = _solids(design)
    if not bodies:
        return {'bodies': 0, 'volume_mm3': 0.0, 'bbox_mm': (0,) * 6}
    lo = [min(getattr(b.boundingBox.minPoint, a) for b in bodies) * 10 for a in 'xyz']
    hi = [max(getattr(b.boundingBox.maxPoint, a) for b in bodies) * 10 for a in 'xyz']
    # high-accuracy integration: the default volume differs by ~5e-4 between a body
    # and its re-imported copy on spline faces, which would swamp the comparison
    acc = adsk.fusion.CalculationAccuracy.VeryHighCalculationAccuracy
    vol = sum(b.getPhysicalProperties(acc).volume for b in bodies) * 1000
    return {'bodies': len(bodies), 'volume_mm3': vol,
            'bbox_mm': tuple(round(v, 4) for v in lo + hi)}


def roundtrip(design, path):
    app = adsk.core.Application.get()
    em = design.exportManager
    if not em.execute(em.createSTEPExportOptions(path, design.rootComponent)):
        raise RuntimeError('STEP export failed')
    with open(path, encoding='utf-8', errors='replace') as f:
        info = inspect_step(f.read())
    orig = _stats(design)
    distinct = len(set(b.nativeObject.entityToken if b.nativeObject else b.entityToken
                       for b in _solids(design)))
    unit = UNIT_NAMES.get(design.unitsManager.defaultLengthUnits)
    rep = {'path': path, 'file': info, 'original': orig,
           'file_problems': file_problems(info, distinct, unit)}

    src_doc = app.activeDocument
    doc = app.importManager.importToNewDocument(app.importManager.createSTEPImportOptions(path))
    try:
        back_design = adsk.fusion.Design.cast(doc.products.itemByProductType('DesignProductType'))
        rep['reimported'] = _stats(back_design)
    finally:
        doc.close(False)   # never saved, so it costs no cloud document slot
        try:
            src_doc.activate()
        except Exception:
            pass
    rep['roundtrip_problems'] = compare(orig, rep['reimported'])
    return rep


def _summary(rep):
    o, b, f = rep['original'], rep['reimported'], rep['file']
    lines = ['Wrote %s' % rep['path'],
             'File: %d solids, unit %s, products %s' % (f['solids'], f['length_unit'], ', '.join(f['products'][:6])),
             'Original:   %d bodies, %.3f mm3, bbox %s' % (o['bodies'], o['volume_mm3'], o['bbox_mm']),
             'Re-imported: %d bodies, %.3f mm3, bbox %s' % (b['bodies'], b['volume_mm3'], b['bbox_mm'])]
    problems = rep['file_problems'] + rep['roundtrip_problems']
    lines.append('PASS' if not problems else 'FAIL:\n- ' + '\n- '.join(problems))
    if any(re.match(r'\d{4}-\d{2}-\d{2}-', p) for p in f['products']):
        lines.append('Note: some PRODUCT names are Fusion auto ids. STEP names come from '
                     'component part numbers; set them before sending the file.')
    return '\n'.join(lines)


def run(context):
    app = adsk.core.Application.get()
    ui = app.userInterface
    try:
        design = adsk.fusion.Design.cast(app.activeProduct)
        if design is None:
            ui.messageBox('Open a design first (the Design workspace must be active).')
            return
        out_dir = OUTPUT_DIR or tempfile.gettempdir()
        name = re.sub(r'[^\w.-]+', '_', app.activeDocument.name) or 'design'
        rep = roundtrip(design, os.path.join(out_dir, name + '.step'))
        text = _summary(rep)
        try:
            app.log(text)
        except Exception:
            print(text)
        ui.messageBox(text, 'STEP round-trip')
    except Exception:
        ui.messageBox('Failed:\n{}'.format(traceback.format_exc()), 'STEP round-trip')
