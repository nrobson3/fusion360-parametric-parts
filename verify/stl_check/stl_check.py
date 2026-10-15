"""Watertight check for STL meshes. Works inside Fusion and from the command line.

A slicer will often print a broken mesh without complaint, and it fails later
as a missing wall or a part that will not seal. This checks the things that
matter for a printable solid:
  - every edge is shared by exactly two triangles (no holes, no fins)
  - each pair of neighbours walks their shared edge in opposite directions
    (consistent winding, so inside and outside are defined)
  - no zero-area triangles
  - signed volume is positive (normals point out) and, in Fusion, matches the B-rep
  - shell count, which catches loose pieces a union silently left behind

In Fusion (Utilities > Add-Ins > Scripts and Add-Ins) it meshes every visible
solid body and checks each one. From a shell:  python stl_check.py part.stl
"""
import math
import os
import struct
import sys
import tempfile
import traceback

try:
    import adsk.core
    import adsk.fusion
except ImportError:
    adsk = None

# ------------------------------------------------------------------ parameters
WELD_DECIMALS = 6        # mm; vertices that agree to this many places are one vertex
VOLUME_TOL = 0.005       # relative mesh-vs-B-rep volume difference allowed (faceting)


# ------------------------------------------------------------------ STL reading
def read_stl(path):
    """Return a list of triangles ((x,y,z), (x,y,z), (x,y,z)). ASCII or binary."""
    with open(path, 'rb') as f:
        data = f.read()
    if len(data) >= 84:
        n = struct.unpack('<I', data[80:84])[0]
        # some writers pad the end of a binary file; an ASCII file starts with 'solid'
        if 84 + 50 * n == len(data) or (84 + 50 * n < len(data) and not data.lstrip()[:5].lower() == b'solid'):
            tris = []
            for i in range(n):
                v = struct.unpack('<12f', data[84 + 50 * i: 84 + 50 * i + 48])
                tris.append((v[3:6], v[6:9], v[9:12]))
            return tris
    tris, cur = [], []
    for line in data.decode('ascii', errors='replace').splitlines():
        parts = line.split()
        if parts and parts[0].lower() == 'vertex':
            cur.append(tuple(float(c) for c in parts[1:4]))
            if len(cur) == 3:
                tris.append(tuple(cur))
                cur = []
    return tris


def write_binary_stl(path, tris):
    with open(path, 'wb') as f:
        f.write(b'\0' * 80)
        f.write(struct.pack('<I', len(tris)))
        for a, b, c in tris:
            f.write(struct.pack('<12fH', 0, 0, 0, *a, *b, *c, 0))


# ------------------------------------------------------------------ checks
def check_mesh(tris, decimals=WELD_DECIMALS):
    """Topology and volume report for a triangle soup. Vertices are welded by
    rounding, which suits exporters that write shared vertices identically."""
    index, faces = {}, []
    for t in tris:
        ids = []
        for p in t:
            key = tuple(round(c, decimals) for c in p)
            ids.append(index.setdefault(key, len(index)))
        faces.append(ids)
    verts = [None] * len(index)
    for k, i in index.items():
        verts[i] = k

    degenerate = 0
    directed = {}
    for a, b, c in faces:
        if a == b or b == c or a == c or _area(verts[a], verts[b], verts[c]) < 1e-12:
            degenerate += 1
            continue
        for u, v in ((a, b), (b, c), (c, a)):
            directed[(u, v)] = directed.get((u, v), 0) + 1

    undirected = {}
    for (u, v), n in directed.items():
        key = (u, v) if u < v else (v, u)
        undirected[key] = undirected.get(key, 0) + n
    boundary = sum(1 for n in undirected.values() if n == 1)
    nonmanifold = sum(1 for n in undirected.values() if n > 2)
    # consistent winding: an edge used twice must be used once each way
    flipped = sum(1 for (u, v), n in undirected.items()
                  if n == 2 and not (directed.get((u, v)) == 1 and directed.get((v, u)) == 1))

    volume = 0.0
    for a, b, c in faces:
        p, q, r = verts[a], verts[b], verts[c]
        volume += (p[0] * (q[1] * r[2] - q[2] * r[1]) - p[1] * (q[0] * r[2] - q[2] * r[0])
                   + p[2] * (q[0] * r[1] - q[1] * r[0])) / 6.0

    used = set(i for f in faces for i in f)
    shells = _components(faces, len(verts))
    v_count, e_count, f_count = len(used), len(undirected), len(faces) - degenerate
    return {
        'triangles': len(faces), 'vertices': v_count, 'edges': e_count,
        'boundary_edges': boundary, 'nonmanifold_edges': nonmanifold,
        'inconsistent_winding': flipped, 'degenerate': degenerate,
        'shells': shells, 'volume': volume,
        'euler': v_count - e_count + f_count,   # 2 per closed genus-0 shell
        'watertight': boundary == 0 and nonmanifold == 0 and flipped == 0 and degenerate == 0,
    }


def _area(p, q, r):
    u = [q[i] - p[i] for i in range(3)]
    v = [r[i] - p[i] for i in range(3)]
    c = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
    return 0.5 * math.sqrt(sum(x * x for x in c))


def _components(faces, nverts):
    parent = list(range(nverts))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for a, b, c in faces:
        for u, v in ((a, b), (b, c)):
            ru, rv = find(u), find(v)
            if ru != rv:
                parent[ru] = rv
    return len(set(find(i) for f in faces for i in f))


def problems(rep, expected_volume=None, vol_tol=VOLUME_TOL):
    if not rep['triangles']:
        return ['no triangles read (empty or unreadable file)']
    p = []
    for key, label in (('boundary_edges', 'open (boundary) edges'), ('nonmanifold_edges', 'non-manifold edges'),
                       ('inconsistent_winding', 'edges with inconsistent winding'),
                       ('degenerate', 'zero-area triangles')):
        if rep[key]:
            p.append('%d %s' % (rep[key], label))
    if rep['volume'] <= 0:
        p.append('signed volume %.3f: normals point inward' % rep['volume'])
    if expected_volume and abs(rep['volume'] - expected_volume) > vol_tol * expected_volume:
        p.append('mesh volume %.3f vs solid %.3f mm3' % (rep['volume'], expected_volume))
    return p


# ------------------------------------------------------------------ Fusion side
def _body_triangles(body):
    """Mesh a body with Fusion's own mesher, in mm. Skips the STL export step, which
    has produced non-manifold edges on some faceted solids."""
    calc = body.meshManager.createMeshCalculator()
    calc.setQuality(adsk.fusion.TriangleMeshQualityOptions.VeryHighQualityTriangleMesh)
    mesh = calc.calculate()
    c = mesh.nodeCoordinatesAsDouble
    idx = mesh.nodeIndices
    pts = [(c[3 * i] * 10, c[3 * i + 1] * 10, c[3 * i + 2] * 10) for i in range(len(c) // 3)]
    return [(pts[idx[k]], pts[idx[k + 1]], pts[idx[k + 2]]) for k in range(0, len(idx), 3)]


def check_design(design, out_dir=None):
    root = design.rootComponent
    bodies = [(b.name, b) for b in root.bRepBodies if b.isSolid and b.isVisible]
    for occ in root.allOccurrences:
        if occ.isVisible:
            bodies += [(occ.name + '/' + b.name, b) for b in occ.bRepBodies if b.isSolid and b.isVisible]
    out = {}
    for name, b in bodies:
        tris = _body_triangles(b)
        rep = check_mesh(tris)
        acc = adsk.fusion.CalculationAccuracy.VeryHighCalculationAccuracy
        rep['solid_volume'] = b.getPhysicalProperties(acc).volume * 1000
        rep['solid_lumps'] = b.lumps.count
        rep['problems'] = problems(rep, rep['solid_volume'])
        # a body with an enclosed cavity has one lump but two shells, so this is a note
        rep['notes'] = [] if rep['shells'] == rep['solid_lumps'] else [
            '%d mesh shells vs %d B-rep lumps (expected if the body has internal voids)'
            % (rep['shells'], rep['solid_lumps'])]
        if out_dir:
            safe = ''.join(ch if ch.isalnum() or ch in '-_' else '_' for ch in name)
            write_binary_stl(os.path.join(out_dir, safe + '.stl'), tris)
        out[name] = rep
    return out


def _summary(res):
    lines = []
    for name, r in res.items():
        state = 'PASS' if not r['problems'] else 'FAIL: ' + '; '.join(r['problems'])
        lines.append('%s: %d triangles, %d shell(s), volume %.3f mm3 -> %s'
                     % (name, r['triangles'], r['shells'], r['volume'], state))
        lines += ['   note: ' + n for n in r.get('notes', [])]
    return '\n'.join(lines) or 'No visible solid bodies.'


def run(context):
    app = adsk.core.Application.get()
    ui = app.userInterface
    try:
        design = adsk.fusion.Design.cast(app.activeProduct)
        if design is None:
            ui.messageBox('Open a design first (the Design workspace must be active).')
            return
        res = check_design(design, tempfile.gettempdir())
        text = _summary(res)
        try:
            app.log(text)
        except Exception:
            print(text)
        ui.messageBox(text + '\n\nChecked meshes written to ' + tempfile.gettempdir(), 'STL check')
    except Exception:
        ui.messageBox('Failed:\n{}'.format(traceback.format_exc()), 'STL check')


if __name__ == '__main__' and adsk is None:
    if len(sys.argv) < 2:
        print('usage: python stl_check.py part.stl [more.stl ...]')
    for path in sys.argv[1:]:
        r = check_mesh(read_stl(path))
        p = problems(r)
        print('%s: %d triangles, %d shell(s), Euler %d, volume %.3f -> %s'
              % (path, r['triangles'], r['shells'], r['euler'], r['volume'], 'PASS' if not p else 'FAIL: ' + '; '.join(p)))
