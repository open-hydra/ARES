#!/usr/bin/env python3
"""Convert the SU2 tutorial mesh (Tecplot FEQUADRILATERAL, one structured patch) into an ARES
single-block mesh.tec: i along the axis (inlet -> outlet), j from the axis (y = 0) to the wall,
k = 2 nodes, either a planar layer (2Dplane) or a wedge rotated about x (2Daxi), as GRIB make3d does.

usage: su2tec2ares.py <mesh_SU2.dat> <mesh.tec> [planar|axi] [width_or_theta]
"""
import sys
import numpy as np


def read_fe(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    head = next(l for l in lines if l.lstrip().upper().startswith('ZONE'))
    N = int(head.split('N=')[1].split(',')[0]); E = int(head.split('E=')[1].split(',')[0])
    i0 = lines.index(head) + 1
    xy = np.array([[float(v) for v in l.split()] for l in lines[i0:i0 + N]])
    q = np.array([[int(v) - 1 for v in l.split()] for l in lines[i0 + N:i0 + N + E]])
    return xy, q


def structure(xy, q):
    edges = {}
    for e, c in enumerate(q):
        for a in range(4):
            key = frozenset((c[a], c[(a + 1) % 4]))
            edges.setdefault(key, []).append(e)
    nbr = {}
    for key, els in edges.items():
        if len(els) == 1:                                   # boundary edge
            a, b = tuple(key)
            nbr.setdefault(a, []).append(b); nbr.setdefault(b, []).append(a)
    count = np.bincount(q.ravel(), minlength=len(xy))
    corners = [n for n in range(len(xy)) if count[n] == 1]
    assert len(corners) == 4, corners
    # start corner: inlet on the axis (min x among y = min)
    c0 = min(corners, key=lambda n: (round(xy[n, 1], 12), xy[n, 0]))
    # walk the boundary along the axis (the neighbour with the smallest |dy|) to the next corner
    first = min(nbr[c0], key=lambda n: abs(xy[n, 1] - xy[c0, 1]))
    row = [c0, first]
    while row[-1] not in corners:
        nxt = [n for n in nbr[row[-1]] if n != row[-2]]
        row.append(nxt[0])
    rows = [row]
    used = set()
    while True:
        new = []
        for a in range(len(row) - 1):
            els = [e for e in edges[frozenset((row[a], row[a + 1]))] if e not in used]
            if not els:
                break
            e = els[0]; used.add(e)
            c = list(q[e]); ia = c.index(row[a]); ib = c.index(row[a + 1])
            # the quad is a b b' a' in cyclic order (either orientation)
            if (ib - ia) % 4 == 1:
                an, bn = c[(ia - 1) % 4], c[(ib + 1) % 4]
            else:
                an, bn = c[(ia + 1) % 4], c[(ib - 1) % 4]
            if a == 0:
                new.append(an)
            new.append(bn)
        if len(new) != len(row):
            break
        rows.append(new); row = new
    idx = np.array(rows).T                                  # (ni, nj)
    assert idx.size == len(xy) and len(np.unique(idx)) == len(xy)
    return xy[idx, 0], xy[idx, 1]


def write_tec(path, X, Y, Z):
    ni, nj, nk = X.shape
    with open(path, 'w') as fh:
        fh.write(' VARIABLES ="x" "y" "z"\n')
        fh.write(f' ZONE  T = Block1, I={ni}, J={nj}, K={nk}, DATAPACKING=BLOCK, VARLOCATION=([1-3]=NODAL)\n')
        for A in (X, Y, Z):
            np.savetxt(fh, A.ravel(order='F'), fmt='%22.15E')


def main():
    src, dst = sys.argv[1], sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) > 3 else 'planar'
    X2, Y2 = structure(*read_fe(src))
    if X2[-1, 0] < X2[0, 0]:
        X2, Y2 = X2[::-1], Y2[::-1]
    print(f'structured {X2.shape[0]} x {X2.shape[1]} nodes, x in [{X2.min():.4f}, {X2.max():.4f}], '
          f'y in [{Y2.min():.4f}, {Y2.max():.4f}]')
    if mode == 'axi':
        theta = np.radians(float(sys.argv[4]) if len(sys.argv) > 4 else 1.0)
        Y2 = np.maximum(Y2, 1e-8)                            # axis clamp, as GRIB (axis = 1e-8)
        z = np.array([-0.5, 0.5]) * theta
        X = np.repeat(X2[:, :, None], 2, axis=2)
        Y = Y2[:, :, None] * np.cos(z)[None, None, :]
        Z = Y2[:, :, None] * np.sin(z)[None, None, :]
    else:
        w = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-3
        X = np.repeat(X2[:, :, None], 2, axis=2); Y = np.repeat(Y2[:, :, None], 2, axis=2)
        Z = np.broadcast_to(np.array([-0.5, 0.5])[None, None, :] * w, X.shape).copy()
    write_tec(dst, X, Y, Z)
    print(f'wrote {dst} ({mode})')


if __name__ == '__main__':
    main()
