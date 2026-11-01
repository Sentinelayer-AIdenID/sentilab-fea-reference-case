"""Conforming P1 tetrahedral conduction on a rectilinear Kuhn mesh (SI units).

Each box is split along 000--111 into six tetrahedra. Exact element
integration reduces the diffusion matrix to axis-edge contributions; this
avoids storing millions of tetrahedron connectivity records.
"""

import time
from itertools import pairwise, permutations

import numpy as np
import pyamg
from scipy import sparse
from scipy.sparse.linalg import cg

PERMS = tuple(permutations(range(3)))


def unique_edges(values):
    """Merge only floating point copies of the same geometric coordinate."""
    out = []
    for v in sorted(values):
        if not out or v - out[-1] > 2e-16:
            out.append(v)
    return np.array(out)


def axis_mesh(edges, lo, hi, level, inner_h=0.003, outer_h=0.012):
    base = unique_edges(edges)
    out = []
    for a, b in pairwise(base):
        h = inner_h if lo <= (a + b) / 2 <= hi else outer_h
        n = max(1, int(np.ceil((b - a) / h))) * 2**level
        out.extend(np.linspace(a, b, n + 1)[:-1])
    return np.r_[out, base[-1]]


def mesh_axes(data, case_id, level):
    case = data["cases"][case_id]
    axes = []
    for ax, width in [("x", "width"), ("y", "height")]:
        edges = []
        for layer in data["stack"]:
            edges += [layer["extent_m"][ax + "0"], layer["extent_m"][ax + "1"]]
        for blocks in case["layers"].values():
            for b in blocks:
                edges += [b[ax + "0_m"], b[ax + "0_m"] + b[width + "_m"]]
        # Put centre readouts on nodes at every level.
        edges += [d[ax + "_centre_m"] for d in case["chiplets"]]
        axes.append(axis_mesh(edges, 0.0, 0.045, level))
    z, layer_ids = [], []
    counts = [1, 1, 1, 1, 2, 1, 2, 6]
    for layer, n in zip(data["stack"], counts):
        n *= 2**level
        z.extend(np.linspace(layer["z_bottom_m"], layer["z_top_m"], n + 1)[:-1])
        layer_ids.extend([layer["index"]] * n)
    z.append(data["stack"][-1]["z_top_m"])
    return (*axes, np.array(z)), np.array(layer_ids), np.array(counts) * 2**level


def rectangle_mask(x, y, b):
    return (
        (x[:, None] > b["x0_m"])
        & (x[:, None] < b["x0_m"] + b["width_m"])
        & (y[None, :] > b["y0_m"])
        & (y[None, :] < b["y0_m"] + b["height_m"])
    )


def layer_fields(data, case_id, axes):
    x, y, _ = axes
    xc, yc = (x[:-1] + x[1:]) / 2, (y[:-1] + y[1:]) / 2
    case = data["cases"][case_id]
    fields = []
    for layer in data["stack"]:
        e = layer["extent_m"]
        inside = (
            (xc[:, None] > e["x0"])
            & (xc[:, None] < e["x1"])
            & (yc[None, :] > e["y0"])
            & (yc[None, :] < e["y1"])
        )
        k = np.zeros(inside.shape)
        covered = np.zeros(inside.shape, dtype=np.int8)
        if layer["index"] < 6:
            for b in case["layers"]["L" + str(layer["index"])]:
                mask = rectangle_mask(xc, yc, b)
                k[mask] = b["k_W_per_mK"]
                covered[mask] += 1
            for gap in (
                case["tiling_check"]
                .get("L" + str(layer["index"]), {})
                .get("uncovered_rectangles_m", [])
            ):
                a, b, c, d = gap
                mask = (
                    (xc[:, None] > a)
                    & (xc[:, None] < b)
                    & (yc[None, :] > c)
                    & (yc[None, :] < d)
                )
                k[mask] = data["materials"]["underfill"]["k_W_per_mK"]
                covered[mask] += 1
            if np.any(covered[inside] != 1) or np.any(covered[~inside]):
                raise ValueError("Material rectangles do not tile the mesh exactly")
        else:
            k[inside] = layer["default_k_W_per_mK"]
        fields.append(k)
    q = np.zeros(fields[0].shape)
    flux = q.copy()
    for die in case["chiplets"]:
        mask = rectangle_mask(xc, yc, die)
        q[mask] = die["volumetric_source_W_m3"]
        flux[mask] = die["power_density_W_m2"]
    return fields, q, flux


def assemble(axes, layer_ids, fields, q, flux, source_layer, h):
    """Exact P1 stiffness, two load vectors, and integrated Robin weights.

    Array storage is (z, x, y). The returned unknown is theta=T-T_ref.
    An edge parallel to axis a receives k*V/(3*da^2) at the two
    equal-bit transverse corners, and k*V/(6*da^2) at the other two.
    """
    x, y, z = axes
    dx, dy, dz = np.diff(x), np.diff(y), np.diff(z)
    shape = (len(z), len(x), len(y))
    gx = np.zeros((len(z), len(dx), len(y)))
    gy = np.zeros((len(z), len(x), len(dy)))
    gz = np.zeros((len(dz), len(x), len(y)))
    load = np.zeros(shape)
    sheet = np.zeros(shape)
    active_boxes = 0
    for iz, lid in enumerate(layer_ids):
        k = fields[lid]
        active_boxes += np.count_nonzero(k)
        vx = k * dz[iz] * dy[None, :] / dx[:, None]
        vy = k * dz[iz] * dx[:, None] / dy[None, :]
        vz = k * dx[:, None] * dy[None, :] / dz[iz]
        for a in (0, 1):
            for b in (0, 1):
                w = 1 / (3 if a == b else 6)
                gx[iz + a, :, b : b + len(dy)] += w * vx
                gy[iz + a, b : b + len(dx), :] += w * vy
                gz[iz, a : a + len(dx), b : b + len(dy)] += w * vz
        if lid == source_layer:
            vq = q * dx[:, None] * dy[None, :] * dz[iz]
            for a in (0, 1):
                for b in (0, 1):
                    for c in (0, 1):
                        w = 1 / (4 if a == b == c else 12)
                        load[iz + c, a : a + len(dx), b : b + len(dy)] += w * vq
            if iz == 0 or layer_ids[iz - 1] != source_layer:
                sf = flux * dx[:, None] * dy[None, :]
                for a in (0, 1):
                    for b in (0, 1):
                        sheet[iz, a : a + len(dx), b : b + len(dy)] += sf / (
                            3 if a == b else 6
                        )
    diag = np.zeros(shape)
    for g, axis in [(gx, 1), (gy, 2), (gz, 0)]:
        low, high = [slice(None)] * 3, [slice(None)] * 3
        low[axis], high[axis] = slice(None, -1), slice(1, None)
        diag[tuple(low)] += g
        diag[tuple(high)] += g
    active = diag > 0
    ids = np.full(shape, -1, dtype=np.int32)
    ids[active] = np.arange(np.count_nonzero(active), dtype=np.int32)
    n = np.count_nonzero(active)
    K = sparse.diags(diag[active], format="csr")
    for g, axis in [(gx, 1), (gy, 2), (gz, 0)]:
        low, high = [slice(None)] * 3, [slice(None)] * 3
        low[axis], high[axis] = slice(None, -1), slice(1, None)
        nz = g > 0
        rows, cols = ids[tuple(low)][nz], ids[tuple(high)][nz]
        edge = sparse.coo_matrix((-g[nz], (rows, cols)), shape=(n, n)).tocsr()
        K += edge + edge.T
    # Consistent triangle mass on the top, diagonal 00--11.
    corners = [ids[-1, :-1, :-1], ids[-1, 1:, :-1], ids[-1, :-1, 1:], ids[-1, 1:, 1:]]
    area = dx[:, None] * dy[None, :]
    rows, cols, vals = [], [], []
    robin = np.zeros(n)
    for tri in [(0, 1, 3), (0, 2, 3)]:
        for a in tri:
            np.add.at(robin, corners[a].ravel(), (h * area / 6).ravel())
            for b in tri:
                rows.append(corners[a].ravel())
                cols.append(corners[b].ravel())
                vals.append((h * area * (2 if a == b else 1) / 24).ravel())
    K += sparse.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
    ).tocsr()
    return K, load[active], sheet[active], robin, ids, int(active_boxes * 6)


def solve_system(K, rhs, preconditioner=None):
    start = time.perf_counter()
    if preconditioner is None:
        hierarchy = pyamg.ruge_stuben_solver(K, max_coarse=50)
        preconditioner = hierarchy.aspreconditioner(cycle="V")
    iters = [0]

    def count(_):
        iters[0] += 1

    u, info = cg(
        K, rhs, M=preconditioner, rtol=1e-10, atol=0.0, maxiter=4000, callback=count
    )
    residual = float(np.linalg.norm(K @ u - rhs) / np.linalg.norm(rhs))
    if info or residual > 2e-8:
        raise RuntimeError(f"Linear solve failed: info={info}, residual={residual}")
    return (
        u,
        {
            "residual_rel": residual,
            "iterations": iters[0],
            "solve_time_s": time.perf_counter() - start,
        },
        preconditioner,
    )


def sample(axes, ids, u, point):
    """P1 evaluation; sort normalized coordinates to select the Kuhn tet."""
    idx, loc = [], []
    for axis, p in zip(axes, point):
        if p < axis[0] - 1e-15 or p > axis[-1] + 1e-15:
            raise ValueError("Point outside mesh")
        i = np.clip(np.searchsorted(axis, p, side="right") - 1, 0, len(axis) - 2)
        idx.append(int(i))
        loc.append((p - axis[i]) / (axis[i + 1] - axis[i]))
    order = np.argsort(-np.array(loc))
    weights = np.r_[1.0, np.array(loc)[order]] - np.r_[np.array(loc)[order], 0.0]
    vertices = [idx.copy()]
    for axis in order:
        idx = idx.copy()
        idx[axis] += 1
        vertices.append(idx)
    node_ids = [ids[v[2], v[0], v[1]] for v in vertices]
    if any(i < 0 and abs(w) > 1e-12 for i, w in zip(node_ids, weights)):
        raise ValueError("Point outside active package")
    return float(sum(w * u[i] for i, w in zip(node_ids, weights) if abs(w) > 1e-15))


def clip_polygon(poly, a, b, c):
    """Clip to a*x+b*y+c >= 0 in normalized box coordinates."""
    result = []
    for p, q in zip(poly, poly[1:] + poly[:1]):
        fp, fq = a * p[0] + b * p[1] + c, a * q[0] + b * q[1] + c
        if fp >= -1e-15:
            result.append(p)
        if (fp > 0 and fq < 0) or (fp < 0 and fq > 0):
            t = fp / (fp - fq)
            result.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return result


def cell_mean(axes, ids, u, bounds, zval):
    """Exact area integral of P1 field over a horizontal rectangular cell."""
    x, y, z = axes
    iz = int(np.clip(np.searchsorted(z, zval, side="right") - 1, 0, len(z) - 2))
    tz = (zval - z[iz]) / (z[iz + 1] - z[iz])
    coeff = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, tz)]
    integral, total_area = 0.0, 0.0
    for ix in range(
        max(0, np.searchsorted(x, bounds["x0_m"], side="right") - 1),
        min(len(x) - 1, np.searchsorted(x, bounds["x1_m"])),
    ):
        for iy in range(
            max(0, np.searchsorted(y, bounds["y0_m"], side="right") - 1),
            min(len(y) - 1, np.searchsorted(y, bounds["y1_m"])),
        ):
            dx, dy = x[ix + 1] - x[ix], y[iy + 1] - y[iy]
            a, b = (
                max(0.0, (bounds["x0_m"] - x[ix]) / dx),
                min(1.0, (bounds["x1_m"] - x[ix]) / dx),
            )
            c, d = (
                max(0.0, (bounds["y0_m"] - y[iy]) / dy),
                min(1.0, (bounds["y1_m"] - y[iy]) / dy),
            )
            for perm in PERMS:
                poly = [(a, c), (b, c), (b, d), (a, d)]
                for first, second in pairwise(perm):
                    aa, bb, cc = np.subtract(coeff[first], coeff[second])
                    poly = clip_polygon(poly, aa, bb, cc)
                    if len(poly) < 3:
                        break
                if len(poly) < 3:
                    continue
                p0 = np.array(poly[0])
                for j in range(1, len(poly) - 1):
                    p1, p2 = np.array(poly[j]), np.array(poly[j + 1])
                    cross = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (
                        p2[0] - p0[0]
                    )
                    ar = abs(cross) * dx * dy / 2
                    centre = (p0 + p1 + p2) / 3
                    if ar > 0:
                        val = sample(
                            axes,
                            ids,
                            u,
                            [x[ix] + centre[0] * dx, y[iy] + centre[1] * dy, zval],
                        )
                        integral += ar * val
                        total_area += ar
    area = (bounds["x1_m"] - bounds["x0_m"]) * (bounds["y1_m"] - bounds["y0_m"])
    if not np.isclose(total_area, area, rtol=1e-10, atol=1e-22):
        raise ValueError(f"Cell integration coverage: {total_area} != {area}")
    return float(integral / area)
