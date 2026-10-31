#!/usr/bin/env python3
"""Rebuild the independent FE study from case.json in isolated worker processes."""

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# Support both `python fe/run_study.py` and module invocation.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fe.solver import assemble, cell_mean, layer_fields, mesh_axes, sample, solve_system

ROOT = Path(__file__).resolve().parent.parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_inputs(path):
    for line in (path.parent / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        if digest(path.parent / name.lstrip("*")) != expected:
            raise ValueError("Reference manifest mismatch: " + name)


def environment():
    cpu = platform.processor() or platform.machine()
    if Path("/proc/cpuinfo").exists():
        cpu = next(
            (
                line.split(":", 1)[1].strip()
                for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")
            ),
            cpu,
        )
    return {
        "python": platform.python_version(),
        **{
            p: importlib.metadata.version(p)
            for p in ["numpy", "scipy", "scikit-fem", "pyamg"]
        },
        "mesher": "rectilinear Kuhn tetrahedra v1 (fe/solver.py)",
        "os": platform.platform(),
        "cpu": cpu,
        "architecture": platform.machine(),
        "worker_threads": {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
    }


def extrema(axes, ids, u, ambient, stack):
    active = ids >= 0
    values = np.full(ids.shape, np.nan)
    values[active] = u[ids[active]] + ambient

    def at(flat):
        iz, ix, iy = np.unravel_index(flat, values.shape)
        z = axes[2][iz]
        return {
            "T_K": float(values[iz, ix, iy]),
            "T_C": float(values[iz, ix, iy] - 273.15),
            "at_m": [float(axes[0][ix]), float(axes[1][iy]), float(z)],
            "layers": [
                s["index"]
                for s in stack
                if s["z_bottom_m"] - 1e-15 <= z <= s["z_top_m"] + 1e-15
                and s["extent_m"]["x0"] - 1e-15
                <= axes[0][ix]
                <= s["extent_m"]["x1"] + 1e-15
                and s["extent_m"]["y0"] - 1e-15
                <= axes[1][iy]
                <= s["extent_m"]["y1"] + 1e-15
            ],
        }

    layers = {}
    for layer in stack:
        e = layer["extent_m"]
        mask = (
            (axes[2][:, None, None] >= layer["z_bottom_m"] - 1e-15)
            & (axes[2][:, None, None] <= layer["z_top_m"] + 1e-15)
            & (axes[0][None, :, None] >= e["x0"] - 1e-15)
            & (axes[0][None, :, None] <= e["x1"] + 1e-15)
            & (axes[1][None, None, :] >= e["y0"] - 1e-15)
            & (axes[1][None, None, :] <= e["y1"] + 1e-15)
        )
        layers["L" + str(layer["index"])] = at(
            np.nanargmax(np.where(mask, values, np.nan))
        )
    return {
        "max": at(np.nanargmax(values)),
        "min": at(np.nanargmin(values)),
        "layer_max": layers,
    }


def readouts(data, cid, axes, ids, u, rhs, robin):
    ambient = data["boundary_conditions"]["convective_face"]["T_ref_K"]
    case = data["cases"][cid]
    temps = {}
    for die in case["chiplets"]:
        tk = {
            plane: sample(
                axes,
                ids,
                u,
                [
                    die["x_centre_m"],
                    die["y_centre_m"],
                    data["metric"]["z_die_" + plane + "_m"],
                ],
            )
            + ambient
            for plane in ["bottom", "mid", "top"]
        }
        temps[die["flp_name"]] = {"K": tk, "C": {k: v - 273.15 for k, v in tk.items()}}
    winner = max(temps, key=lambda d: temps[d]["K"]["mid"])
    cells = {}
    for lid, blocks in case["layers"].items():
        layer = data["stack"][int(lid[1:])]
        z = (layer["z_bottom_m"] + layer["z_top_m"]) / 2
        cache = {}
        cells[lid] = {}
        for block in blocks:
            value = 0.0
            for cell in block["hotspot_centre_cells_64"]:
                key = (cell["row_i"], cell["col_j"])
                if key not in cache:
                    cache[key] = cell_mean(axes, ids, u, cell, z) + ambient
                value += cell["weight"] * cache[key]
            cells[lid][block["name"]] = {"K": value, "C": value - 273.15}
    selector = max(
        ((lid, name, val) for lid in ["L3", "L4"] for name, val in cells[lid].items()),
        key=lambda item: item[2]["K"],
    )
    pin, pout = float(rhs.sum()), float(robin @ u)
    metric = temps[winner]["C"]["mid"]
    return {
        "M_C": metric,
        "M_die": winner,
        "M_minus_paper_K": metric - case["paper_target_C"],
        "T_i": temps,
        "T_cell": cells,
        "cell_selector_L3_L4": {
            "layer": selector[0],
            "block": selector[1],
            **selector[2],
        },
        "energy_balance": {
            "P_in_W": pin,
            "Q_out_W": pout,
            "rel_error": (pout - pin) / pin,
            "mean_excess_convective_face_K": pout / float(robin.sum()),
            "analytic_mean_excess_K": case["closed_form_checks"][
                "mean_excess_on_convective_face_K"
            ],
        },
        "field_extremes": extrema(axes, ids, u, ambient, data["stack"]),
    }


def worker(data, cid, level, include_sheet):
    start = time.perf_counter()
    axes, lids, counts = mesh_axes(data, cid, level)
    fields, q, flux = layer_fields(data, cid, axes)
    K, primary, sheet, robin, ids, elements = assemble(
        axes,
        lids,
        fields,
        q,
        flux,
        4,
        data["boundary_conditions"]["convective_face"]["h_W_per_m2K"],
    )
    assembly_time = time.perf_counter() - start
    minimum = min(float(np.min(np.diff(a))) for a in axes)
    # All tetrahedra share the box body diagonal (their diameter).
    max_diameter = max(
        float(
            np.sqrt(
                np.max(np.diff(axes[0])[np.any(f > 0, axis=1)]) ** 2
                + np.max(np.diff(axes[1])[np.any(f > 0, axis=0)]) ** 2
                + np.max(np.diff(axes[2])[lids == lid]) ** 2
            )
        )
        for lid, f in enumerate(fields)
    )
    preconditioner = None
    records = []
    for variant, rhs in [("primary", primary)] + (
        [("S", sheet)] if include_sheet else []
    ):
        u, stats, preconditioner = solve_system(K, rhs, preconditioner)
        result = {
            "case": cid,
            "source_variant": variant,
            "level": level,
            "refinement_factor": 2**level,
            "element": "P1 tetrahedron",
            "element_order": 1,
            "elements": elements,
            "nodes": K.shape[0],
            "dofs": K.shape[0],
            "h_min_edge_m": minimum,
            "h_max_diameter_m": max_diameter,
            "elements_through_layer": {
                "L" + str(i): int(n) for i, n in enumerate(counts)
            },
            "assembly_time_s": assembly_time,
            **stats,
            **readouts(data, cid, axes, ids, u, rhs, robin),
        }
        result["process_elapsed_s"] = time.perf_counter() - start
        result["assembly_plus_solve_time_s"] = assembly_time + stats["solve_time_s"]
        result["peak_memory_MiB"] = (
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        )
        result["memory_scope"] = (
            "isolated case/mesh process high-water RSS; S shares primary assembly and AMG setup"
        )
        records.append(result)
        print(
            f"{cid} m{level} {variant}: M={result['M_C']:.8f} C; DOFs={K.shape[0]}; residual={stats['residual_rel']:.3g}; RSS={result['peak_memory_MiB']:.1f} MiB",
            flush=True,
        )
    return records


def convergence(records):
    a, b, c = records[-3:]
    d1, d2 = b["M_C"] - a["M_C"], c["M_C"] - b["M_C"]
    ratio = b["refinement_factor"] / a["refinement_factor"]
    if d1 * d2 <= 0 or abs(d2) >= abs(d1) or len({v["M_die"] for v in [a, b, c]}) != 1:
        return {
            "status": "inconclusive",
            "reason": "non-asymptotic differences or changing maximizing die",
            "finest_change_K": abs(d2),
        }
    order = float(np.log(abs(d1 / d2)) / np.log(ratio))
    correction = d2 / (ratio**order - 1)
    error = 1.25 * abs(correction)
    return {
        "status": "estimated",
        "observed_order": order,
        "refinement_ratio": ratio,
        "M_extrapolated_C": c["M_C"] + correction,
        "Richardson_correction_K": correction,
        "GCI_absolute_K": error,
        "GCI_safety_factor": 1.25,
        "accuracy_goal_met": error <= 0.1,
        "caveat": "Three-level asymptotic estimate, not a rigorous error bound.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-json", type=Path, default=ROOT / "case.json")
    parser.add_argument("--output", type=Path, default=ROOT / "fe/results")
    parser.add_argument("--levels", type=int, nargs="+", default=[0, 1, 2])
    parser.add_argument("--worker-case", choices=["A", "B"])
    parser.add_argument("--sheet", action="store_true")
    args = parser.parse_args()
    verify_inputs(args.case_json)
    data = json.loads(args.case_json.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    if args.worker_case:
        records = worker(data, args.worker_case, args.levels[0], args.sheet)
        name = f"{args.worker_case}-m{args.levels[0]}.json"
        (args.output / name).write_text(
            json.dumps(records, indent=2, allow_nan=False) + "\n"
        )
        return
    if len(args.levels) < 3 or any(
        b != a + 1 for a, b in zip(args.levels[:-1], args.levels[1:])
    ):
        parser.error("Provide at least three successive uniform refinement levels")
    for cid in ["A", "B"]:
        for level in args.levels:
            cmd = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--case-json",
                str(args.case_json.resolve()),
                "--output",
                str(args.output.resolve()),
                "--worker-case",
                cid,
                "--levels",
                str(level),
            ]
            if level == args.levels[-1]:
                cmd.append("--sheet")
            env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
            subprocess.run(cmd, env=env, check=True)
    cases = {}
    for cid in ["A", "B"]:
        runs = [
            json.loads((args.output / f"{cid}-m{level}.json").read_text())
            for level in args.levels
        ]
        cases[cid] = {
            "paper_target_C": data["cases"][cid]["paper_target_C"],
            "primary": [r[0] for r in runs],
            "variant_S": runs[-1][1],
            "convergence": convergence([r[0] for r in runs]),
        }
    report = {
        "schema": "independent-tap25d-fe-v1",
        "case_json_sha256": digest(args.case_json),
        "script_sha256": digest(__file__),
        "solver_sha256": digest(ROOT / "fe/solver.py"),
        "report_script_sha256": digest(ROOT / "fe/report.py"),
        "versions": environment(),
        "cases": cases,
    }
    (args.output / "study.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps({c: cases[c]["convergence"] for c in cases}, indent=2))
    from fe.report import write_report

    write_report(report, args.output)
    # Input manifests remain unchanged; this manifest covers implementation and output.
    paths = [
        args.case_json,
        Path(__file__).resolve(),
        ROOT / "fe/solver.py",
        ROOT / "fe/requirements.txt",
        ROOT / "fe/report.py",
        ROOT / "fe/test_solver.py",
        ROOT / "fe/README.md",
    ]
    paths += sorted(args.output.glob("*.json")) + [
        args.output / "REPORT.md",
        args.output / "convergence.svg",
    ]
    (args.output / "SHA256SUMS").write_text(
        "".join(
            digest(p) + "  " + os.path.relpath(p, args.output) + "\n" for p in paths
        )
    )


if __name__ == "__main__":
    main()
