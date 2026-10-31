"""Generate the human-readable tables and convergence plot from study.json."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def write_report(study, output):
    output = Path(output)
    lines = [
        "# Independent steady-state FE comparison",
        "",
        (
            "The input is the configured TAP-2.5D package, at reference commit "
            "`af66bcfc21220ce3d8e92e89273979241e54b04a`. Temperatures below are computed "
            "with conforming P1 tetrahedral FE. The primary metric is the maximum die-centre "
            "temperature at the chip mid-plane with uniform volumetric power."
        ),
        "",
        (
            "The paper uses a different HotSpot block/cell selector. Differences here are "
            "numerical comparisons of the configured model, not silicon validation."
        ),
        "",
        "| Case | Primary M (°C) | Die | Paper (°C) | M − paper (K) | S mid-plane M (°C) |",
        "|---|---:|---|---:|---:|---:|",
    ]
    for cid, case in study["cases"].items():
        fine, sheet = case["primary"][-1], case["variant_S"]
        lines.append(
            f"| {cid} | {fine['M_C']:.6f} | {fine['M_die']} | {case['paper_target_C']:.2f} | "
            f"{fine['M_minus_paper_K']:+.6f} | {sheet['M_C']:.6f} |"
        )
    lines += [
        "",
        "## Method and reproducibility",
        "",
        (
            "The weak form is ∫ k ∇θ·∇v dV + ∫top h θ v dA = ∫ q v dV, with θ=T−318.15 K. "
            "Variant S replaces the volume term with a single internal die-bottom sheet load. "
            "Shared nodes impose perfect contact; all other exterior faces are natural zero-flux boundaries. "
            "Exact element integration is used for diffusion, volume loads, sheet loads and the consistent Robin matrix."
        ),
        "",
        (
            "Every box is split into six Kuhn tetrahedra. All block edges and all layer interfaces conform "
            "to the mesh. Floating duplicates within 2e−16 m are merged; the 5.1 µm A gap is retained "
            "and filled with k=1.6 W/(m K). Meshes halve every x/y/z interval, including thin layers "
            "and halos. The listed minimum size is the shortest edge, and the maximum is the largest "
            "tetrahedron diameter. Centre points are nodes at every level. Horizontal cell means are "
            "integrated exactly by clipping the P1 pieces against the requested rectangles."
        ),
        "",
        (
            "Classical Ruge–Stuben AMG preconditions conjugate gradients (relative tolerance 1e−10). "
            "The reported residual is ||Kθ−f||₂/||f||₂. Assembly and primary solver times include AMG setup "
            "where applicable; S reuses the same matrix and preconditioner. Peak RSS is measured in "
            "a fresh process per case/mesh; the S measurement includes the primary run high-water mark."
        ),
        "",
        (
            "Run `python fe/run_study.py` with the pinned dependencies in `fe/requirements.txt`. "
            "It verifies the original four-file manifest, regenerates all six meshes and eight solutions, "
            "and writes these tables, the plot, JSON outputs and an output SHA256 manifest. "
            "Timing/RSS and platform strings vary by machine; numeric fields should reproduce within solver tolerance. "
            "The full command and validation instructions are in [../README.md](../README.md)."
        ),
        "",
        f"Input SHA256: `{study['case_json_sha256']}`. Script SHA256: `{study['script_sha256']}`.",
        "",
        "## Convergence",
        "",
        (
            "The estimate uses the last three nested meshes: p=log₂(|Δcoarse/Δfine|), "
            "M∞=Mfine+Δfine/(2ᵖ−1), and absolute GCI=1.25|M∞−Mfine|. "
            "It is an asymptotic estimate, not a rigorous bound. No estimate is accepted for oscillating "
            "differences or a changing maximizing die."
        ),
        "",
        "| Case | Observed p | Extrapolated M (°C) | Absolute GCI (K) | ≤0.1 K aim |",
        "|---|---:|---:|---:|---|",
    ]
    for cid, case in study["cases"].items():
        c = case["convergence"]
        if c["status"] == "estimated":
            lines.append(
                f"| {cid} | {c['observed_order']:.4f} | {c['M_extrapolated_C']:.6f} | "
                f"{c['GCI_absolute_K']:.6f} | {'met' if c['accuracy_goal_met'] else 'not met'} |"
            )
        else:
            lines.append(f"| {cid} | inconclusive | — | — | not established |")
    lines += ["", "![Mesh convergence](convergence.svg)", ""]
    for cid, case in study["cases"].items():
        lines += [
            f"## Case {cid}",
            "",
            "| Mesh | P1 tetrahedra | DOFs | Min edge (m) | Max diameter (m) | M (°C) | Residual | Solver s | Peak MiB |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for r in case["primary"] + [case["variant_S"]]:
            lines.append(
                f"| m{r['level']} {r['source_variant']} | {r['elements']:,} | {r['dofs']:,} | "
                f"{r['h_min_edge_m']:.4g} | {r['h_max_diameter_m']:.4g} | {r['M_C']:.6f} | "
                f"{r['residual_rel']:.3g} | {r['solve_time_s']:.3f} | {r['peak_memory_MiB']:.1f} |"
            )
        lines += [
            "",
            "Through-layer interval counts (L0…L7): "
            + "; ".join(
                "m"
                + str(r["level"])
                + " = "
                + str(list(r["elements_through_layer"].values()))
                for r in case["primary"]
            )
            + ".",
            "",
            (
                "Finest-mesh die temperatures below are °C; `study.json` also contains Kelvin values, "
                "every mesh, all L0–L5 block cell averages, and extrema with coordinates."
            ),
            "",
            "| Source | Die | Bottom | Mid | Top | L4 HotSpot-cell mean |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for r in [case["primary"][-1], case["variant_S"]]:
            for die, temp in r["T_i"].items():
                t = temp["C"]
                lines.append(
                    f"| {r['source_variant']} | {die} | {t['bottom']:.6f} | {t['mid']:.6f} | "
                    f"{t['top']:.6f} | {r['T_cell']['L4'][die]['C']:.6f} |"
                )
        lines += [
            "",
            "| Source | Plane | Max die-centre (°C) | Difference from paper (K) |",
            "|---|---|---:|---:|",
        ]
        for r in [case["primary"][-1], case["variant_S"]]:
            for plane in ["bottom", "mid", "top"]:
                m = max(t["C"][plane] for t in r["T_i"].values())
                lines.append(
                    f"| {r['source_variant']} | {plane} | {m:.6f} | {m - case['paper_target_C']:+.6f} |"
                )
        lines += [
            "",
            "| Source | Input W | Convective W | Relative imbalance | Mean top excess K |",
            "|---|---:|---:|---:|---:|",
        ]
        for r in [case["primary"][-1], case["variant_S"]]:
            e = r["energy_balance"]
            lines.append(
                f"| {r['source_variant']} | {e['P_in_W']:.9f} | {e['Q_out_W']:.9f} | "
                f"{e['rel_error']:.3g} | {e['mean_excess_convective_face_K']:.9f} |"
            )
        lines += [
            "",
            "| Source | Global min (°C), (x,y,z) m | Global max (°C), (x,y,z) m | L3/L4 cell selector (°C) |",
            "|---|---|---|---|",
        ]
        for r in [case["primary"][-1], case["variant_S"]]:
            mn, mx = r["field_extremes"]["min"], r["field_extremes"]["max"]
            s = r["cell_selector_L3_L4"]
            fmt = lambda point: "(" + ", ".join(f"{v:.9g}" for v in point) + ")"
            lines.append(
                f"| {r['source_variant']} | {mn['T_C']:.6f}, {fmt(mn['at_m'])} | "
                f"{mx['T_C']:.6f}, {fmt(mx['at_m'])} | {s['C']:.6f}, {s['layer']}/{s['block']} |"
            )
        lines += [
            "",
            "| Layer | Primary field max °C | S field max °C |",
            "|---|---:|---:|",
        ]
        for lid in case["primary"][-1]["field_extremes"]["layer_max"]:
            lines.append(
                f"| {lid} | {case['primary'][-1]['field_extremes']['layer_max'][lid]['T_C']:.6f} | "
                f"{case['variant_S']['field_extremes']['layer_max'][lid]['T_C']:.6f} |"
            )
        lines.append("")
    lines += ["## Environment and limits", ""]
    lines += [f"- {key}: `{value}`" for key, value in study["versions"].items()]
    lines += [
        "",
        (
            "scikit-fem provides the independent assembly oracle in the tests; production assembly uses "
            "the algebraically identical compressed P1 element contributions. The mesher is the versioned "
            "Python implementation in `fe/solver.py`."
        ),
        "",
        (
            "The GCI estimates concern the primary centre metric. S is computed on the finest mesh only; "
            "it has no separate convergence estimate. Field extrema and cell selectors are different metrics "
            "and do not inherit the primary metric error estimate. The paper solver tolerance/build, exact "
            "layout-to-figure fidelity and the full HotSpot package-node selector are not reproduced. "
            "Decoded coordinate uncertainty is not propagated. No comparison to withheld evaluator values is made."
        ),
        "",
        "## Attribution",
        "",
        (
            "Inputs: SentiLab reference case (CC BY 4.0), preserving its LICENSE and provenance. "
            "Paper: Y. Ma et al., “TAP-2.5D: A Thermally-Aware Chiplet Placement Methodology for 2.5D Systems”, "
            "DATE 2021, pp. 1246–1251, doi:10.23919/DATE51398.2021.9474011. "
            "TAP-2.5D and HotSpot retain their own licences; no source from either solver is vendored."
        ),
        "",
    ]
    (output / "REPORT.md").write_text("\n".join(lines))
    plt.rcParams["svg.hashsalt"] = "independent-fe-convergence"
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.5), layout="constrained")
    for ax, (cid, case) in zip(axs, study["cases"].items()):
        records = case["primary"]
        levels = [r["level"] for r in records]
        ax.plot(levels, [r["M_C"] for r in records], "o-", label="Primary FE")
        c = case["convergence"]
        if c["status"] == "estimated":
            ax.axhline(
                c["M_extrapolated_C"],
                color="tab:orange",
                linestyle="--",
                label="Richardson estimate",
            )
        ax.set(
            title=f"Case {cid}",
            xlabel="Uniform refinement level",
            ylabel="Centre metric M (°C)",
            xticks=levels,
        )
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.savefig(output / "convergence.svg", metadata={"Date": None})
    plt.close(fig)
    svg = output / "convergence.svg"
    svg.write_text(
        "\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n"
    )
