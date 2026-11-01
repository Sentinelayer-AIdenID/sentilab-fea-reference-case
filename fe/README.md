# Independent 3D steady-state finite-element study

The Python code in `fe/` is [MIT-licensed](LICENSE); the reference input data remain [CC BY 4.0](../LICENSE).

This directory adds a reproducible numerical solution for both layouts in the
unchanged, pinned `case.json`. Read [results/REPORT.md](results/REPORT.md) for the
computed comparison and [results/study.json](results/study.json) for every mesh,
three-plane die-centre readout, L0–L5 cell mean, extrema and conservation result.
The original README, schema, input JSON, licence and manifest are unchanged.

## Rebuild

Python 3.12 was used. From the repository root:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r fe/requirements.txt
python fe/run_study.py
(cd fe/results && sha256sum -c SHA256SUMS)
```

`run_study.py` verifies the original input manifest, constructs all geometry and
materials from `case.json`, runs three nested meshes per case, solves the finest
mesh again with variant S, and regenerates the JSON, Markdown report, SVG plot,
and SHA256 manifest. It runs one case/mesh process at a time and limits numerical
libraries to one thread. The measured peak process RSS was below 1 GiB in the
recorded run; allow 2 GiB available RAM. No external mesher or proprietary solver
is needed. The algorithm does not call TAP-2.5D or HotSpot.

The output manifest binds the input, implementation, dependencies, tests, this
README and outputs. It is relative to `fe/results`. Timing, memory and platform
strings are observations and will vary on reproduction, so regenerated output
hashes can differ. Temperature comparisons should use the reported solver and
mesh tolerances rather than output-byte equality. The SVG contains no timestamp.
The original input manifest is always checked byte-for-byte.
`--case-json` must select the manifest-verified `case.json`; unlisted custom inputs
are refused. The driver parses the verified bytes, binds their resolved path and
SHA256 before launching workers, and requires each worker to verify that same hash.

## Scientific implementation

`solver.py` partitions each rectilinear box into six conforming P1 tetrahedra.
All physical block and layer boundaries are mesh-aligned. Every axis interval is
bisected at each refinement, so the factor two applies in-plane and through every
layer. Only duplicate coordinates within 2e-16 m are merged; the narrow declared
underfill strip is retained. Exact element integration admits a compact
axis-edge diffusion assembly without retaining the full tetrahedral connectivity.
This is a finite-element stiffness matrix, independently compared with explicit
scikit-fem element assembly in the tests.

The solve is for excess temperature relative to the JSON Robin reference.
Volumetric and sheet source terms use consistent FE loads. Variant S is assembled
once on the shared internal die-bottom plane, preserving thermal continuity.
The entire sink top receives a consistent Robin boundary matrix, and every other
exterior boundary remains adiabatic. Conductivity uses each block's JSON value;
stack defaults do not replace block overrides. No parameters are fitted to the
paper values.

Centre values use the actual P1 field. Cell means integrate the P1 pieces exactly
on the horizontal rectangle at the layer mid-plane; they are not samples at a
cell centroid. Global P1 extrema occur at nodes. Shared-interface extrema are
attributed to all incident layers containing that point. Three-mesh Richardson
extrapolation and absolute GCI with safety factor 1.25 estimate the primary
metric's discretization error. A changing hottest die or non-asymptotic sequence
returns an inconclusive estimate rather than an accuracy claim.

## Validation gates

The pinned repository was data-only, with no Node manifest or TypeScript/build
configuration. These Python and scientific gates apply to this implementation;
no artificial npm project is added to run unrelated commands.

```sh
python -m compileall -q fe
ruff check fe
python -m pytest -q fe/test_solver.py
python -m pip check
python fe/run_study.py
sha256sum -c SHA256SUMS
(cd fe/results && sha256sum -c SHA256SUMS)
python -m pip_audit -r fe/requirements.txt
```

The tests compare heterogeneous diffusion, volume loads, internal sheet loads and
Robin terms against independent scikit-fem assembly; check exact layered flux
solutions and convergence to a volumetrically heated analytic slab; integrate
affine and piecewise-linear fields; verify geometric volumes, both exterior area
classes, power, material overrides, the gap and nested meshes; and exercise the
convergence reporting guards. Production assembly uses NumPy/SciPy/PyAMG;
scikit-fem is the independent test oracle. See `results/study.json` for versions,
iterations, measured residuals, wall time and memory.

## Scope and attribution

This is the configured model from the published layout, not a physical package,
CAD reconstruction, measured temperature or manufacturer validation. The paper
selector includes compact-model block/cell and package-node values with no unique
continuum point. We report the prescribed primary continuum metric, its three
planes and variant S, and a separate L3/L4 cell-selector analogue. Withheld
reference-evaluator values are neither used nor inferred. Layout uncertainty is
not propagated; the primary GCI does not certify field extrema or variant S.

The SentiLab input data retain CC BY 4.0 and their source attribution in the
repository LICENSE and README. Ma et al., TAP-2.5D, DATE 2021,
doi:10.23919/DATE51398.2021.9474011 is the paper source. TAP-2.5D and HotSpot retain
their own licences; neither codebase is copied into this implementation.
