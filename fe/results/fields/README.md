# Finest primary temperature fields

`A-m2-primary.npz` and `B-m2-primary.npz` contain the complete level-2,
volumetric-source (primary) solutions underlying the die-mid-plane readout.
They cover the full model stack. Variant S is not included.

The arrays retain float64 precision. These are lossless NumPy NPZ archives with
ZIP LZMA compression, readable by `numpy.load(..., allow_pickle=False)` using
Python's standard `lzma` support. The two archives total 17,931,691 bytes, below
20,000,000 bytes. Their hashes are in [../SHA256SUMS](../SHA256SUMS).

| Array | Meaning |
|---|---|
| `x_m`, `y_m`, `z_m` | Increasing mesh-node coordinates, in metres. |
| `active_layer_map` | int8 cell map indexed `(iz, ix, iy)`, shape `(nz-1, nx-1, ny-1)`; layer index 0–7 in active cells, -1 outside the model. |
| `active_nodes` | Boolean node mask, shape `(nz, nx, ny)`. |
| `temperature_excess_K` | Float64 nodal vector `T - T_ref`, in the C-order traversal of `active_nodes == True`. No rounding or quantization. |
| `T_ref_K` | Scalar reference temperature; add it to the nodal vector to obtain absolute kelvin. |
| `case`, `level`, `source_variant` | Case A or B, refinement level 2, and `primary`. |
| `case_json_sha256`, `solver_sha256` | SHA256 of the input bytes and unchanged solver used to generate the field. |
| `reference_record_sha256` | SHA256 of the corresponding unchanged `../A-m2.json` or `../B-m2.json`, whose first record is the primary solution. |

Each active rectilinear cell is divided into the six Kuhn P1 tetrahedra used by
`fe/solver.py`: for each permutation of the x, y and z directions, follow the
vertices 000, the first unit step, the second unit step, and 111. Shared nodes
have one temperature. The layer map labels cells; nodes on layer interfaces
belong to all incident cells. The Boolean mask avoids storing redundant node
indices and makes the vector ordering explicit.

From the repository root, load the field for sampling with the existing solver:

```python
import numpy as np
from fe.solver import sample

with np.load("fe/results/fields/A-m2-primary.npz", allow_pickle=False) as packet:
    axes = tuple(packet[key] for key in ("x_m", "y_m", "z_m"))
    active = packet["active_nodes"]
    ids = np.full(active.shape, -1, dtype=np.int32)
    ids[active] = np.arange(active.sum(), dtype=np.int32)
    theta = packet["temperature_excess_K"]
    ambient = float(packet["T_ref_K"])

# point_m must lie in an active tetrahedron:
# T_K = sample(axes, ids, theta, point_m) + ambient
```

The excess-temperature vector preserves the actual solver unknown bit for bit,
including its original readout arithmetic. Retaining it instead of subtracting
the reference from an already rounded absolute-temperature vector avoids an
extra floating-point rounding step.

## Regeneration and verification

The exports used the pinned `fe/requirements.txt` environment, with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`, and the unchanged numerical routines:
`mesh_axes(data, case, 2)`, `layer_fields`, `assemble` with source layer 4 and the
JSON convection coefficient, and `solve_system(K, primary)`. The nodal map came
from `assemble`; the cell map is `where(fields[layer] > 0, layer, -1)` on each
z interval. Archives contain ordinary `.npy` members compressed with
`zipfile.ZIP_LZMA`; each member uses the fixed `ZipInfo` default date.

Before export, all 14 checked readout, mesh and solver result groups matched the
published primary JSON exactly: metric, winning die, paper difference, all
three-plane die values, cell means and selector, energy balance, field extrema,
nodes, degrees of freedom, elements, layer counts, residual and iterations.
Every archive array then passed exact equality after loading with NumPy.
Loading the final archives and recomputing their readouts also reproduced the
published values exactly. No tolerance was used for these comparisons.

The full driver was separately rerun in scratch space: all six meshes and eight
solutions matched every published scientific, mesh and solver value, including
both convergence records. Runtime and memory observations differed and were
excluded from that comparison. The original result JSON, report and plot remain
byte for byte unchanged; the script hashes inside the original `study.json`
still identify its original generation. The enclosing SHA256 manifest now also
records the revised driver, tests, README, code licence, and these field packets.
