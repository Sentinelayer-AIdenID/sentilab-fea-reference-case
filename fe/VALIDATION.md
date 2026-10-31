# Validation receipt

Run on 2026-10-31 with Python 3.12.3 and the dependencies pinned in
`fe/requirements.txt`. Commands below were run from the repository root in the
installed virtual environment. All final gates exited 0.

| Gate command | Observed result |
|---|---|
| `python -m compileall -q fe` | Python sources compile |
| `ruff check fe` | All checks passed |
| `python -m pytest -q fe/test_solver.py` | 24 passed |
| `python -m pip check` | No broken requirements found |
| `python fe/run_study.py` | Six meshes, eight solutions; complete report and plot generated |
| `sha256sum -c SHA256SUMS` | Four original files verified |
| `(cd fe/results && sha256sum -c SHA256SUMS)` | Sixteen input/code/test/document/output entries verified |
| `python -m pip_audit -r fe/requirements.txt` | No known vulnerabilities found |

The final rebuild reproduced both primary sequences and both finest S metrics
from the earlier run. The maximum measured relative linear residual over all
solutions is below 1e-10; the absolute relative energy imbalance is below 2e-11.
The observed orders and GCI estimates are reported in `results/REPORT.md` and
`results/study.json`. The convergence plot was visually inspected.

For base-failure evidence, a separate detached worktree was created at
`af66bcfc21220ce3d8e92e89273979241e54b04a`. Only the new test file was copied
into it. `python -m pytest -q fe/test_solver.py` exited 2 during collection:
`ModuleNotFoundError: No module named 'fe.run_study'`. The reference base contains
no FE implementation. This is a new-feature availability failure, not evidence
of a pre-existing numerical defect. On the implementation, all 24 tests pass,
including an independent scikit-fem numerical assembly oracle and analytic
solutions.

The paper PDF was downloaded from its public HTTPS source. Its measured SHA256
is `b209f51aec188e2bcdbf6f4acc3723fb30fbc5c14cca7e70141d478edadec234`, matching
the input provenance. No source parameters were fitted or changed.

This receipt documents local execution, not external review or publication.
The original data-only repository has no Node manifest, TypeScript project or
npm build. No such files were introduced. The Python/scientific gates above are
the implementation's applicable gates.
