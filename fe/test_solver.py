"""Numerical checks against independent assembly and analytic solutions."""

import hashlib
import json
import shutil
import sys
from itertools import permutations
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose
from skfem import (
    Basis,
    BilinearForm,
    ElementTetP1,
    FacetBasis,
    InteriorFacetBasis,
    LinearForm,
    MeshTet,
    asm,
)
from skfem.helpers import dot, grad

from fe import run_study
from fe.run_study import convergence, verify_inputs
from fe.solver import (
    assemble,
    cell_mean,
    layer_fields,
    mesh_axes,
    sample,
    solve_system,
)


def nodal_grid(axes):
    x, y, z = axes
    zz, xx, yy = np.meshgrid(z, x, y, indexing="ij")
    ids = np.arange(xx.size).reshape(xx.shape)
    return ids, np.array([xx.ravel(), yy.ravel(), zz.ravel()])


def independent_tetrahedra(axes):
    """Explicit connectivity, independent of the compressed edge assembly."""
    ids, points = nodal_grid(axes)
    cells, boxes = [], []
    nx, ny, nz = (len(a) - 1 for a in axes)
    for iz in range(nz):
        for ix in range(nx):
            for iy in range(ny):
                for order in permutations(range(3)):
                    vertex = [ix, iy, iz]
                    tet = [ids[iz, ix, iy]]
                    for direction in order:
                        vertex[direction] += 1
                        i, j, k = vertex
                        tet.append(ids[k, i, j])
                    cells.append(tet)
                    boxes.append((ix, iy, iz))
    return MeshTet(points, np.array(cells).T), np.array(boxes)


def test_assembly_matches_scikit_fem_on_heterogeneous_irregular_mesh():
    axes = (
        np.array([-0.4, 0.0, 0.7]),
        np.array([0.1, 0.35, 1.2]),
        np.array([0.0, 0.2, 0.55, 0.9]),
    )
    layers = np.array([0, 1, 1])
    fields = [np.array([[2.0, 6.0], [4.0, 3.0]]), np.array([[9.0, 1.5], [7.0, 5.0]])]
    q = np.array([[11.0, 7.0], [3.0, 9.0]])
    flux = np.array([[8.0, 4.0], [5.0, 2.0]])
    h = 3.7
    matrix, load, sheet, robin, ids, count = assemble(
        axes, layers, fields, q, flux, 1, h
    )
    mesh, boxes = independent_tetrahedra(axes)
    ix, iy, iz = boxes.T
    conductivity = np.array([fields[layers[k]][i, j] for i, j, k in boxes])
    source = np.where(layers[iz] == 1, q[ix, iy], 0.0)
    basis = Basis(mesh, ElementTetP1(), intorder=2)

    @BilinearForm
    def diffusion(u, v, w):
        return w.k * dot(grad(u), grad(v))

    @BilinearForm
    def convection(u, v, w):
        return h * u * v

    @LinearForm
    def volume_source(v, w):
        return w.q * v

    @LinearForm
    def face_source(v, w):
        i = np.searchsorted(axes[0], w.x[0], side="right") - 1
        j = np.searchsorted(axes[1], w.x[1], side="right") - 1
        return flux[i, j] * v

    @LinearForm
    def robin_weights(v, w):
        return h * v

    top = mesh.facets_satisfying(lambda p: np.isclose(p[2], axes[2][-1]))
    source_face = mesh.facets_satisfying(lambda p: np.isclose(p[2], 0.2))
    top_basis = FacetBasis(mesh, ElementTetP1(), facets=top, intorder=2)
    sheet_basis = InteriorFacetBasis(
        mesh, ElementTetP1(), facets=source_face, intorder=2
    )
    reference = asm(diffusion, basis, k=conductivity[:, None]) + asm(
        convection, top_basis
    )
    assert np.array_equal(ids.ravel(), np.arange(basis.N))
    assert count == mesh.nelements
    assert_allclose(matrix.toarray(), reference.toarray(), rtol=2e-13, atol=2e-13)
    assert_allclose(load, asm(volume_source, basis, q=source[:, None]), atol=2e-14)
    assert_allclose(sheet, asm(face_source, sheet_basis), atol=2e-14)
    assert_allclose(robin, asm(robin_weights, top_basis), atol=2e-14)
    assert_allclose(matrix @ np.ones(basis.N), robin, atol=2e-13)
    assert_allclose(matrix.toarray(), matrix.toarray().T, atol=2e-14)
    assert np.linalg.eigvalsh(matrix.toarray())[0] > 0.0


@pytest.mark.parametrize("source_layer", [0, 1])
def test_uniform_sheet_flux_matches_exact_layered_slab(source_layer):
    axes = (
        np.array([0.0, 0.3, 1.0]),
        np.array([0.0, 0.2, 0.7]),
        np.array([0.0, 0.1, 0.3, 0.55, 0.8, 1.0]),
    )
    layers = np.array([0, 0, 1, 1, 1])
    fields = [np.full((2, 2), 2.0), np.full((2, 2), 5.0)]
    zeros = np.zeros((2, 2))
    matrix, _, sheet, robin, ids, _ = assemble(
        axes, layers, fields, zeros, np.full((2, 2), 7.0), source_layer, 3.0
    )
    u, stats, _ = solve_system(matrix, sheet)
    source_z = 0.0 if source_layer == 0 else 0.3

    def exact(z):
        z = np.maximum(z, source_z)
        return 7.0 / 3.0 + 7.0 * (
            np.maximum(0.3 - z, 0.0) / 2.0 + (1.0 - np.maximum(z, 0.3)) / 5.0
        )

    _, points = nodal_grid(axes)
    assert_allclose(u, exact(points[2]), rtol=2e-10, atol=2e-10)
    assert_allclose(sample(axes, ids, u, [0.24, 0.46, 0.19]), exact(0.19), atol=2e-10)
    assert_allclose(robin @ u, 7.0 * 0.7, rtol=2e-10)
    assert_allclose(sheet.sum(), 7.0 * 0.7, rtol=2e-14)
    assert stats["residual_rel"] < 2e-9


def test_volumetric_slab_converges_to_analytic_parabola():
    errors = []
    for n in [2, 4, 8]:
        axis = np.linspace(0.0, 1.0, n + 1)
        axes = (axis, axis, axis)
        shape = (n, n)
        matrix, load, _, robin, _, _ = assemble(
            axes,
            np.zeros(n, dtype=int),
            [np.full(shape, 2.0)],
            np.full(shape, 3.0),
            np.zeros(shape),
            0,
            4.0,
        )
        u, _, _ = solve_system(matrix, load)
        _, points = nodal_grid(axes)
        # -2 T'' = 3, T'(0)=0, -2 T'(1)=4 theta(1).
        exact = 0.75 + 0.75 * (1.0 - points[2] ** 2)
        errors.append(np.max(np.abs(u - exact)))
        assert_allclose(load.sum(), 3.0, rtol=2e-14)
        assert_allclose(robin @ u, 3.0, rtol=2e-9)
    # P1 interpolation and the consistent source converge to the continuum
    # solution; a source at the bottom face instead would not pass this test.
    assert errors[1] < errors[0] / 2.0
    assert errors[2] < errors[1] / 2.0
    assert errors[2] < 0.01


@pytest.mark.parametrize("zval", [0.0, 0.11, 0.2, 0.41, 0.55, 1.0])
def test_cell_mean_integrates_affine_field_over_irregular_mesh(zval):
    axes = (
        np.array([0.0, 0.13, 0.61, 1.0]),
        np.array([0.0, 0.27, 0.8, 1.0]),
        np.array([0.0, 0.2, 0.55, 1.0]),
    )
    ids, points = nodal_grid(axes)
    u = 17.0 + np.array([2.0, -3.0, 4.0]) @ points
    bounds = {"x0_m": 0.07, "x1_m": 0.93, "y0_m": 0.14, "y1_m": 0.87}
    expected = 17.0 + 2.0 * 0.5 - 3.0 * 0.505 + 4.0 * zval
    assert_allclose(cell_mean(axes, ids, u, bounds, zval), expected, atol=2e-12)
    assert_allclose(
        sample(axes, ids, u, [0.22, 0.62, zval]),
        17.0 + 2.0 * 0.22 - 3.0 * 0.62 + 4.0 * zval,
        atol=2e-13,
    )


@pytest.mark.parametrize("zval", [0.0, 0.2, 0.5, 0.93, 1.0])
def test_cell_mean_integrates_piecewise_linear_kuhn_field(zval):
    axes = tuple(np.array([0.0, 1.0]) for _ in range(3))
    ids, points = nodal_grid(axes)
    # max(x,y,z) is affine on each of the six Kuhn tetrahedra, but is not
    # globally affine; this distinguishes tetrahedral from trilinear sampling.
    u = np.max(points, axis=0)
    bounds = {"x0_m": 0.0, "x1_m": 1.0, "y0_m": 0.0, "y1_m": 1.0}
    assert_allclose(
        cell_mean(axes, ids, u, bounds, zval), 2.0 / 3.0 + zval**3 / 3.0, atol=2e-14
    )
    assert_allclose(
        sample(axes, ids, u, [0.32, 0.67, zval]), max(0.32, 0.67, zval), atol=2e-14
    )


@pytest.mark.parametrize("case_id", ["A", "B"])
def test_case_mesh_preserves_materials_power_gap_and_exterior_areas(case_id):
    data = json.loads((Path(__file__).parents[1] / "case.json").read_text())
    axes, layers, counts = mesh_axes(data, case_id, 0)
    fine_axes, _, fine_counts = mesh_axes(data, case_id, 1)
    for coarse, fine in zip(axes, fine_axes):
        assert_allclose(fine[::2], coarse, atol=2e-16)
    assert np.array_equal(fine_counts, 2 * counts)
    fields, q, flux = layer_fields(data, case_id, axes)
    dx, dy, dz = (np.diff(a) for a in axes)
    area = dx[:, None] * dy[None, :]
    occupied = np.array([fields[lid] > 0.0 for lid in layers])
    volume = np.sum(occupied * dz[:, None, None] * area[None, :, :])
    expected_volume = 0.0
    for layer, field in zip(data["stack"], fields):
        extent = layer["extent_m"]
        expected_area = (extent["x1"] - extent["x0"]) * (extent["y1"] - extent["y0"])
        assert_allclose(np.sum(area[field > 0.0]), expected_area, rtol=2e-14)
        expected_volume += expected_area * layer["thickness_m"]
    assert_allclose(volume, expected_volume, rtol=2e-14)
    assert_allclose(fields[1][fields[1] > 0.0], 70.76666666666667, rtol=2e-14)
    assert_allclose(fields[2][fields[2] > 0.0], 112.0, rtol=2e-14)
    power = data["cases"][case_id]["total_power_W"]
    assert_allclose(
        np.sum(q * area) * data["stack"][4]["thickness_m"], power, rtol=2e-14
    )
    assert_allclose(np.sum(flux * area), power, rtol=2e-14)
    if case_id == "A":
        xc, yc = ((a[:-1] + a[1:]) / 2 for a in axes[:2])
        gap = (
            (xc[:, None] > 0.001355)
            & (xc[:, None] < 0.0013601)
            & (yc[None, :] > 0.0152603)
            & (yc[None, :] < 0.0299403)
        )
        assert_allclose(np.sum(area[gap]), 7.4868e-8, rtol=2e-12)
        for layer in [3, 4]:
            assert np.all(fields[layer][gap] == 1.6)
    # Count each exposed Cartesian box face once, including underside annuli.
    exterior = np.sum(area * (occupied[0].astype(int) + occupied[-1]))
    exterior += np.sum(area[None, :, :] * (occupied[1:] != occupied[:-1]))
    exterior += np.sum(
        dz[:, None] * dy[None, :] * (occupied[:, 0, :].astype(int) + occupied[:, -1, :])
    )
    exterior += np.sum(
        dz[:, None, None]
        * dy[None, None, :]
        * (occupied[:, 1:, :] != occupied[:, :-1, :])
    )
    exterior += np.sum(
        dz[:, None] * dx[None, :] * (occupied[:, :, 0].astype(int) + occupied[:, :, -1])
    )
    exterior += np.sum(
        dz[:, None, None]
        * dx[None, :, None]
        * (occupied[:, :, 1:] != occupied[:, :, :-1])
    )
    convective = np.sum(area[occupied[-1]])
    assert_allclose(convective, 0.0324, rtol=2e-14)
    assert_allclose(exterior - convective, 0.0378288, rtol=2e-14)


def test_convergence_recovers_known_quadratic_limit_and_gci():
    # M(h) = 10 + h^2, with successive mesh widths 1, 1/2, 1/4.
    records = [
        {
            "M_C": 10.0 + 4.0 ** (-level),
            "M_die": "Chiplet_0",
            "refinement_factor": 2**level,
        }
        for level in range(3)
    ]
    result = convergence(records)
    assert result["status"] == "estimated"
    assert_allclose(result["observed_order"], 2.0, atol=1e-14)
    assert_allclose(result["M_extrapolated_C"], 10.0, atol=1e-14)
    assert_allclose(result["Richardson_correction_K"], -0.0625, atol=1e-14)
    assert_allclose(result["GCI_absolute_K"], 1.25 * 0.0625, atol=1e-14)
    assert result["accuracy_goal_met"]


@pytest.mark.parametrize(
    ("values", "winners"),
    [
        ([11.0, 10.0, 10.5], ["Chiplet_0"] * 3),
        ([10.0, 11.0, 12.0], ["Chiplet_0"] * 3),
        ([10.0, 11.0, 13.0], ["Chiplet_0"] * 3),
        ([10.0, 10.0, 10.0], ["Chiplet_0"] * 3),
        ([11.0, 10.25, 10.0625], ["Chiplet_0", "Chiplet_1", "Chiplet_1"]),
    ],
    ids=[
        "oscillating",
        "constant-differences",
        "growing-differences",
        "zero-differences",
        "changing-winner",
    ],
)
def test_convergence_rejects_sequences_without_a_valid_estimate(values, winners):
    records = [
        {"M_C": value, "M_die": winner, "refinement_factor": 2**level}
        for level, (value, winner) in enumerate(zip(values, winners))
    ]
    result = convergence(records)
    assert result["status"] == "inconclusive"
    assert "accuracy_goal_met" not in result
    assert "M_extrapolated_C" not in result
    assert "GCI_absolute_K" not in result


def test_verify_inputs_rejects_changed_reference_file(tmp_path):
    root = Path(__file__).parents[1]
    for name in ["case.json", "README.md", "CASE-SCHEMA.md", "LICENSE", "SHA256SUMS"]:
        shutil.copyfile(root / name, tmp_path / name)
    verify_inputs(tmp_path / "case.json")
    readme = tmp_path / "README.md"
    readme.write_bytes(readme.read_bytes() + b"\nChanged reference content.\n")
    with pytest.raises(ValueError, match="Reference manifest mismatch: README.md"):
        verify_inputs(tmp_path / "case.json")


@pytest.mark.parametrize("worker_case", [None, "A"])
def test_driver_rejects_unlisted_case_before_work(tmp_path, monkeypatch, worker_case):
    root = Path(__file__).parents[1]
    for name in ["case.json", "README.md", "CASE-SCHEMA.md", "LICENSE", "SHA256SUMS"]:
        shutil.copyfile(root / name, tmp_path / name)
    selected = tmp_path / "unlisted.json"
    shutil.copyfile(tmp_path / "case.json", selected)
    output = tmp_path / "results"
    argv = ["run_study.py", "--case-json", str(selected), "--output", str(output)]
    if worker_case:
        argv += ["--worker-case", worker_case]
    monkeypatch.setattr(sys, "argv", argv)

    def unexpected_work(*args, **kwargs):
        pytest.fail("Unverified selected case reached a worker")

    monkeypatch.setattr(run_study, "worker", unexpected_work)
    monkeypatch.setattr(run_study.subprocess, "run", unexpected_work)
    with pytest.raises(ValueError, match="Selected case is not the verified case.json"):
        run_study.main()
    assert not output.exists()


def test_verify_inputs_requires_case_manifest_entry(tmp_path):
    shutil.copyfile(Path(__file__).parents[1] / "case.json", tmp_path / "case.json")
    (tmp_path / "SHA256SUMS").write_text("")
    with pytest.raises(ValueError, match="case.json is missing from the reference manifest"):
        verify_inputs(tmp_path / "case.json")


def test_verify_inputs_rejects_replaced_case_with_new_manifest(tmp_path):
    case = tmp_path / "case.json"
    original = b'{"version": 1}'
    replacement = b'{"version": 2}'
    original_sha256 = hashlib.sha256(original).hexdigest()
    case.write_bytes(replacement)
    (tmp_path / "SHA256SUMS").write_text(
        hashlib.sha256(replacement).hexdigest() + "  case.json\n"
    )
    with pytest.raises(ValueError, match="changed after driver verification"):
        verify_inputs(case, original_sha256)


def test_verify_inputs_rejects_manifest_path_outside_input_directory(tmp_path):
    root = Path(__file__).parents[1]
    case = tmp_path / "case.json"
    shutil.copyfile(root / "case.json", case)
    (tmp_path / "SHA256SUMS").write_text(
        hashlib.sha256(case.read_bytes()).hexdigest() + "  " + str(root / "case.json") + "\n"
    )
    with pytest.raises(ValueError, match="Reference manifest path escapes input directory"):
        verify_inputs(case)
