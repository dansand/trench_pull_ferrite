"""tests/test_quick.py — fast checks (seconds).  Run from the repository root:   pytest tests/

No solve and no figure design is exercised here: synthetic data where a known answer exists, the shipped reference
model where a real field is needed, and the two SI benchmark scripts (which assert their own misfits).  The Julia
constitutive kernel is checked by tests/test_quick.jl, run from here when `julia` is on the PATH.
"""
import glob, hashlib, os, shutil, subprocess, sys
import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
PY = sys.executable
REF = "data/suite1_strength/tresca_deep_150_60km_V4"


@pytest.fixture(autouse=True)
def _run_from_root(monkeypatch):
    monkeypatch.chdir(ROOT)


def _sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _synthetic_model(tmp_path, x, z, ux, uz, S):
    """Write a minimal gpe_model.vtu (rectangular node grid, displacement u, constant 2nd-PK stress S) and load it."""
    import pyvista as pv
    from gpe_analysis import Model
    X, Z = np.meshgrid(x, z, indexing="ij")
    pts = np.column_stack([X.ravel(), Z.ravel(), np.zeros(X.size)])
    g = pv.PolyData(pts).cast_to_unstructured_grid()
    g.point_data["u"] = np.column_stack([ux(X, Z).ravel(), uz(X, Z).ravel(), np.zeros(X.size)])
    for name, val in (("sigma_xx [Pa]", S[0, 0]), ("sigma_zz [Pa]", S[1, 1]), ("sigma_xz [Pa]", S[0, 1])):
        g.point_data[name] = np.full(X.size, val)
    d = tmp_path / "synthetic"; d.mkdir()
    g.save(str(d / "gpe_model.vtu"))
    return Model(str(d))


def test_cauchy_pushforward_rigid_rotation(tmp_path):
    """A rigid rotation R of the reference grid (u = RX − X, J = 1) must push S forward to R S Rᵀ exactly."""
    th = 0.05; c, s = np.cos(th), np.sin(th)
    R = np.array([[c, -s], [s, c]]); S = np.array([[3.0e6, 0.7e6], [0.7e6, -1.2e6]])
    x = np.linspace(0, 4000.0, 41); z = np.linspace(0, 100.0, 11)
    m = _synthetic_model(tmp_path, x, z, lambda X, Z: (c - 1) * X - s * Z, lambda X, Z: s * X + (c - 1) * Z, S)
    Cxx, Czz, Cxz = m.cauchy_fields()
    C = R @ S @ R.T
    assert np.allclose(Cxx, C[0, 0], rtol=1e-9, atol=1.0)
    assert np.allclose(Czz, C[1, 1], rtol=1e-9, atol=1.0)
    assert np.allclose(Cxz, C[0, 1], rtol=1e-9, atol=1.0)
    assert not np.allclose(Cxz, S[0, 1], rtol=1e-3)          # the rotation genuinely mixes components (O(θ) on the shear)
    assert m.stress_frame() == "massless"                     # no stamp, no manifest ⇒ the σzz heuristic (small stress)


def test_signed_centroid_sign_changing_profile():
    """Signed centroid of τ(z) = z/L − 0.4 on [0, L]: ∫zτ/∫τ = (L²/3 − 0.2L²)/(L/2 − 0.4L) = 4L/3 — outside the
    support, unlike any magnitude-weighted centre, which is the point of keeping it signed."""
    from gpe_analysis import signed_centroid
    L = 60e3; z = np.linspace(0, L, 6001); tau = z / L - 0.4
    zc = signed_centroid(z, tau)
    assert abs(zc / (4 * L / 3) - 1) < 1e-4
    zabs = np.trapz(z * np.abs(tau), z) / np.trapz(np.abs(tau), z)
    assert 0 < zabs < L and abs(zabs - zc) > 0.3 * L


def test_scalar_vs_batched_deformed_resultants():
    """The vectorised many-line extractor must reproduce the scalar one on interior lines of the reference model."""
    from gpe_analysis import Model
    m = Model(REF)
    xs = [40.0, 100.0, 250.0]
    sz, nd = m.batch_deformed_resultants(xs)
    for i, xk in enumerate(xs):
        s1, n1 = m.deformed_resultants(xk)
        assert abs(sz[i] / s1 - 1) < 1e-9 and abs(nd[i] / n1 - 1) < 1e-9


def test_reference_model_reads_manifest_and_headline_numbers():
    """The shipped reference model resolves its stress frame from data/DATA_MANIFEST.json (no heuristic), and the
    paper's headline numbers hold: ΔGPE* = 2.542 TN/m (±0.5 %), identity < 0.05 %, arm 35.4 km (±0.2; w_T on the trench column)."""
    from gpe_analysis import Model, trench_pull, trench_deflection
    m = Model(REF)
    assert m.manifest_entry().get("stress_frame") == "massless" and m.stress_frame() == "massless"
    dG, dN, x_I = trench_pull(m)
    assert abs(dG / 2.542e12 - 1) < 0.005, dG
    assert abs(dG - dN) / abs(dG) < 5e-4
    arm = dG / ((3300.0 - 1000.0) * 9.81 * trench_deflection(m))
    assert abs(arm / 1e3 - 35.4) < 0.2, arm


def test_write_read_table_round_trip(tmp_path, monkeypatch):
    """tables/: write_table writes a PLAIN csv (header + rows) plus a provenance sidecar json; read_table returns the
    sidecar's meta and typed rows; floats are 6-s.f. formatted so a rerun is byte-stable; sequence rows and dict rows
    write identically."""
    import json
    from gpe_analysis import write_table, read_table
    monkeypatch.chdir(tmp_path)
    rows_seq = [("a", 1.23456789, 3), ("b", -2.0e-7, 4)]
    rows_dict = [{"name": "a", "value": 1.23456789, "n": 3}, {"name": "b", "value": -2.0e-7, "n": 4}]
    p1 = write_table("t1", ["name", "value", "n"], rows_seq, script="s.py", figure="figures/f.png", models=["data/m"], meta={"slope": 0.6476543})
    p2 = write_table("t2", ["name", "value", "n"], rows_dict, script="s.py", figure="figures/f.png", models=["data/m"], meta={"slope": 0.6476543})
    assert open(p1).read() == open(p2).read()
    assert open(p1).read().splitlines()[0] == "name,value,n"                      # nothing before the header
    prov = json.load(open("tables/t1.json"))
    assert prov == {"written_by": "s.py", "figure": "figures/f.png", "models": ["data/m"], "meta": {"slope": 0.647654}}
    meta, rows = read_table(p1)
    assert meta == {"slope": 0.647654}
    assert rows == [{"name": "a", "value": 1.23457, "n": 3.0}, {"name": "b", "value": -2e-07, "n": 4.0}]
    assert open(p1).read() == open(write_table("t1", ["name", "value", "n"], rows_seq, script="s.py", figure="figures/f.png", models=["data/m"], meta={"slope": 0.6476543})).read()


def test_every_committed_table_has_provenance_sidecar():
    """Every tables/*.csv is a plain CSV (header first) with a provenance sidecar naming the script that wrote it."""
    import json
    from gpe_analysis import read_table
    files = sorted(glob.glob("tables/*.csv"))
    assert len(files) >= 14, files
    for p in files:
        assert not open(p).readline().startswith("#"), p
        side = os.path.splitext(p)[0] + ".json"
        assert os.path.isfile(side), side
        assert json.load(open(side))["written_by"], side
        meta, rows = read_table(p)
        assert rows and all(isinstance(r, dict) for r in rows), p


def test_manifest_check_passes():
    r = subprocess.run([PY, "analysis/make_manifest.py", "--check"], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_convergence_script_refuses_without_data():
    """Without data/convergence/ the script must exit nonzero and leave the recorded tables untouched."""
    if os.path.isdir("data/convergence"):
        pytest.skip("data/convergence/ present — the refusal path is not exercised")
    before = _sha("tables/convergence.md") + _sha("tables/convergence.csv")
    r = subprocess.run([PY, "scripts/render_convergence.py"], capture_output=True, text=True)
    assert r.returncode != 0
    assert _sha("tables/convergence.md") + _sha("tables/convergence.csv") == before


@pytest.mark.parametrize("script", sorted(glob.glob("scripts/render_*.py")))
def test_import_has_no_side_effects(script):
    """Importing a figure module must not render, write a figure or table, or change the matplotlib backend."""
    mod = os.path.basename(script)[:-3]
    snap = {p: os.path.getmtime(p) for p in glob.glob("figures/*") + glob.glob("tables/*") + glob.glob("data/*.md")}
    code = ("import sys, matplotlib; sys.path[:0] = ['analysis', 'scripts']; b0 = matplotlib.get_backend()\n"
            f"import {mod}\n"
            "assert matplotlib.get_backend() == b0, (b0, matplotlib.get_backend())\n"
            "import matplotlib.pyplot as plt; assert plt.get_fignums() == [], plt.get_fignums()\n")
    r = subprocess.run([PY, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    assert {p: os.path.getmtime(p) for p in glob.glob("figures/*") + glob.glob("tables/*") + glob.glob("data/*.md")} == snap


@pytest.mark.parametrize("script", ["scripts/render_benchmark.py", "scripts/render_mp_benchmark.py"])
def test_benchmark_misfits_assert(script):
    """The SI benchmarks assert their own misfits in-script (S1: V ≤ 1 %, end resultant 4.000 ± 1 %, parabola ≤ 0.5 %;
    S2: 148 sections, mean ≤ 0.15 %, max ≤ 1 %) — the guard against the 2026-09-11 stress-frame regression."""
    r = subprocess.run([PY, script], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-1500:] + r.stderr[-1500:]
    assert "benchmark assertions passed" in r.stdout


def test_julia_constitutive_kernel():
    """Plastic return lands on the yield surface; the consistent tangent matches a finite difference (tests/test_quick.jl)."""
    if shutil.which("julia") is None:
        pytest.skip("julia not on PATH")
    r = subprocess.run(["julia", "--project=.", "tests/test_quick.jl"], capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
