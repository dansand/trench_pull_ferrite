# The loaded-edge boundary artifact — investigation record

> **Historical record, kept verbatim.** This was written in the private development archive (`ferrite_plate_flexure`,
> 2026-07) and is shipped because the SI quotes its numbers (the 29 % force- vs displacement-controlled split, the
> 12–14 % trench over-read, convergence within 40–60 km of the loaded edge). Names in it belong to the archive, not to
> this release: `render_np_diag.py`, `render_mid_decomp.py`, `render_paper_figures.py`, `render_gradient_compare.py`,
> the `dispload`/`bctest`/`set4` driver branches and the `out/`/`set1_rheology` paths do not exist here. In this release
> the corresponding figures are `scripts/render_hero.py` (Fig. 3, Branch B gradient) and `scripts/render_core_profiles.py`
> (Fig. S4), the models live in `data/suite1_strength/`, and "pseudo-density" is the equivalent density ρ̂ = τzx,x/g
> (`README.md`, `REPRODUCE.md`). The physics and numbers below are unchanged.

Record of a long forensic session on the loaded-trench-edge boundary layer and its effect on the pseudo-density
centroid / dipole arm / trench pull. Written before a context compaction. Baseline model throughout:
`set1_rheology/elastic_deep_60km_V4` (elastic) and `tresca_deep_150_60km_V4` (Tresca), massless, h=60 km.

## Headline (the thing not to lose)

Across **every** hero (elastic, Tresca, both DD-VMs) and **both** loading modes, the pseudo-density centroid /
dipole arm sits **near the mid-plane (~h/2 ≈ 30 km)**. That is the result — the trench-pull dipole compensates
at ~half the mechanical thickness, not at a deep mantle datum. Everything below is a **bounded, ~½-thickness
loaded-edge perturbation** on that ~30 km arm — it matters for rigour and for a referee, but it does **not** move
the headline.

## What the artifact is

- **Elastic, not plastic.** It appears in the pure-elastic run (σ_Y→∞), and the trench column is 0% yielded in
  every case (the edge yield-ramp keeps it elastic). So it is NOT plastic breakdown.
- **A loaded-edge / corner effect.** It's the σxx (bending-stress) boundary layer at the loaded free edge plus a
  corner stress-concentration where the applied traction meets the free surfaces. In clean bending ∂σxz/∂x is a
  positive parabola (zero at top/base, peak mid); the artifact is a departure from that near the corners.
- **Shape (Branch B / FE gradient):** a rotated-L — a strong **negative pseudo-density lobe hugging the top-left
  vertical boundary** (the dominant part, ~5 km scale) plus a fainter tongue along the top surface (~25 km). The
  bottom corner is near-clean. It grows with a reinforcing moment (moment traction is largest at the top corner).
- Node-by-node at the trench, top-10km pseudo-density: elastic −0.33, Tresca −0.44, V2/M+1.0 −0.83 kPa/m.

## Effect on the centroid — SIGNED ONLY (the |·| measure has been purged)

> **UPDATE 2026-07-09 — signed-centroid purge.** Two figures were secretly drawing a MAGNITUDE-weighted
> centroid ∫z|ρ̂|/∫|ρ̂| — `render_hero.py` panel (d) and `render_np_diag.py` `zarm`. That is non-physical: the
> *sign* of ∂σxz/∂x is the sign of the charge each point contributes to σzz below it, so abs-weighting throws
> away the dipole itself and, near the trench, gives the OPPOSITE trend to the real one. Both are now fixed to
> the **signed** centroid; the |·| measure is not used anywhere (repo-wide audit confirmed only those two).
> Panel (d) now correctly DEEPENS toward the trench. The paragraph below is kept for the record of the sign trap.

The **same negative top lobe** would move a signed vs a |·| centroid **oppositely**:
- **signed centroid** ∫zρ/∫ρ (this sets the dipole moment / arm / PULL, and is the ONLY one we plot) → **deeper**.
- **|ρ̂|-weighted centroid** ∫z|ρ|/∫|ρ| (PURGED — was the old panel-D line) → shallower (a spurious uptick).

Synthetic check (parabola centroid 30 + negative top lobe): signed → 34.7, |·| → 28.5. The old panel-D "uptick"
toward mid-plate near the trench was the |·| artifact; the signed pull is *inflated* (reads HIGH) by the same
lobe. Inflation is modest for V-only (arm ~34 vs interior ~30, ~13%) and large near the (V,M) breakdown limit.

## The displacement-based (bend-to-shape) test — AND AN OPEN CONCERN

Swapped the force load (parabolic shear traction) for a displacement load (prescribe trench-face u_z = the known
2899 m, `load_uz` in `solve_plastic`, driver branch `dispload`). Same deflection, both elastic:

| loading | w | ΔGPE* | arm | z_ps | top lobe |
|---|---|---|---|---|---|
| FORCE (parabola) | 2899 | **2.234** | 34.1 | 37.3 | −0.33 |
| DISPLACEMENT (flat block) | 2899 | **1.586** | 24.2 | 27.9 | **+3.55** |

- **The pull differs by 29% at the same deflection** — this is the OPEN CONCERN. The flat-block imposes its OWN
  artifact (a big *positive* top lobe, opposite sign), so it is a *worse* loading, not a clean one. The two
  **bracket** the interior value; neither edge read is right.
- **They converge inboard** (Saint-Venant): force−displacement pull difference is +0.65 at the trench, +0.05 at
  x=20, <0.02 (<1.5%) by x=40, <1% by 60 km. So the loaded-edge ambiguity is confined to ~½ a plate thickness.
- **Flat block does NOT change the curvature** (trench κ −1.81 vs −1.87 ×10⁻⁶, identical inboard) — only the
  through-thickness shear skin. So the bending shape is loading-robust; only the near-edge stress is not.
- Interior-arm-corrected pull ≈ g·Δρ·w·(interior arm ~30 km) ≈ **1.96 TN/m** (elastic), ~12% below the force read.

**Why the concern matters:** the headline pull uses the force (parabola) trench-column value, which over-reads
the interior by ~12–14%. The displacement test proves the trench read is loading-mode-dependent. We have NOT yet
explained the 29% split from first principles, only bracketed it.

## Gradient method: Branch A vs Branch B (RESOLVED)

Two gradients of ∂σxz/∂x lived in the repo: **A** = finite-difference `grad_x_field`; **B** = FE-exported
`dsxz_dx`. They diverged, especially at corners. Decisive check — compare both to the **gradient-free** arm
(pull/(gΔρw), from σzz only, cannot be cooked):

| x [km] | pull-based arm | branch A \|·\| | branch B \|·\| | branch B signed |
|---|---|---|---|---|
| 20 | 32.5 | 29.9 | 32.1 | 32.4 |
| 60 | 30.8 | 29.8 | 30.6 | 30.6 |
| 140 | 28.4 | 29.9 | 27.4 | 26.8 |

**Branch B tracks the gradient-free ground truth to ~1 km; Branch A is flat at h/2 and MISSES the real decline.**
So **A was the artifact (over-smoothing), B is correct.** The pull-based arm genuinely declines ~32→28 inboard —
it is NOT a flat h/2. The isostatic-column "jump" in B is a real local degeneracy (net pseudo-density→0 at the
sign reversal; any centroid there is ill-defined) — cosmetic, does not contaminate the x=20–140 decline.

**Audit done:** every CURRENT manuscript figure is on Branch B (`render_hero.py` via `pgrad`, `render_mid_decomp.py`).
Remaining Branch A: legacy `render_paper_figures.py`/`stress_depth_profiles` (old `out/production` models, not a
submitted figure) and deliberate comparison tools (`render_gradient_compare.py`, the `d_old` metric in
`render_core_profiles.py` whose *plotted* line is B). A brief isostatic mask was added to the hero and then
**reverted** at the user's instruction — hero centroid is now raw Branch B.

## Recourse / options on the table

1. **Move the nominal trench reference a little inboard** (past the ~5 km corner skin) — principled cutoff.
2. **Report the pull from the interior arm** (`ΔGPE*/gΔρw` at the interior plateau, ~30 km) rather than the
   boundary-layer trench column — Branch B confirms this is the honest value.
3. **Remove the artifact at source with a cleaner load** — kinematic/rotation BC (`u_x = θ(z−h/2)` = moment's
   displacement conjugate, Q3), or a down-dip loading buffer (most physical; needs geometry change).
4. **Reviewer framing:** ΔGPE* = ΔP_T·d with ΔP_T (topography) and d (interior neutral-plane depth) both
   loading-invariant; the identity ΔN_D=ΔGPE* holds to 0.02% for both loadings; the ~½-thickness edge skin
   brackets, not defines, the answer.

## Code touched this session

- `gpe_plastic.jl`: `solve_plastic` gained `load_uz` (displacement-controlled trench, ramped Dirichlet u_z).
- `paper_models.jl`: `run_tresca_deep` gained `load_uz`, `shear_parabolic`, `springs` kwargs; branches
  `dispload`, `bctest`, `set4`, `convergence`.
- `render_hero.py`: `pgrad` + `HERO_BRANCH` env (A/B, default B); mask added then reverted; **panel (d) centroid
  → SIGNED** `∫zρ̂/∫ρ̂` (was `np.abs`-weighted), clipped to trench limb `x<0.85·x_iso`.
- `render_np_diag.py`: **`zarm` → SIGNED centroid** (was "|τ|-weighted, hero style"), clipped to `x<0.7·x_iso`,
  M_max-only reference lines on panel A (V_max apparatus dropped so nothing points at the isostatic degeneracy).
- `render_mid_decomp.py`: pseudo-density → Branch B.
- Heroes re-rendered signed: paper (tresca, elastic, dd_vm_asym, dd_vm_sym, diff) + combined (s4_V2_Mp10, _Mp5).
