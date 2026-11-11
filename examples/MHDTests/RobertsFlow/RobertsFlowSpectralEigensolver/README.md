# Roberts-flow dynamo eigensolver

Reference solutions for the kinematic Roberts-flow dynamos I–IV in a periodic box, for comparison with
the SWIFT and Pencil runs (Shchutskyi, Brandenburg & Schaller, *Convergence of Roberts flow dynamos in
Eulerian and Lagrangian codes*).

| File | Purpose |
|---|---|
| `solver.py` | Velocity fields (`FLOWS`); growth rate λ, frequency ω, wavenumber kz of the fastest-growing mode, and the Table-6 quantities, for many η in parallel |
| `eigen_patterns.py` | Magnetic-field patterns r₁–r₃, θ₁–θ₃ of the eigensolution (as in Figures 6–9), and comparison with a code slice (eqs. 2.30–2.32) |

Requirements: Python ≥ 3.8, `numpy`, `scipy`, `pandas`, `matplotlib`; `h5py` for writing and reading slice files.

## Method

The flows (eq. 2.1, k₀ = 1, Lx = Ly = 2π)

    u = (v0 sin x cos y, −v0 cos x sin y, w0 vz(x, y)),
    vz = sin x sin y (I), cos x cos y (II), (cos 2x + cos 2y)/2 (III), sin x (IV)

do not depend on z, so every growing solution has the form B = Re[ b(x, y) e^{i kz z + p t} ] with
p = λ + iω.  Following Roberts (1972, *Phil. Trans. R. Soc. A* 271, 411, §§ 8 and 10), only bx and by are
solved for (bz follows from ∇·B = 0, so the divergence constraint holds exactly), expanded in Fourier
modes e^{i(mx + ny)} with |m|, |n| ≤ s.  This turns the induction equation into a sparse matrix eigenvalue
problem for each kz.  In a box of height Lz only kz = 2πj/Lz is allowed.

For each η, `solver.growth()`:

1. scans kz = 0.5, 1, 1.5, … (for Lz = 4π) with a coarse eigensolve, until the growth rate has clearly
   passed its maximum in kz (the optimal kz increases as η decreases);
2. refines the most promising (kz, branch) candidates, increasing s until the eigenvalue changes by less
   than `tol` between successive s **and** the resolution diagnostic z(s) < √tol;
3. checks at the final s that no other branch grows faster;
4. computes the Table-6 quantities from the eigenvector.

Small problems use dense LAPACK; large ones use ARPACK shift-invert with shifts placed at the frequencies
where growing modes occur.

## `solver.py`

### Run the parameter table

```bash
python solver.py
```

reads `./test_run_parameters.csv` (tab- or comma-separated, columns `Run #`, `eta`, `Flow_kind` (1–4 or
I–IV), `Status`; only rows with `Status == done` are used), runs all of them in parallel and writes
`results.csv` (tab-separated).  A line is printed as each run finishes:

```
[3/40] elapsed 0:12:31  run 17 flow II eta=0.01: lambda=+0.00496 omega=+0.48351 kz=0.5 s=43 ok (55 s)
```

Rows are appended to `results.csv` as they finish, so an interrupted run keeps its completed results.

### From Python

```python
from solver import growth, run_jobs, read_parameter_csv
import numpy as np

r = growth('I', eta=0.05, w0=np.sqrt(2), tol=1e-12, k1=1.0)   # one point
print(r['p'].real, r['kz'], r['Bmean_rms/Brms'])

if __name__ == "__main__":                                     # required on macOS for parallel runs
    jobs = [(i, 'I', eta) for i, eta in enumerate(np.logspace(-3, np.log10(0.6), 40))]
    rows = run_jobs(jobs, out='flowI.csv', tol=1e-4, workers=8)
```

### Main parameters of `growth()` (also accepted by `run_jobs`, `scan`)

| Parameter | Default | Meaning |
|---|---|---|
| `tol` | 1e-4 | absolute accuracy of p = λ + iω (units v₀k₀); also sets the resolution requirement z(s) < √tol |
| `w0` | 1.0 | vertical flow amplitude (w₀ = v₀ in Figure 2; w₀ = √2 v₀ gives v_rms = 1 for flows I–III) |
| `Lz` | 4π | box height (2π for a cubic box) |
| `kzmax` | 10 | hard limit of the kz scan |
| `c_coarse`, `s_coarse` | 2.5, – | resolution of the kz scan, s = c_coarse/√η unless `s_coarse` is given |
| `kz_margin`, `n_kz`, `n_branch` | 0.02, 4, 2 | which coarse candidates are refined |
| `smax` | 200 | largest truncation |
| `k1` | kz | wavenumber normalising J_rms |

Recommendation: `tol=1e-4` or smaller wherever λ is small (near onset, flow II's dynamo window);
`tol=1e-3` is a 20 % error when λ ≈ 0.005.

### Output columns (`results.csv`)

| Column | Meaning |
|---|---|
| `growth_rate`, `omega` | λ and ω of the fastest mode (the sign of ω only gives the travel direction) |
| `omega_mod` | 2|ω|, the B_rms modulation frequency when a ±ω pair interferes (§ 2.4) |
| `kz` | vertical wavenumber of the fastest mode |
| `kz_scanned` | largest kz examined; should be well above `kz` |
| `s`, `z` | final truncation and edge diagnostic; eigenvalue error ≈ z² |
| `converged` | False if `smax` was reached before `tol` was met |
| `Rm` | v_rms/(η k₀) |
| Table-6 columns | `Bmean_rms/Brms`, `Jrms/(k1 Brms)`, `<J.B>/(Jrms Brms)`, `<(J.B)^2>^1/2/(Jrms Brms)`, `<(JxB)^2>^1/2/(Jrms Brms)`, `-W_L/(urms Jrms Brms)` with W_L = ⟨u·(J×B)⟩, `<J^2 B^2>/(Jrms^2 Brms^2)`, `-<A.B>/Brms^2` (Coulomb gauge); volume averages, μ₀ = 1 |

### Run time and parallelisation

Each η runs in its own process (one BLAS thread each), smallest η first.  Indicative single-core times:
η = 0.05: seconds; η = 0.01: about 1 min; η = 2.6×10⁻³: about 20 min; η = 10⁻³: about an hour (estimated).

## `eigen_patterns.py`

```bash
python eigen_patterns.py -flow I  -eta 0.05 -kz 1.0 -res_xy 64 -outim eig_I.png  -out eig_I.h5
python eigen_patterns.py -flow IV -eta 0.05 -kz 0.5 -res_xy 64 -outim eig_IV.png -out eig_IV.h5
python eigen_patterns.py -flow I  -eta 0.05 -kz 1.0 -res_xy 64 -compare slice.h5
```

| Option | Default | Meaning |
|---|---|---|
| `-flow`, `-eta`, `-kz`, `-w0` | I, 0.05, 1.0, 1.0 | flow, diffusivity, vertical wavenumber of the mode, w₀ |
| `-res_xy` | 64 | pixels per direction; use the resolution of the code slice |
| `-x0` | 0 | coordinate of the first pixel centre (change it if the code grid starts elsewhere, e.g. −π) |
| `-pol` | x | which mode to show when two modes have the same growth rate (see below) |
| `-tol` | 1e-8 | accuracy of the eigenvalue |
| `-outim`, `-out` | | figure and slice file (`-out ""` skips the file) |
| `-compare` | | slice `.h5` from `FT_slice_pixel_show.py`; writes `<outim>_diff.png` and prints ⟨|δB̃|⟩_rms (eq. 2.32) |

**Convention**, identical to `FT_slice_pixel_show.py`: B is normalised by its box rms (eq. 2.6);
b(x, y) is the coefficient of e^{+i kz z} (numpy FFT divided by N, i.e. b̃/(2 B_rms)); one global phase
makes Im⟨b_x⟩ = 0 and Re⟨b_x⟩ > 0 (eq. 2.8); pixels are at xᵢ = i·2π/N; the `.h5` file has the same
datasets (`r1`–`r3`, `theta1`–`theta3`, `ReB1`–`ReB3`, `ImB1`–`ImB3`).  The map is evaluated exactly (a
direct Fourier sum), so there is no aliasing at low pixel resolution.

**Modes with equal growth rates.**  The fastest growth rate can belong to two modes:
flows II and III have a ±ω pair of travelling waves, and flow IV (e.g. η = 0.05, kz = 0.5) has a
degenerate real eigenvalue.  The script then returns the mode that the paper's initial condition (2.12),
B ∝ e_x sin(kz z), excites: the one with non-zero ⟨B_x⟩ and, for a degenerate eigenvalue, ⟨B_y⟩ = 0.
`-pol y` gives the other one.  The choice is printed as a note.

**When comparing with `FT_slice_pixel_show.py`:** its FFT index `modez` must correspond to the same kz,
`modez = kz·Lz/2π` (kz = 1 in Lz = 4π → `modez = 2`).

## Other flows and the original Roberts (1972) setup

### Changing the velocity field

All flows are defined in one place, the `FLOWS` dictionary near the top of `solver.py`; every other part
of the code (operator, v_rms and Rm, Lorentz work, patterns) uses it.  To study another flow, add one
entry:

```python
FLOWS['myflow'] = lambda x, y, v0, w0: (v0*np.sin(y), v0*np.sin(x), w0*np.cos(x)*np.cos(y))
```

and pass `'myflow'` as the flow name: `growth('myflow', ...)`, `converge('myflow', ...)`, or
`python eigen_patterns.py -flow myflow ...`.  Requirements:

* u depends on x and y only and is 2π-periodic in both (z is the invariant direction);
* each component is a trigonometric polynomial of degree ≤ 4 (the Fourier coefficients are computed
  exactly with a 16-point FFT; raise `n` and `amax` in `_coeffs` for higher harmonics);
* ∇·u = 0 is assumed by the induction operator, as for all Roberts flows.

`v0` and `w0` scale the horizontal and vertical parts (defaults 1).  Nothing else needs to change.

### Original Roberts (1972) setup

Roberts' motions are invariant in **x**, with the field ∝ exp(pt + ijx).  Renaming his x → z (his y, z →
our x, y) gives four predefined flows; his j is our `kz` and his λ is our `eta`:

| Flow name | Roberts 1972 | u in this code |
|---|---|---|
| `R72_I` | motion I (5.1): (cos y − cos z, sin z, sin y) | (sin y, sin x, cos x − cos y) |
| `R72_II` | motion II (6.1): (cos y + cos z, sin z, sin y) | (sin y, sin x, cos x + cos y) |
| `R72_III` | motion III (6.2): (2 cos y cos z, sin z, sin y) | (sin y, sin x, 2 cos x cos y) |
| `R72_IV` | motion IV (6.3): (sin(y + z), sin 2z, sin 2y) | (sin 2y, sin 2x, sin(x + y)) |

Roberts' domain is infinite in the invariant direction, so any j is allowed.  For a given j, use
`converge` (fixed kz, no box scan) instead of `growth`:

```python
from solver import converge
r = converge('R72_I', eta=1/8, kz=0.5, tol=1e-9)
print(r['p'])            # growth rate (Roberts' Re p) + i * frequency, converged in s
```

To find the fastest j, scan `kz` yourself, or use `growth('R72_I', eta, Lz=...)` with a box tall enough
that the allowed kz = 2πj/Lz resolve the maximum (e.g. `Lz=40*np.pi` gives kz steps of 0.05).

**Comparison with Roberts' Table 4** (same λ, j and truncation s as in his table):

| Motion | λ, j | s | Roberts | This code | |
|---|---|---|---|---|---|
| I | 1/8, 1/2 | 3, 5, 7 | 0.166, 0.169165, 0.1691710 | 0.16635, 0.169165, 0.1691709 | ✓ |
| I | 1/32, 1 and 1/64, 1 | 9, 13 | 0.14613, 0.13297 | 0.14613, 0.13297 | ✓ |
| I | 1, 1/2 | 3–13 | +0.0840389 | −0.0840379 | sign differs (see below) |
| II | 1/8, 1/4 | 3, 5 | 0.023347 ± 0.4676i, 0.023345 ± 0.463825i | 0.023347 + 0.46760i, 0.023345 + 0.463825i | ✓ |
| II | 1/16, 1/4 | 7 | 0.01593 ± 0.4768i | 0.015930 − 0.47681i | ✓ |
| III | 1/32, 0.294 | 6, 8 | 0.0638 − 0.5517i, 0.06415 − 0.5521i | 0.06380 − 0.55171i, 0.06415 − 0.55198i | ✓ |
| III | 1/32, 1 | 8 | 0.0497 − 1.934i | 0.04973 − 1.93354i | ✓ |
| III | 1/2, 1/2 | 4–7 | −0.0983255646 − 0.426276464i | −0.1183313912 − 0.426276464i | Im ✓, Re differs by 0.0200 |
| IV | 1/4, 1 | 6, 7, 8 | −0.1556805, −0.15567850, −0.155678597 | −0.155686, −0.155677, −0.155679 | ✓ (5–6 digits) |

The sign of the frequency only gives the travel direction (Roberts lists ±).  Three rows do not agree,
and the reason is not known: motion I at λ = 1 has the same magnitude but the opposite sign, even though
λ = 1/8, 1/32 and 1/64 agree to all digits; motion III at λ = 1/2 has the frequency right to 9 digits but
a growth rate 0.0200 lower; and at the smallest resistivities (motion III, λ = 1/128; motion IV,
λ = 1/64 and 1/256) Roberts' own resolution diagnostic is large (z = 0.025–0.47), and the converged
values here differ from his (e.g. motion IV, λ = 1/64, j = 1/16: 0.000168 here, converged at s = 39,
against 0.0001987 at his s = 25).  These are worth checking against the printed paper before relying on
those rows.

**Relation to the paper's flows.**  `R72_I` is the paper's flow I with w₀ = √2 v₀, rotated by 45° in
the xy plane and with lengths and velocities scaled by √2.  Growth rates map as
λ_I(η, kz) = λ_R72(2η, kz/√2); for example `converge('I', 1/16, 0.5*np.sqrt(2), w0=np.sqrt(2))` and
`converge('R72_I', 1/8, 0.5)` both give 0.16917096627.  The same rotation relates the paper's flows II
and III to `R72_II` and `R72_III`, and every eigenvalue of `R72_II` / `R72_III` also appears in the
spectrum of the paper's flow.  However, the paper's 2π × 2π box contains two of Roberts' cells and
therefore allows additional field patterns, and for flows II and III these grow fastest (e.g. at
η = 1/16, w₀ = √2: λ = −0.0107 for the paper's flow II against −0.0760 for `R72_II` at η = 1/8).
The paper's flow IV (vz = sin x) is a different flow from Roberts' motion IV.

## Validation

| Check | Reference | Eigensolver |
|---|---|---|
| Roberts 1972 Table 4, λ = 1/8, j = 1/2, s = 7 | 0.1691710 | 0.1691709 |
| η_crit, w₀ = v₀, Lz = 4π (I, II, III, IV) | Pencil 64²×128: 0.495, 0.145, 0.255, 0.220 | 0.501, 0.146, 0.254, 0.220 |
| Flow IV, η = 0.011: λ, ω_mod | Table 3: 0.0336, 0.1397 | 0.0337, 0.1398 |
| Flow I, η = 0.05, w₀ = √2 v₀: λ and Table-6 quantities | Pencil 256²×512 | agree to 10–12 digits |
| Patterns, η = 0.05 (flows I, II, IV) | Figures 6, 8, 9 | agree visually |

## Known limitations

* The Table-6 quantities in `solver.py` are computed for one eigenvector.  Where two modes share the
  fastest growth rate (flow IV's degenerate eigenvalue, the ±ω pairs of flows II and III), that is an
  arbitrary member of the pair; λ and ω are unaffected, but ⟨B̄⟩-related quantities may not match a
  simulation seeded with a specific polarisation.  `eigen_patterns.py` makes the selection explicitly.
* The search for the fastest mode can still miss a branch that is only resolved at high s.  If a result
  looks suspicious, increase `c_coarse`, `n_kz` or `n_branch` and check whether it changes.
* `growth()` only considers the box modes kz = 2πj/Lz; for an infinite domain use `converge()` at
  chosen kz, or a tall box (see above).

## Credits

Method: G. O. Roberts (1972), *Phil. Trans. R. Soc. Lond. A* 271, 411–454.
Initial version developed with the assistance of Claude (Opus 5.5, Anthropic), November 2026;
reviewed and validated by M. Schaller and N. Shchutskyi.
