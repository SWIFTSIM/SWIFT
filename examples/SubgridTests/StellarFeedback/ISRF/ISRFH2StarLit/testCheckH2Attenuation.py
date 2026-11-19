#!/usr/bin/env python3
"""Run checkH2Attenuation.py on synthetic snapshots and print the exit code of each case.

Usage: testCheckH2Attenuation.py WORK_DIR

WORK_DIR is a scratch directory. The script runs the three makeIC.py files in it (random
placement, level 4), builds synthetic snapshots on their gas, and runs the check on a good
case and on every planted error. A planted error must give the exit class named in the table.
The exit code of this script is 0 when every row matches.
"""

import contextlib
import copy
import io
import json
import os
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import checkH2Attenuation as gate  # noqa: E402
import h2_fixture_common as gate_common  # noqa: E402

ROOT = os.path.dirname(HERE)
GOOD_PARAMS = """GEARFeedback:
  ISRF_H2_line_absorption: 1
  with_photoionization: 0
  HII_couple_ionization_rate: 0
  with_stellar_wind_feedback: 0
  supernovae_efficiency: 0
  ISRF_packet_light_speed: unlimited
"""
START_UP_LOG = ""
GATES = ("starlit", "twosource", "velocity")
ROWS = []
U_L = 3.0856775814913673e21


def run_gate(argv):
    if argv[0] in GATES and "--startup-log" not in argv:
        argv = [*argv, "--startup-log", START_UP_LOG]
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = gate.main(argv)
    return code, out.getvalue()


def record(name, expect, argv, key):
    code, text = run_gate(argv)
    lines = [
        ln for ln in text.splitlines() if ln.startswith(("FAIL", "PASS", "INFO N2"))
    ]
    detail = next((ln for ln in lines if key in ln), lines[0] if lines else "")
    ROWS.append((name, expect, code, detail[:110]))


def make_ic(work, name, level=4):
    d = os.path.join(work, name)
    os.makedirs(d, exist_ok=True)
    script = os.path.join(ROOT, f"ISRFH2{name}", "makeIC.py")
    subprocess.run(
        [sys.executable, script, "--level", str(level)],
        cwd=d,
        check=True,
        capture_output=True,
    )
    ic = [f for f in os.listdir(d) if f.startswith("ICs")][0]
    with h5py.File(os.path.join(d, ic), "r") as f:
        pos = f["PartType0/Coordinates"][:]
        rho = f["PartType0/Densities"][:].astype(float)
        h = f["PartType0/SmoothingLength"][:].astype(float)
        box = float(f["Header"].attrs["BoxSize"][0])
    return d, pos, rho, h, box


def write_snap(
    path,
    pos,
    rho,
    h,
    box,
    temp,
    x_h2,
    lw,
    ratio,
    ids=None,
    drop=(),
    time=8.0,
    column=None,
):
    """Write a synthetic snapshot; fields in drop are left out."""
    n = len(pos)
    ids = np.arange(n) if ids is None else ids
    species = {k: np.full(n, 1e-20) for k in gate.SPECIES_MASS_NUMBER}
    species.update(HI=np.full(n, 0.76), HeI=np.full(n, 0.24), HII=np.full(n, 7.6e-9))
    species["H2I"] = x_h2
    n_per_mass = sum(species[k] / a for k, a in gate.SPECIES_MASS_NUMBER.items())
    n_per_mass = n_per_mass + sum(
        c * species[k] for k, c in gate.ELECTRONS_PER_MASS.items()
    )
    energy = (
        1.5
        * gate.K_B
        * temp
        * n_per_mass
        / gate.M_H
        / (3.0856775814913673e21 / 3.0856775814913673e16) ** 2
    )
    fields = {
        **species,
        "ParticleIDs": ids,
        "Coordinates": pos,
        "SmoothingLengths": h,
        "Densities": rho,
        "InternalEnergies": energy,
        "H2I": x_h2,
        "LWSpecificEnergies": lw,
        "LWDissociatingSpecificEnergies": lw * ratio,
        "LWLineTransmissions": ratio,
        "LWOwnSpecificEnergies": lw,
        "LWOwnDissociatingSpecificEnergies": lw * ratio,
        "LWLineH2ColumnDensities": np.full(n, 1e17) if column is None else column,
    }
    with h5py.File(path, "w") as f:
        hdr = f.create_group("Header")
        hdr.attrs["BoxSize"] = [box, box, box]
        hdr.attrs["Time"] = [time]
        u = f.create_group("Units")
        u.attrs["Unit length in cgs (U_L)"] = 3.0856775814913673e21
        u.attrs["Unit mass in cgs (U_M)"] = 1.988409870698051e43
        u.attrs["Unit time in cgs (U_t)"] = 3.0856775814913673e16
        g = f.create_group("PartType0")
        for k, v in fields.items():
            if k not in drop:
                g.create_dataset(k, data=v)


def tool_args(cmd, d, last, first, params, extra=()):
    return [
        cmd,
        "--first",
        first,
        "--last",
        last,
        "--fixture",
        os.path.join(d, "fixture.json"),
        "--params",
        params,
        *extra,
    ]


def starlit_cases(work):
    d, pos, rho, h, box = make_ic(work, "StarLit")
    fix = json.load(open(os.path.join(d, "fixture.json")))
    snap = {
        "Units": {
            "unit_length": 3.0856775814913673e21,
            "unit_mass": 1.988409870698051e43,
        }
    }
    star = np.array(fix["star_kpc"])
    dd = pos - star
    dd -= box * np.round(dd / box)
    r = np.linalg.norm(dd, axis=1) * 3.0856775814913673e21
    rho_cgs = rho * 1.988409870698051e43 / 3.0856775814913673e21**3
    temp = np.full(len(pos), 50.0)
    rng = np.random.default_rng(7)
    x_h2 = np.full(len(pos), 2.0 * fix["h2_ratio"] * 0.76) * (
        1.0 + 0.03 * rng.uniform(-1.0, 1.0, len(pos))
    )
    n_h2 = rho_cgs * x_h2 / (2.0 * gate.M_H)
    dens = rho_cgs / gate.M_H
    lw = 1.0 / np.maximum(r, 1e18) ** 2 * 1e36
    good = gate.transmission(n_h2 * r, temp, dens)
    params = os.path.join(d, "used_parameters.yml")
    open(params, "w").write(GOOD_PARAMS)
    counters = os.path.join(d, "counters.json")
    json.dump(
        {"mean_hop_path_cm": float(np.mean(h) * 3.0856775814913673e21)},
        open(counters, "w"),
    )
    first = os.path.join(d, "s0.hdf5")
    write_snap(first, pos, rho, h, box, temp, x_h2, lw, good, time=1.0)

    def snap_with(name, ratio, **kw):
        p = os.path.join(d, name)
        kw["column"] = kw.get("column", n_h2 * r) * U_L**2  # stored in U_L^-2
        write_snap(
            p,
            pos,
            rho,
            h,
            box,
            kw.pop("temp", temp),
            kw.pop("x_h2", x_h2),
            kw.pop("lw", lw),
            ratio,
            **kw,
        )
        return p

    def args(last, p=params, extra=("--counters", counters)):
        return tool_args("starlit", d, last, first, p, extra)

    record(
        "S1 good run", 0, args(snap_with("good.hdf5", good)), "PASS  own transmission"
    )
    h_cm = h * 3.0856775814913673e21
    mask = (r > 0.0) & (r < 0.5 * box * 3.0856775814913673e21 - 2.0 * h_cm)
    nudged = good * np.exp(0.4 * gate.column_spread_bar(r, n_h2, mask, temp, dens))
    record(
        "S2 inside the n_H2 spread bar (0.4 of it)",
        0,
        args(snap_with("nudge.hdf5", np.minimum(nudged, 1.0))),
        "transmission",
    )
    record(
        "C1 carried column field doubled",
        1,
        args(snap_with("c1.hdf5", good, column=2.0 * n_h2 * r)),
        "FAIL  carried column",
    )
    record(
        "C2 carried column field halved",
        1,
        args(snap_with("c2.hdf5", good, column=0.5 * n_h2 * r)),
        "FAIL  carried column",
    )
    n_lo, n_hi = float(n_h2[mask].min()), float(n_h2[mask].max())
    record(
        "C3 carried column field 0.5 per cent off, clipped to the n_H2 range",
        0,
        args(snap_with("c3.hdf5", good, column=np.clip(1.005 * n_h2, n_lo, n_hi) * r)),
        "PASS  carried column",
    )
    bad_col = n_h2 * r
    bad_col[5] = np.nan
    record(
        "C4 carried column field non-finite",
        3,
        args(snap_with("c4.hdf5", good, column=bad_col)),
        "non-finite",
    )
    record(
        "N1 switch off: line fields absent",
        2,
        args(
            snap_with(
                "n1.hdf5",
                good,
                drop=(
                    "LWLineTransmissions",
                    "LWDissociatingSpecificEnergies",
                    "LWOwnSpecificEnergies",
                    "LWOwnDissociatingSpecificEnergies",
                ),
            )
        ),
        "lacks",
    )
    n2 = gate.transmission(2.0 * n_h2 * r, temp, dens)
    record(
        "N2 doubled n_H2 in the hop",
        1,
        args(snap_with("n2.hdf5", n2)),
        "FAIL  own transmission",
    )
    bad = good.copy()
    bad[5] = np.nan
    record("NaN in the ratio", 3, args(snap_with("nan.hdf5", bad)), "non-finite")
    record(
        "L3 t = 1 everywhere on lit probes (no line absorption)",
        1,
        args(snap_with("void.hdf5", np.ones(len(pos)))),
        "FAIL  count of lit",
    )
    # Light front: every probe at r < (hops - 1) h_min is lit; the snapshots are 8 hops apart.
    in_dom = np.where(mask)[0]
    r_full = 7.0 * float(h_cm[mask].min())
    inside = in_dom[r[in_dom] < r_full]
    r_full3 = 2.0 * float(h_cm[mask].min())
    reach2 = 2.0 * gate.KERNEL_GAMMA * float(h_cm[mask].max())
    ROWS.append(
        (
            "L0 the domain holds probes beyond r_full (3 hops) and beyond the reach (2 hops)",
            0,
            0 if (r[in_dom] >= r_full3).any() and (r[in_dom] > reach2).any() else 1,
            f"r_full {r_full3:.3g}, reach {reach2:.3g}, max r {r[in_dom].max():.3g} cm",
        )
    )
    few = np.where(r < 0.5 * r_full, lw, 0.0)
    record(
        "L1 lit count too small (only r < r_full / 2 lit)",
        1,
        args(snap_with("few.hdf5", good, lw=few)),
        "FAIL  lit count",
    )
    hole = lw.copy()
    hole[inside[len(inside) // 2]] = 0.0
    record(
        "L2 one probe inside the front unlit",
        1,
        args(snap_with("hole.hdf5", good, lw=hole)),
        "FAIL  lit set complete",
    )
    record(
        "L4 3 hops, probes beyond the front unlit: not read as t = 1",
        0,
        args(
            snap_with("beyond.hdf5", good, lw=np.where(r < r_full3, lw, 0.0), time=3.0)
        ),
        "PASS  lit set complete",
    )
    record(
        "L5 2 hops, a probe lit beyond the reach",
        1,
        args(snap_with("reach.hdf5", good, time=2.0)),
        "FAIL  lit count",
    )
    record(
        "void run: no LW light",
        4,
        args(snap_with("nolw.hdf5", good, lw=np.zeros(len(pos)))),
        "no gas particle",
    )
    p_hii = os.path.join(d, "hii.yml")
    open(p_hii, "w").write(
        GOOD_PARAMS.replace("with_photoionization: 0", "with_photoionization: 1")
    )
    record(
        "void run: HII feedback on",
        4,
        args(snap_with("good.hdf5", good), p_hii),
        "premise",
    )
    x_wide = x_h2 * 10.0 ** rng.uniform(-1.0, 1.0, len(pos))
    n_wide = rho_cgs * x_wide / (2.0 * gate.M_H)
    good_wide = gate.transmission(n_wide * r, temp, dens)
    first_wide = os.path.join(d, "s0_wide.hdf5")
    write_snap(first_wide, pos, rho, h, box, temp, x_wide, lw, good_wide, time=1.0)
    record(
        "N2 not detectable (n_H2 spread of 100, bar far larger)",
        4,
        tool_args(
            "starlit",
            d,
            snap_with("wide.hdf5", good_wide, x_h2=x_wide),
            first_wide,
            params,
            ("--counters", counters),
        ),
        "doubled",
    )
    moved = snap_with("moved.hdf5", good, x_h2=x_h2 * 1.01)
    record("H2 not frozen (1 per cent change)", 4, args(moved), "H2 changed")
    p_t = os.path.join(d, "t50.yml")
    open(p_t, "w").write(GOOD_PARAMS + "SPH:\n  initial_temperature: 50\n")
    record(
        "S3 temperature scale from SPH:initial_temperature",
        0,
        args(snap_with("good.hdf5", good), p_t),
        "PASS  own transmission",
    )
    edge = pos.copy()
    edge[1948 % len(pos), 0] = 0.0
    face = os.path.join(d, "face0.hdf5")
    write_snap(face, edge, rho, h, box, temp, x_h2, lw, good, time=1.0)
    record(
        "E1 a gas particle exactly on a box face",
        4,
        [
            "starlit",
            "--first",
            face,
            "--last",
            os.path.join(d, "good.hdf5"),
            "--fixture",
            os.path.join(d, "fixture.json"),
            "--params",
            params,
            "--counters",
            counters,
        ],
        "box face",
    )
    try:
        gate_common.check_away_from_faces(edge, box, "gas")
        raised = 0
    except RuntimeError:
        raised = 1
    ROWS.append(
        (
            "E2 makeIC refuses a coordinate on a face",
            1,
            raised,
            "RuntimeError from check_away_from_faces",
        )
    )
    startup = _write(
        d, "startup.log", "h2_transported_lw_rate_shielding = 0\n", text=True
    )
    record(
        "start-up log says shielded_sum",
        4,
        [
            "starlit",
            "--first",
            first,
            "--last",
            os.path.join(d, "good.hdf5"),
            "--fixture",
            os.path.join(d, "fixture.json"),
            "--params",
            params,
            "--startup-log",
            startup,
            "--counters",
            counters,
        ],
        "start-up log",
    )
    record(
        "start-up log missing",
        2,
        [
            "starlit",
            "--first",
            first,
            "--last",
            os.path.join(d, "good.hdf5"),
            "--fixture",
            os.path.join(d, "fixture.json"),
            "--params",
            params,
            "--startup-log",
            os.path.join(d, "none.log"),
            "--counters",
            counters,
        ],
        "start-up log",
    )
    return d, (pos, rho, h, box, lw, good, temp, x_h2)


def _write(d, name, content, text=False):
    p = os.path.join(d, name)
    with open(p, "w") as f:
        f.write(content if text else json.dumps(content))
    return p


def stale_key_cases(work):
    d = os.path.join(work, "stale_key")
    os.makedirs(d, exist_ok=True)
    line = (
        "[0000.1] cooling.c:cooling_init():799: h2_transported_lw_rate_shielding = {}\n"
    )
    ignored = _write(d, "ignored.log", line.format(1), text=True)
    honoured = _write(d, "honoured.log", line.format(0), text=True)
    silent = _write(d, "silent.log", "main: done. Bye.\n", text=True)
    record(
        "N3 stale key ignored: starts, library reports 1",
        0,
        ["stale_key", "--log", ignored, "--exit-status", "0"],
        "library reports",
    )
    record(
        "N3 stale key honoured: library reports 0",
        1,
        ["stale_key", "--log", honoured, "--exit-status", "0"],
        "library reports",
    )
    record(
        "N3 run stopped at start-up",
        1,
        ["stale_key", "--log", ignored, "--exit-status", "1"],
        "cleanly",
    )
    record(
        "N3 no library line in the log",
        1,
        ["stale_key", "--log", silent, "--exit-status", "0"],
        "library reports",
    )
    record(
        "N3 log missing",
        2,
        ["stale_key", "--log", os.path.join(d, "none.log"), "--exit-status", "0"],
        "missing",
    )


def parity_cases(work, d, data):
    pos, rho, h, box, lw, good, temp, x_h2 = data
    k = 12.0
    n = len(pos)
    rng = np.random.default_rng(3)
    pa = os.path.join(d, "p_a.hdf5")
    write_snap(pa, pos, rho, h, box, temp, x_h2, lw, good)
    dis = np.exp(rng.uniform(-0.5, 0.5, n) * k * gate.FLOAT32_EPS)
    pb = os.path.join(d, "p_b.hdf5")
    write_snap(pb, pos, rho, h, box, temp, x_h2, lw * dis, good)
    perm = rng.permutation(n)
    pc = os.path.join(d, "p_perm.hdf5")
    write_snap(
        pc,
        pos[perm],
        rho[perm],
        h[perm],
        box,
        temp[perm],
        x_h2[perm],
        lw[perm],
        good[perm],
        ids=perm,
    )
    i = 7
    drop = lw.copy()
    drop[i] *= 1.0 - 1.0 / k
    pd = os.path.join(d, "p_drop.hdf5")
    write_snap(pd, pos, rho, h, box, temp, x_h2, drop, good)
    # Previous snapshot: the faintest tenth is dark (light front), particle i is lit.
    front = lw <= np.quantile(lw, 0.1)
    front[i] = False
    pp = os.path.join(d, "p_prev.hdf5")
    write_snap(
        pp, pos, rho, h, box, temp, x_h2, np.where(front, 0.0, lw), good, time=7.0
    )
    # A front particle 1e-6 of its neighbours, off by 1e-3 of its own value.
    j = int(np.flatnonzero(front)[0])
    fa = lw.copy()
    fa[j] *= 1e-6
    pfa = os.path.join(d, "p_front_a.hdf5")
    write_snap(pfa, pos, rho, h, box, temp, x_h2, fa, good)
    fb = fa.copy()
    fb[j] *= 1.0 + 1e-3
    pf = os.path.join(d, "p_front.hdf5")
    write_snap(pf, pos, rho, h, box, temp, x_h2, fb, good)
    # The same 1e-3 on a front particle as bright as its neighbours.
    fl = lw.copy()
    fl[j] *= 1.0 + 1e-3
    pfl = os.path.join(d, "p_front_lit.hdf5")
    write_snap(pfl, pos, rho, h, box, temp, x_h2, fl, good)
    pn = os.path.join(d, "p_nan.hdf5")
    bad = lw.copy()
    bad[3] = np.nan
    write_snap(pn, pos, rho, h, box, temp, x_h2, bad, good)
    pm = os.path.join(d, "p_nocol.hdf5")
    write_snap(
        pm, pos, rho, h, box, temp, x_h2, lw, good, drop=("LWLineH2ColumnDensities",)
    )
    pt = os.path.join(d, "p_time.hdf5")
    write_snap(pt, pos, rho, h, box, temp, x_h2, lw, good, time=9.0)
    cnt = _write(
        d,
        "pc.json",
        {
            "probes_on_other_rank": 100,
            "foreign_bookings": 500,
            "steps": 9,
            "held_lw_slots_at_restart": 30,
            "refits_after_restart": 20,
        },
    )
    zero = _write(
        d,
        "pz.json",
        {"probes_on_other_rank": 100, "foreign_bookings": 0, "steps": 9},
    )
    nofit = _write(
        d, "pf.json", {"held_lw_slots_at_restart": 30, "refits_after_restart": 0}
    )

    def restart(b, c=cnt, stop=pp):
        return [
            "parity",
            "--a",
            pa,
            "--b",
            b,
            "--mode",
            "restart",
            "--counters",
            c,
            "--stop",
            stop,
        ]

    def ranks(b, c=cnt, a=pa):
        return [
            "parity",
            "--a",
            a,
            "--b",
            b,
            "--mode",
            "ranks",
            "--counters",
            c,
            "--previous",
            pp,
        ]

    p0 = os.path.join(d, "p_first.hdf5")
    write_snap(p0, pos, rho, h, box, temp, x_h2, lw, good, time=1.0)

    def reduced_ranks(b, first=p0, a=pa):
        return [
            "parity",
            "--a",
            a,
            "--b",
            b,
            "--mode",
            "ranks",
            "--previous",
            pp,
            *(("--first", first) if first else ()),
        ]

    record(
        "R1 restart, no counters: reduced",
        5,
        ["parity", "--a", pa, "--b", pa, "--mode", "restart", "--stop", pp],
        "bitwise",
    )
    record("R2 2 ranks, no counters: reduced", 5, reduced_ranks(pb), "round-off")
    record("R3 2 ranks, no counters, dropped booking", 1, reduced_ranks(pd), "FAIL")
    record(
        "R4 2 ranks, no counters, no --first", 2, reduced_ranks(pb, first=None), "needs"
    )
    record("P1 restart: identical", 0, restart(pa), "bitwise")
    record("P2 restart: reordered rows, same values", 0, restart(pc), "bitwise")
    record("P3 restart: one ulp-level change", 1, restart(pb), "bitwise")
    record("P4 2 ranks: round-off within the bar", 0, ranks(pb), "round-off")
    record("D1 drop one booking: restart", 1, restart(pd), "FAIL")
    record("D2 drop one booking: 2 ranks, lit before", 1, ranks(pd), "FAIL")
    record(
        "F1 2 ranks: front particle off by 1e-3 of its own value",
        0,
        ranks(pf, a=pfa),
        "front",
    )
    record("F2 2 ranks: bright front particle off by 1e-3", 1, ranks(pfl), "front")
    record("P5 2 ranks: no foreign booking (void)", 4, ranks(pa, zero), "foreign")
    record(
        "P6 restart: nothing lit after the stop (void)", 4, restart(pa, stop=pa), "lit"
    )
    record(
        "P7 restart: no refit after the stop (void)", 4, restart(pa, nofit), "refits"
    )
    record("P8 2 ranks: column field missing", 2, ranks(pm), "lacks")
    record("P9 2 ranks: non-finite value", 3, ranks(pn), "non-finite")
    record("P10 restart: snapshots at different times (void)", 4, restart(pt), "times")
    # The ratio alone cannot see the dropped booking: both energies drop together.
    ratio_a = good
    ratio_d = (lw * (1 - 1 / k) * good)[i] / drop[i]
    ROWS.append(
        (
            "D3 drop one booking: ratio-only compare (defect A)",
            "pass",
            0 if abs(ratio_d - ratio_a[i]) < 1e-15 else 1,
            f"ratio {ratio_d:.17g} against {ratio_a[i]:.17g}",
        )
    )


def twosource_cases(work):
    d, pos, rho, h, box = make_ic(work, "TwoSources", level=5)
    fix = json.load(open(os.path.join(d, "fixture.json")))
    cm = 3.0856775814913673e21
    rho_cgs = rho * 1.988409870698051e43 / cm**3
    x_h2 = np.full(len(pos), 2.0 * 3e-4 * 0.76)
    dens = rho_cgs / gate.M_H
    lw = np.full(len(pos), 1.0)
    t_iso = np.full(len(pos), fix["temperature_K"])
    ii = np.ones(len(pos))
    held_off = np.ones(len(pos))
    for p in fix["probes"]:
        k = p["id"]
        rules = gate.merge_rules(p["shares"], t_iso[k], dens[k], p["downstream_N_cm2"])
        ii[k] = rules["ii"]
        held_off[k] = float(
            gate.transmission(
                p["shares"][1]["N_cm2"] + p["downstream_N_cm2"], t_iso[k], dens[k]
            )
        )
    params = os.path.join(d, "used_parameters.yml")
    open(params, "w").write(GOOD_PARAMS)
    first = os.path.join(d, "s0.hdf5")
    write_snap(first, pos, rho, h, box, t_iso, x_h2, lw, ii)
    first_cold = os.path.join(d, "s0_cold.hdf5")
    write_snap(first_cold, pos, rho, h, box, np.full(len(pos), 50.0), x_h2, lw, ii)

    def go(name, ratio, expect, counters, key, label, first_snap=first, temp=t_iso):
        last = os.path.join(d, name)
        write_snap(last, pos, rho, h, box, temp, x_h2, lw, ratio)
        record(
            label,
            expect,
            tool_args(
                "twosource",
                d,
                last,
                first_snap,
                params,
                ("--counters", counters) if counters else (),
            ),
            key,
        )

    go(
        "good.hdf5",
        ii,
        0,
        None,
        "PASS  ratio equals",
        "T1 two sources: form (ii) at 163.4 K",
    )
    go("n4.hdf5", held_off, 1, None, "FAIL  ratio equals", "N4 merge held term removed")
    cold = np.full(len(pos), 50.0)
    go(
        "cold.hdf5",
        ii,
        4,
        None,
        "fit not uniform",
        "T2 density contrast at 50 K (fit not uniform)",
        first_snap=first_cold,
        temp=cold,
    )


def falsifier_cases(work):
    d = os.path.join(work, "falsifier")
    os.makedirs(d, exist_ok=True)
    cm = 3.0856775814913673e21
    n_part = 4
    pos = np.full([n_part, 3], 0.1) + np.arange(n_part)[:, None] * 0.01
    rho_int = 100.0 * gate.M_H * cm**3 / 1.988409870698051e43
    rho = np.full(n_part, rho_int)
    h = np.full(n_part, 1e-6)
    temp = np.full(n_part, 100.0)
    x_h2 = np.full(n_part, 1e-3)
    lw = np.ones(n_part)
    params = os.path.join(d, "used_parameters.yml")
    open(params, "w").write(GOOD_PARAMS)
    cnt = _write(d, "c.json", {"mean_hop_path_cm": 1.0})
    first = os.path.join(d, "s0.hdf5")
    write_snap(first, pos, rho, h, 1.0, temp, x_h2, lw, np.ones(n_part))
    for down, tag in ((1e18, "D = 1e18"), (1e14, "D = 1e14")):
        shares = [{"w": 0.5, "N_cm2": 0.0}, {"w": 0.5, "N_cm2": 1e20}]
        fixture = {
            "probes": [
                {"id": 0, "r_b_kpc": 1.0, "downstream_N_cm2": down, "shares": shares}
            ]
        }
        fx = _write(d, "fix.json", fixture)
        rules = gate.merge_rules(shares, 100.0, 100.0, down)
        for rule, expect, label in (
            ("ii", 0, "form (ii)"),
            ("i", 1, "form (i)"),
            ("r3", 1, "r3 rule"),
        ):
            ratio = np.ones(n_part)
            ratio[0] = rules[rule]
            last = os.path.join(d, f"{rule}_{int(down)}.hdf5")
            write_snap(last, pos, rho, h, 1.0, temp, x_h2, lw, ratio)
            want = expect if down == 1e18 else 4
            key = "ratio equals" if down == 1e18 else "discriminate"
            record(
                f"F {label} at {tag}: R = {rules[rule] / rules['true']:.4f}",
                want,
                [
                    "twosource",
                    "--first",
                    first,
                    "--last",
                    last,
                    "--fixture",
                    fx,
                    "--params",
                    params,
                    "--counters",
                    None,
                ],
                key,
            )


def velocity_cases(work):
    d, pos, rho, h, box = make_ic(work, "VelocityShift", level=5)
    fix = json.load(open(os.path.join(d, "fixture.json")))
    cm = 3.0856775814913673e21
    rho_cgs = rho * 1.988409870698051e43 / cm**3
    temp = np.full(len(pos), 50.0)
    x_h2 = np.full(len(pos), 2.0 * 3e-4 * 0.76)
    dens = rho_cgs / gate.M_H
    n_h2 = rho_cgs * x_h2 / (2.0 * gate.M_H)
    star = np.array(fix["star_kpc"])
    dd = pos - star
    dd -= box * np.round(dd / box)
    r = np.linalg.norm(dd, axis=1) * cm
    h_cm = h * cm
    face = fix["cloud_radius_kpc"][0] * cm
    depth = np.maximum(r - face, 0.0)
    lw = np.ones(len(pos))

    def ratio_of(column):
        return np.array(
            [
                float(gate.transmission(column[k], temp[k], dens[k]))
                for k in range(len(pos))
            ]
        )

    # A correct code carries the in-cloud path plus part of the resetting hop: 0.74 h is
    # the median excess measured on the real run.
    base = n_h2 * (depth + 0.74 * h_cm)
    own = ratio_of(base)
    full = ratio_of(n_h2 * r)
    params = os.path.join(d, "used_parameters.yml")
    open(params, "w").write(GOOD_PARAMS)
    first = os.path.join(d, "s0.hdf5")
    write_snap(first, pos, rho, h, box, temp, x_h2, lw, own)
    slow = dict(fix, cloud_speed_widths=1.0)
    json.dump(slow, open(os.path.join(d, "fixture_slow.json"), "w"))

    def go(name, ratio, expect, fixture, label):
        last = os.path.join(d, name)
        write_snap(last, pos, rho, h, box, temp, x_h2, lw, ratio)
        a = tool_args("velocity", d, last, first, params)
        a[a.index("--fixture") + 1] = os.path.join(d, fixture)
        record(label, expect, a, "cloud probes")

    go("v_reset.hdf5", own, 0, "fixture.json", "V1 ten widths: own column (reset)")
    go(
        "v_noreset.hdf5",
        full,
        1,
        "fixture.json",
        "V2 ten widths: column carried (no reset)",
    )
    go(
        "v_double.hdf5",
        ratio_of(2.0 * base),
        1,
        "fixture.json",
        "V3 ten widths: doubled column",
    )
    go(
        "v_late.hdf5",
        ratio_of(np.maximum(base - n_h2 * h_cm, 0.0)),
        1,
        "fixture.json",
        "V4 ten widths: reset one hop late",
    )
    go(
        "v_x15.hdf5",
        ratio_of(1.5 * base),
        1,
        "fixture.json",
        "V5 ten widths: column 1.5 times too large",
    )
    nan = own.copy()
    nan[int(np.argmax(r * (r > face) * (r < 0.3 * box * cm)))] = np.nan
    go("v_nan.hdf5", nan, 3, "fixture.json", "V6 ten widths: NaN in the ratio")
    go(
        "v_slow.hdf5",
        full,
        0,
        "fixture_slow.json",
        "V7 one width: full column (control)",
    )
    go(
        "v_slow_bad.hdf5",
        own,
        1,
        "fixture_slow.json",
        "V8 one width: reset seen (control broken)",
    )


def pi_l_case():
    """The bound is 1 / pi_l: a share whose history ran in other gas makes pi_l
    smaller than its energy share w_l, and correct code then exceeds 1 / w_l."""
    shares = [{"w": 0.5, "N_cm2": 0.0}, {"w": 0.5, "N_cm2": 1e18, "t": 0.01}]
    rules = gate.merge_rules(shares, 100.0, 100.0, 1e16)
    R = rules["ii"] / rules["true"]
    ok = 1.0 / rules["w_l"] < R <= 1.0 / rules["pi_l"]
    ROWS.append(
        (
            f"F carried column, other gas before: R = {R:.3f} above 1/w_l = "
            f"{1.0 / rules['w_l']:.1f}, below 1/pi_l = {1.0 / rules['pi_l']:.1f}",
            0,
            0 if ok else 1,
            f"pi_l {rules['pi_l']:.4g}, w_l {rules['w_l']:.4g}",
        )
    )


def main(work):
    global START_UP_LOG
    os.makedirs(work, exist_ok=True)
    START_UP_LOG = _write(
        work, "start_up.log", "h2_transported_lw_rate_shielding = 1\n", text=True
    )
    anchor_ok = (
        abs(float(gate.shielding(0, 100, 100)) - 0.99997) < 2e-5
        and abs(float(gate.log_slope(1e18, 100, 100)) - 0.52) < 0.005
        and abs(
            float(gate.shielding(2e18, 100, 100) / gate.shielding(1e18, 100, 100))
            - 0.696
        )
        < 0.002
    )
    ROWS.append(
        (
            "A0 f_sh anchors: f(0)=0.99997, slope 0.52, f(2N)/f(N)=0.696 (design 5.2)",
            0,
            0 if anchor_ok else 1,
            f"slope {float(gate.log_slope(1e18, 100, 100)):.4f}",
        )
    )
    d, data = starlit_cases(work)
    stale_key_cases(work)
    parity_cases(work, d, data)
    twosource_cases(work)
    falsifier_cases(work)
    pi_l_case()
    velocity_cases(work)
    width = max(len(r[0]) for r in ROWS)
    print(f"{'case':<{width}}  want  got  verdict  line the check printed")
    bad = 0
    for name, want, got, detail in ROWS:
        ok = (want == "pass" and got == 0) or want == got
        bad += not ok
        print(
            f"{name:<{width}}  {str(want):>4}  {got:>3}  {'ok' if ok else 'WRONG':<7}  {detail}"
        )
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
