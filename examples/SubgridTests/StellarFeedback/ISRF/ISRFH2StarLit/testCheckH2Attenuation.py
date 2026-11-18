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

ROOT = os.path.dirname(HERE)
GOOD_PARAMS = """GEARFeedback:
  ISRF_H2_line_absorption: 1
  with_photoionization: 0
  HII_couple_ionization_rate: 0
  with_stellar_wind_feedback: 0
  supernovae_efficiency: 0
  ISRF_packet_light_speed: unlimited
GeajayCooling:
  h2_transported_lw_rate_shielding: background_only
"""
ROWS = []


def run_gate(argv):
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


def write_snap(path, pos, rho, h, box, temp, x_h2, lw, ratio, ids=None, drop=()):
    """Write a synthetic snapshot; fields in drop are left out."""
    n = len(pos)
    ids = np.arange(n) if ids is None else ids
    fields = {
        "ParticleIDs": ids,
        "Coordinates": pos,
        "SmoothingLengths": h,
        "Densities": rho,
        "Temperatures": temp,
        "H2I": x_h2,
        "LWSpecificEnergies": lw,
        "LWDissociatingSpecificEnergies": lw * ratio,
        "LWLineTransmissions": ratio,
    }
    with h5py.File(path, "w") as f:
        f.create_group("Header").attrs["BoxSize"] = [box, box, box]
        u = f.create_group("Units")
        u.attrs["Unit length in cgs (U_L)"] = 3.0856775814913673e21
        u.attrs["Unit mass in cgs (U_M)"] = 1.988409870698051e43
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
    x_h2 = np.full(len(pos), 2.0 * fix["h2_ratio"] * 0.76)
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
    write_snap(first, pos, rho, h, box, temp, x_h2, lw, good)

    def snap_with(name, ratio, **kw):
        p = os.path.join(d, name)
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

    record("S1 good run", 0, args(snap_with("good.hdf5", good)), "PASS  transmission")
    pos_fix = np.abs(np.log(good))
    nudged = good * np.exp(
        0.4
        * gate.log_slope(np.maximum(n_h2 * r, 1.0), temp, dens)
        * np.mean(h)
        * 3.0856775814913673e21
        / np.maximum(r, 1.0)
    )
    record(
        "S2 inside the hop bar (0.4 of it)",
        0,
        args(snap_with("nudge.hdf5", np.minimum(nudged, 1.0))),
        "transmission",
    )
    record(
        "N1 switch off: line fields absent",
        2,
        args(
            snap_with(
                "n1.hdf5",
                good,
                drop=("LWLineTransmissions", "LWDissociatingSpecificEnergies"),
            )
        ),
        "lacks",
    )
    n2 = gate.transmission(2.0 * n_h2 * r, temp, dens)
    record(
        "N2 doubled n_H2 in the hop",
        1,
        args(snap_with("n2.hdf5", n2)),
        "FAIL  transmission",
    )
    bad = good.copy()
    bad[5] = np.nan
    record("NaN in the ratio", 3, args(snap_with("nan.hdf5", bad)), "non-finite")
    record(
        "void run: no line absorption (count 0)",
        4,
        args(snap_with("void.hdf5", np.ones(len(pos)))),
        "void",
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
    record(
        "N2 not detectable (bar far larger)",
        4,
        args(
            snap_with("good.hdf5", good),
            extra=("--counters", _write(d, "wide.json", {"mean_hop_path_cm": 1e30})),
        ),
        "doubled",
    )
    moved = snap_with("moved.hdf5", good, x_h2=x_h2 * 1.01)
    record("H2 not frozen (1 per cent change)", 4, args(moved), "H2 changed")
    startup = _write(
        d, "startup.log", "h2_transported_lw_rate_shielding = 0\n", text=True
    )
    record(
        "option c: log says shielded_sum",
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
            "--mode",
            "c",
            "--startup-log",
            startup,
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


def refusal_cases(work):
    d = os.path.join(work, "refusal")
    os.makedirs(d, exist_ok=True)
    msg = (
        "[0000.1] cooling_io.h:cooling_read():95: GEARFeedback:ISRF_H2_line_absorption: 1 is not "
        "supported: the geajay library shields the transported LW rate\n"
    )
    ok = _write(d, "refused.log", msg, text=True)
    plain = _write(d, "ran.log", "main: done. Bye.\n", text=True)
    word = _write(
        d, "word.log", "note: ISRF_H2_line_absorption: 1 is not supported\n", text=True
    )
    record(
        "N3 refused: exit 1 and the text",
        0,
        ["refusal", "--log", ok, "--exit-status", "1"],
        "refusal text",
    )
    record(
        "N3 not refused: exit 0, no text",
        1,
        ["refusal", "--log", plain, "--exit-status", "0"],
        "FAIL",
    )
    record(
        "N3 text only in a comment line (not error())",
        1,
        ["refusal", "--log", word, "--exit-status", "1"],
        "refusal text",
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
    drop = lw.copy()
    drop_ratio = good.copy()
    i = 7
    drop[i] *= 1.0 - 1.0 / k
    pd = os.path.join(d, "p_drop.hdf5")
    write_snap(pd, pos, rho, h, box, temp, x_h2, drop, good)
    cnt = _write(
        d,
        "pc.json",
        {
            "probes_on_other_rank": 100,
            "foreign_bookings": 500,
            "max_terms_per_particle": k,
            "held_lw_slots_at_restart": 30,
        },
    )
    zero = _write(
        d,
        "pz.json",
        {
            "probes_on_other_rank": 100,
            "foreign_bookings": 0,
            "max_terms_per_particle": k,
        },
    )
    record(
        "P1 restart: identical",
        0,
        ["parity", "--a", pa, "--b", pa, "--mode", "restart", "--counters", cnt],
        "bitwise",
    )
    record(
        "P2 restart: reordered rows, same values",
        0,
        ["parity", "--a", pa, "--b", pc, "--mode", "restart", "--counters", cnt],
        "bitwise",
    )
    record(
        "P3 restart: one ulp-level change",
        1,
        ["parity", "--a", pa, "--b", pb, "--mode", "restart", "--counters", cnt],
        "bitwise",
    )
    record(
        "P4 2 ranks: round-off within k u",
        0,
        ["parity", "--a", pa, "--b", pb, "--mode", "ranks", "--counters", cnt],
        "round-off",
    )
    record(
        "D1 drop one booking: restart",
        1,
        ["parity", "--a", pa, "--b", pd, "--mode", "restart", "--counters", cnt],
        "FAIL",
    )
    record(
        "D2 drop one booking: 2 ranks",
        1,
        ["parity", "--a", pa, "--b", pd, "--mode", "ranks", "--counters", cnt],
        "FAIL",
    )
    record(
        "P5 2 ranks: no foreign booking (void)",
        4,
        ["parity", "--a", pa, "--b", pa, "--mode", "ranks", "--counters", zero],
        "foreign",
    )
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
    cnt = _write(d, "held.json", {"held_lw_slots_at_merge": 40})
    cnt0 = _write(d, "held0.json", {"held_lw_slots_at_merge": 0})

    def go(name, ratio, expect, counters, key, label, first_snap=first, temp=t_iso):
        last = os.path.join(d, name)
        write_snap(last, pos, rho, h, box, temp, x_h2, lw, ratio)
        record(
            label,
            expect,
            tool_args(
                "twosource", d, last, first_snap, params, ("--counters", counters)
            ),
            key,
        )

    go(
        "good.hdf5",
        ii,
        0,
        cnt,
        "PASS  ratio equals",
        "T1 two sources: form (ii) at 163.4 K",
    )
    go("n4.hdf5", held_off, 1, cnt, "FAIL  ratio equals", "N4 merge held term removed")
    go(
        "void.hdf5",
        held_off,
        4,
        cnt0,
        "held_lw",
        "N4 void: no slot holds light at the merge",
    )
    cold = np.full(len(pos), 50.0)
    go(
        "cold.hdf5",
        ii,
        4,
        cnt,
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
    cnt = _write(d, "c.json", {"held_lw_slots_at_merge": 5, "mean_hop_path_cm": 1.0})
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
            key = {"ii": "main", "i": "sign", "r3": "bound"}[rule]
            key = (
                {"ii": "ratio equals", "i": "sign", "r3": "bound"}[rule]
                if down == 1e18
                else "discriminate"
            )
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
                    cnt,
                ],
                key,
            )


def velocity_cases(work):
    d, pos, rho, h, box = make_ic(work, "VelocityShift")
    fix = json.load(open(os.path.join(d, "fixture.json")))
    cm = 3.0856775814913673e21
    rho_cgs = rho * 1.988409870698051e43 / cm**3
    temp = np.full(len(pos), 50.0)
    x_h2 = np.full(len(pos), 2.0 * 3e-4 * 0.76)
    dens = rho_cgs / gate.M_H
    star = np.array(fix["star_kpc"])
    dd = pos - star
    dd -= box * np.round(dd / box)
    r = np.linalg.norm(dd, axis=1) * cm
    x_rel = dd[:, 0] * cm
    slab = [(v - star[0]) * cm for v in fix["slab_x_kpc"]]
    cloud = [(v - star[0]) * cm for v in fix["cloud_x_kpc"]]
    n_h2 = fix["n_h2_cm3"]
    lw = np.ones(len(pos))
    own = np.ones(len(pos))
    full = np.ones(len(pos))
    for k in range(len(pos)):
        pc_ = gate.path_through(cloud[0], cloud[1], 0.0, x_rel[k], r[k])
        ps_ = gate.path_through(slab[0], slab[1], 0.0, x_rel[k], r[k])
        own[k] = float(gate.transmission(n_h2 * pc_, temp[k], dens[k]))
        full[k] = float(gate.transmission(n_h2 * (pc_ + ps_), temp[k], dens[k]))
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
        "v_slow.hdf5",
        full,
        0,
        "fixture_slow.json",
        "V3 one width: full column (control)",
    )
    go(
        "v_slow_bad.hdf5",
        own,
        1,
        "fixture_slow.json",
        "V4 one width: reset seen (control broken)",
    )


def main(work):
    os.makedirs(work, exist_ok=True)
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
    refusal_cases(work)
    parity_cases(work, d, data)
    twosource_cases(work)
    falsifier_cases(work)
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
