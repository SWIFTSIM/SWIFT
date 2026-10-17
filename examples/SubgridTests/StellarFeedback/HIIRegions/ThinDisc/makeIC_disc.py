"""
Initial conditions for a thin, small gas disc hosting young massive stars.

The gas follows an exponential radial profile and a Gaussian vertical
profile. Single stars of 10 to 40 Msun are placed in a few compact clusters
and as field stars drawn from the gas distribution. The disc is at rest and
starts neutral and cold.
"""

import argparse

import h5py
import numpy as np

M_P_G = 1.67262192e-24
MSUN_G = 1.98841e33
KPC_CM = 3.0856775814913673e21
K_B = 1.380649e-16


def parse_options():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--gas_mass", type=float, default=20.0, help="Gas particle mass [Msun]"
    )
    p.add_argument("--disc_mass", type=float, default=4e6, help="Total gas mass [Msun]")
    p.add_argument(
        "--scale_length_pc", type=float, default=150.0, help="Radial scale length [pc]"
    )
    p.add_argument(
        "--scale_height_pc",
        type=float,
        default=20.0,
        help="Gaussian vertical sigma [pc]",
    )
    p.add_argument(
        "--truncation", type=float, default=5.0, help="Disc radius in scale lengths"
    )
    p.add_argument(
        "--temperature", type=float, default=100.0, help="Initial gas temperature [K]"
    )
    p.add_argument("--n_clusters", type=int, default=1, help="Number of star clusters")
    p.add_argument(
        "--stars_per_cluster", type=int, default=62, help="Stars per cluster"
    )
    p.add_argument(
        "--cluster_sigma_pc",
        type=float,
        default=4.0,
        help="Cluster Gaussian radius [pc]",
    )
    p.add_argument("--n_field_stars", type=int, default=1, help="Isolated stars")
    p.add_argument(
        "--min_star_mass", type=float, default=10.0, help="Lowest star mass [Msun]"
    )
    p.add_argument(
        "--max_star_mass", type=float, default=40.0, help="Highest star mass [Msun]"
    )
    p.add_argument("--boxsize", type=float, default=4.0, help="Box side [kpc]")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("-o", dest="out", type=str, default="ICs_thin_disc.hdf5")
    return p.parse_args()


def main():
    opt = parse_options()
    rng = np.random.default_rng(opt.seed)

    # Internal units: 1e10 Msun, kpc, km/s.
    unit_mass = 1.988409870698051e43
    unit_length = KPC_CM
    unit_vel = 1e5
    unit_time = unit_length / unit_vel
    msun = MSUN_G / unit_mass
    pc = 1e-3

    L = opt.boxsize
    centre = np.full(3, 0.5 * L)
    Rd = opt.scale_length_pc * pc
    hz = opt.scale_height_pc * pc
    R_max = opt.truncation * Rd
    if R_max > 0.5 * L:
        raise SystemExit(
            "The disc does not fit in the box: lower --truncation or raise --boxsize"
        )

    n_gas = int(round(opt.disc_mass / opt.gas_mass))

    # Radial CDF of an exponential disc, truncated: M(<R) ~ 1 - (1 + R/Rd) exp(-R/Rd).
    r_grid = np.linspace(0.0, R_max, 20001)
    cdf = 1.0 - (1.0 + r_grid / Rd) * np.exp(-r_grid / Rd)
    cdf /= cdf[-1]

    def sample_disc(n):
        R = np.interp(rng.random(n), cdf, r_grid)
        phi = 2.0 * np.pi * rng.random(n)
        z = rng.normal(0.0, hz, n)
        return np.stack([R * np.cos(phi), R * np.sin(phi), z], axis=1)

    pos_gas = sample_disc(n_gas)
    # Remove the vertical tail beyond 4 sigma so the gas stays thin.
    keep = np.abs(pos_gas[:, 2]) < 4.0 * hz
    pos_gas = pos_gas[keep]
    n_gas = len(pos_gas)

    # Local density from the analytic profile, for the initial smoothing length.
    R = np.hypot(pos_gas[:, 0], pos_gas[:, 1])
    sigma0 = opt.disc_mass * msun / (2.0 * np.pi * Rd**2)
    rho = (
        sigma0
        * np.exp(-R / Rd)
        / (np.sqrt(2.0 * np.pi) * hz)
        * np.exp(-0.5 * (pos_gas[:, 2] / hz) ** 2)
    )
    m_gas = opt.gas_mass * msun
    h_gas = 1.2348 * (m_gas / np.maximum(rho, 1e-3 * rho.max())) ** (1.0 / 3.0)

    x_h = 0.752
    mu = 4.0 / (1.0 + 3.0 * x_h)
    u_cgs = K_B * opt.temperature / (M_P_G * (5.0 / 3.0 - 1.0) * mu)
    u_gas = np.full(n_gas, u_cgs / unit_vel**2)

    n_h_centre = (
        sigma0 / (np.sqrt(2.0 * np.pi) * hz) * unit_mass / unit_length**3 * x_h / M_P_G
    )

    # Stars: compact clusters inside 2 scale lengths plus field stars from the gas.
    star_pos = []
    for _ in range(opt.n_clusters):
        r_c = Rd * rng.uniform(0.3, 1.5)
        phi_c = 2.0 * np.pi * rng.random()
        c = np.array([r_c * np.cos(phi_c), r_c * np.sin(phi_c), 0.0])
        star_pos.append(
            c
            + rng.normal(0.0, opt.cluster_sigma_pc * pc, (opt.stars_per_cluster, 3))
            * np.array([1.0, 1.0, 0.5])
        )
    field = sample_disc(20 * max(opt.n_field_stars, 1))
    field = field[np.hypot(field[:, 0], field[:, 1]) < 2.0 * Rd][: opt.n_field_stars]
    star_pos.append(field)
    star_pos = np.vstack(star_pos)
    n_star = len(star_pos)

    # Salpeter slope between the mass limits (dN/dM ~ M^-2.35).
    a = 1.0 - 2.35
    u = rng.random(n_star)
    m_star_msun = (
        opt.min_star_mass**a + u * (opt.max_star_mass**a - opt.min_star_mass**a)
    ) ** (1.0 / a)
    m_star = m_star_msun * msun
    h_star = np.full(n_star, float(np.median(h_gas)))

    pos_gas = pos_gas + centre
    star_pos = star_pos + centre
    ids_gas = np.arange(1, n_gas + 1)
    ids_star = np.arange(n_gas + 1, n_gas + n_star + 1)

    print(f"Gas particles: {n_gas}, mass {opt.gas_mass} Msun each")
    print(f"Stars: {n_star}, total {m_star_msun.sum():.1f} Msun")
    print(f"Central hydrogen number density: {n_h_centre:.2f} cm^-3")

    with h5py.File(opt.out, "w") as f:
        g = f.create_group("/Header")
        g.attrs["BoxSize"] = [L, L, L]
        g.attrs["NumPart_Total"] = [n_gas, 0, 0, 0, n_star, 0]
        g.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
        g.attrs["NumPart_ThisFile"] = [n_gas, 0, 0, 0, n_star, 0]
        g.attrs["Time"] = 0.0
        g.attrs["NumFileOutputsPerSnapshot"] = 1
        g.attrs["MassTable"] = [0.0] * 6
        g.attrs["Flag_Entropy_ICs"] = [0] * 6
        g.attrs["Dimension"] = 3

        g = f.create_group("/Units")
        g.attrs["Unit length in cgs (U_L)"] = unit_length
        g.attrs["Unit mass in cgs (U_M)"] = unit_mass
        g.attrs["Unit time in cgs (U_t)"] = unit_time
        g.attrs["Unit current in cgs (U_I)"] = 1.0
        g.attrs["Unit temperature in cgs (U_T)"] = 1.0

        g = f.create_group("/PartType0")
        g.create_dataset("Coordinates", data=pos_gas, dtype="d")
        g.create_dataset("Velocities", data=np.zeros((n_gas, 3)), dtype="f")
        g.create_dataset("Masses", data=np.full(n_gas, m_gas), dtype="f")
        g.create_dataset("SmoothingLength", data=h_gas, dtype="f")
        g.create_dataset("InternalEnergy", data=u_gas, dtype="f")
        g.create_dataset("ParticleIDs", data=ids_gas, dtype="L")
        g.create_dataset("Densities", data=rho, dtype="f")

        g = f.create_group("/PartType4")
        g.create_dataset("Coordinates", data=star_pos, dtype="d")
        g.create_dataset("Velocities", data=np.zeros((n_star, 3)), dtype="f")
        g.create_dataset("Masses", data=m_star, dtype="f")
        g.create_dataset("ParticleIDs", data=ids_star, dtype="L")
        g.create_dataset("SmoothingLength", data=h_star, dtype="f")
        g.create_dataset("BirthMass", data=m_star, dtype="f")
        g.create_dataset("BirthTime", data=np.zeros(n_star), dtype="f")
        g.create_dataset(
            "StellarParticleType", data=np.zeros(n_star, dtype=int), dtype="i"
        )


if __name__ == "__main__":
    main()
