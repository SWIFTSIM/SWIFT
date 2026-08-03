#!/usr/bin/env python3
################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
"""Generate initial conditions for the Boss & Bodenheimer (1979) collapse test.

A uniform-density, rigidly-rotating sphere of gas is perturbed with an m=2
azimuthal density perturbation:

    rho(phi) = rho_0 * (1 + amplitude * cos(m * phi))

which is imposed by displacing each particle's azimuth rather than by
resampling the density field:

    phi' = phi + amplitude * cos(m * phi) / m

Combined with the rigid rotation, this perturbation grows into the classic
two-arm bar/fragmentation pattern used to validate sink particle formation
from gravitational collapse and fragmentation.
"""
import argparse

import h5py
import numpy as np
from astropy import constants as cte
from astropy import units as u

################################################################################
# Sphere model
################################################################################


class HomogeneousSphere:
    """Sphere of constant density."""

    def __init__(self, R_s, M_tot=None, rho_0=None, G=1.0):
        if M_tot is not None:
            self.M_tot = M_tot
            self.rho_0 = M_tot / (4.0 / 3.0 * np.pi * R_s**3)
        elif rho_0 is not None:
            self.rho_0 = rho_0
            self.M_tot = 4.0 / 3.0 * np.pi * R_s**3 * rho_0
        else:
            raise ValueError("Provide either M_tot or rho_0.")
        self.G = G
        self.R_s = R_s

    def free_fall_time(self):
        return np.sqrt(3 * np.pi / (32 * self.G * self.rho_0))

    def potential_energy(self):
        """Analytic gravitational self-energy of a uniform-density sphere."""
        return -3.0 / 5.0 * self.G * self.M_tot**2 / self.R_s

    def moment_of_inertia(self):
        """Analytic moment of inertia of a uniform-density sphere (about any axis
        through its centre)."""
        return 2.0 / 5.0 * self.M_tot * self.R_s**2

    def generate_particles(self, n, seed):
        """Draw n particles uniformly distributed in the volume of the sphere."""
        rng = np.random.default_rng(seed=seed)

        random1 = rng.random(n)
        random2 = rng.random(n)
        random3 = rng.random(n)

        r = self.R_s * random1 ** (1.0 / 3.0)
        phi = random2 * 2.0 * np.pi
        costh = 1.0 - 2.0 * random3
        sinth = np.sqrt(1.0 - costh**2)

        x = r * sinth * np.cos(phi)
        y = r * sinth * np.sin(phi)
        z = r * costh

        pos = np.transpose(np.array([x, y, z]))
        vel = np.zeros([n, 3])
        mass = np.ones(n) * self.M_tot / n

        return pos, vel, mass


################################################################################
# Perturbations
################################################################################


def apply_azimuthal_density_perturbation(pos, amplitude, m_pert):
    """Impose rho(phi) = rho_0 * (1 + amplitude * cos(m_pert * phi)) by shifting
    each particle's azimuth. r and theta are left unchanged."""
    r = np.linalg.norm(pos, axis=1)
    phi = np.arctan2(pos[:, 1], pos[:, 0])
    with np.errstate(invalid="ignore"):
        cos_theta = np.where(r > 0, pos[:, 2] / np.where(r > 0, r, 1.0), 1.0)
    theta = np.arccos(cos_theta)

    phi_prime = phi + amplitude * np.cos(m_pert * phi) / m_pert

    pos[:, 0] = r * np.cos(phi_prime) * np.sin(theta)
    pos[:, 1] = r * np.sin(phi_prime) * np.sin(theta)
    # pos[:, 2] (= r * cos(theta)) is unchanged: the perturbation is azimuthal only.
    return pos


def boss_bodenheimer_density(pos, rho_0, amplitude, m_pert):
    """Local gas density at each particle, following the same perturbation used
    to displace the particles. Used only to populate the IC's Densities field."""
    phi = np.arctan2(pos[:, 1], pos[:, 0])
    return rho_0 * (1.0 + amplitude * np.cos(m_pert * phi))


def apply_rigid_rotation(pos, Omega):
    """v = Omega x r, with Omega along the z-axis."""
    axis = np.array([0.0, 0.0, 1.0])
    return Omega * np.cross(axis, pos)


################################################################################
# CLI
################################################################################


def parse_options():
    description = """Generate a rotating, azimuthally-perturbed homogeneous
sphere for the Boss & Bodenheimer (1979) collapse and fragmentation test."""
    epilog = """
Examples
--------
python3 makeIC.py -o boss_bodenheimer.hdf5
python3 makeIC.py -o boss_bodenheimer.hdf5 -N 200000 --density_amplitude 0.5 --m 2
"""
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-o",
        dest="outputfilename",
        default="boss_bodenheimer.hdf5",
        help="Name of the output file.",
    )
    parser.add_argument(
        "-N",
        dest="N",
        type=int,
        default=50000,
        help="Number of gas particles. The reference Boss & Bodenheimer papers "
        "use resolutions well beyond what is practical for a quick example; "
        "raise this for production runs.",
    )
    parser.add_argument(
        "--R_s",
        type=float,
        default=2406.0,
        help="Radius of the sphere, in AU.",
    )
    parser.add_argument(
        "--M_tot",
        type=float,
        default=1.0,
        help="Total gas mass, in solar masses.",
    )
    parser.add_argument(
        "--density_amplitude",
        type=float,
        default=0.5,
        help="Amplitude of the m=2 azimuthal density perturbation.",
    )
    parser.add_argument(
        "--m",
        type=int,
        default=2,
        help="Azimuthal wavenumber of the density perturbation (cos(m*phi)).",
    )
    parser.add_argument(
        "--Omega",
        type=float,
        default=1.6e-12,
        help="Rigid rotation rate about the z-axis, in s^-1.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.26,
        help="Ratio of thermal to gravitational energy, alpha = U_therm / |E_pot|. "
        "0.26 is the standard Boss & Bodenheimer value.",
    )
    parser.add_argument(
        "--boxsize",
        type=float,
        default=None,
        help="Side length of the (cubic) box, in AU. Defaults to 2.5 times the "
        "maximal particle distance from the centre.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed.",
    )

    return parser.parse_args()


################################################################################
# Main
################################################################################

if __name__ == "__main__":
    opt = parse_options()

    # Units: AU / solar mass / km/s (matching the InternalUnitSystem in params.yml)
    u_Length = 1 * u.AU
    u_Mass = 1 * u.M_sun
    u_Velocity = 1 * u.km / u.s
    u_Time = u_Length.to(u.km) / u_Velocity

    UnitLength_in_cgs = u_Length.to(u.cm).value
    UnitMass_in_cgs = u_Mass.to(u.g).value
    UnitVelocity_in_cgs = u_Velocity.to(u.cm / u.s).value
    UnitTime_in_cgs = u_Time.to(u.s).value

    G = cte.G.to(u_Length**3 / (u_Mass * u_Time**2)).value

    R_s = opt.R_s
    M_tot = (opt.M_tot * u.M_sun).to(u_Mass).value
    N = opt.N

    sphere = HomogeneousSphere(R_s=R_s, M_tot=M_tot, G=G)
    t_ff = sphere.free_fall_time()
    E_pot = sphere.potential_energy()

    print("Boss & Bodenheimer collapse test")
    print("---------------------------------")
    print(f"  N          = {N}")
    print(f"  R_s        = {R_s:.3e} AU")
    print(f"  M_tot      = {M_tot:.3e} M_sun")
    print(f"  rho_0      = {sphere.rho_0:.3e} (code units)")
    print(f"  t_ff       = {t_ff:.3e} (code units)")
    print(f"  E_pot      = {E_pot:.3e} (code units, analytic uniform-sphere value)")

    pos, vel, mass = sphere.generate_particles(N, seed=opt.seed)

    # Azimuthal density perturbation (displaces particles; recomputed density is
    # only used to populate the IC's Densities field).
    pos = apply_azimuthal_density_perturbation(
        pos, amplitude=opt.density_amplitude, m_pert=opt.m
    )
    rho = boss_bodenheimer_density(
        pos, rho_0=sphere.rho_0, amplitude=opt.density_amplitude, m_pert=opt.m
    )

    # Rigid rotation
    Omega_code = (opt.Omega / u.s).to(1 / u_Time).value
    vel = apply_rigid_rotation(pos, Omega_code)

    I = sphere.moment_of_inertia()
    E_rot = 0.5 * I * Omega_code**2
    print(f"  Omega      = {Omega_code:.3e} (code units)")
    print(f"  E_rot      = {E_rot:.3e} (code units)")
    print(f"  beta       = E_rot/|E_pot| = {E_rot / abs(E_pot):.3e}")

    # Thermal energy from the target virial ratio alpha = U_therm / |E_pot|
    U_int_total = opt.alpha * abs(E_pot)
    U_int_per_part = np.full(N, U_int_total / M_tot)
    print(f"  alpha      = U_therm/|E_pot| = {opt.alpha:.3e}")
    print(
        f"  This run uses the isothermal EoS: set params.yml's "
        f"EoS:isothermal_internal_energy to {U_int_per_part[0]:.6e} to match."
    )

    # Boxsize and smoothing length
    r_max = np.max(np.linalg.norm(pos, axis=1))
    L = opt.boxsize if opt.boxsize is not None else 2.5 * r_max
    h = np.full(N, R_s / N ** (1.0 / 3.0))

    print(f"  boxsize    = {L:.3e} AU")
    print(f"  shift      = {L / 2:.3e} AU (set this in params.yml)")

    ids = np.arange(N, dtype=np.uint64)

    with h5py.File(opt.outputfilename, "w") as f:
        grp = f.create_group("/Header")
        grp.attrs["BoxSize"] = [L, L, L]
        grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
        grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
        grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
        grp.attrs["Time"] = 0.0
        grp.attrs["NumFileOutputsPerSnapshot"] = 1
        grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        grp.attrs["Flag_Entropy_ICs"] = 0
        grp.attrs["Dimension"] = 3

        grp = f.create_group("/Units")
        grp.attrs["Unit length in cgs (U_L)"] = UnitLength_in_cgs
        grp.attrs["Unit mass in cgs (U_M)"] = UnitMass_in_cgs
        grp.attrs["Unit time in cgs (U_t)"] = UnitTime_in_cgs
        grp.attrs["Unit current in cgs (U_I)"] = 1.0
        grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

        grp = f.create_group("/PartType0")
        grp.create_dataset("Coordinates", data=pos, dtype="d")
        grp.create_dataset("Velocities", data=vel, dtype="f")
        grp.create_dataset("Masses", data=mass, dtype="f")
        grp.create_dataset("SmoothingLength", data=h, dtype="f")
        grp.create_dataset("InternalEnergy", data=U_int_per_part, dtype="f")
        grp.create_dataset("ParticleIDs", data=ids, dtype="L")
        grp.create_dataset("Densities", data=rho, dtype="f")

    print(f"\n{opt.outputfilename} saved.")
