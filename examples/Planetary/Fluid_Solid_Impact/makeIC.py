import h5py
import numpy as np
from ic_config import Config


class CubicLattice:
    @staticmethod
    def positions(delta: float, Lx: float, Ly: float, Lz: float) -> np.ndarray:
        """Generate particle positions on a simple cubic lattice.

        The domain is centred in x and y, with the top face at z = 0.
        """
        if not (
            (np.isclose(Lx % delta, 0) or np.isclose(Lx % delta, delta))
            and (np.isclose(Ly % delta, 0) or np.isclose(Ly % delta, delta))
            and (np.isclose(Lz % delta, 0) or np.isclose(Lz % delta, delta))
        ):
            raise ValueError(
                f"Lattice dimensions ({Lx}, {Ly}, {Lz}) must be divisible by delta ({delta})"
            )

        # Generate coords
        x = (np.arange(int(Lx / delta)) + 0.5) * delta - 0.5 * Lx
        y = (np.arange(int(Ly / delta)) + 0.5) * delta - 0.5 * Ly
        z = (np.arange(int(Lz / delta)) + 0.5) * delta - Lz

        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        return np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    @staticmethod
    def separation(distance: float, points_per_distance: int) -> float:
        return distance / points_per_distance


class BCCLattice:
    @staticmethod
    def positions(delta: float, Lx: float, Ly: float, Lz: float) -> np.ndarray:
        """Generate particle positions on a BCC lattice.

        Constructed from two offset cubic sub-lattices, centred in the domain.
        """
        base = CubicLattice.positions(delta, Lx, Ly, Lz)
        return np.vstack([base - 0.25 * delta, base + 0.25 * delta])

    @staticmethod
    def separation(distance: float, points_per_distance: int) -> float:
        return np.cbrt(2) * distance / points_per_distance


LATTICE_TYPES = {"cubic": CubicLattice, "bcc": BCCLattice}


class ImpactICs:
    def __init__(self, cfg: Config):
        self.cfg = cfg

        if cfg.lattice_type not in LATTICE_TYPES:
            raise ValueError(
                f"Unknown lattice type '{cfg.lattice_type}'. Choose from {list(LATTICE_TYPES)}"
            )
        self.lattice = LATTICE_TYPES[cfg.lattice_type]

        self.deltas = self._compute_deltas()
        self.widths, self.depths = self._compute_region_dimensions()

        self.surface: dict | None = None
        self.impactor: dict | None = None

    def _compute_deltas(self) -> list[float]:
        """Particle separations for each resolution level, doubling each time."""
        delta0 = self.lattice.separation(
            self.cfg.impactor_radius, self.cfg.particles_per_projectile_radius
        )
        return [delta0 * 2**i for i in range(self.cfg.N_transitions + 1)]

    def _compute_region_dimensions(self) -> tuple[list[float], list[float]]:
        """Compute width and depth of each nested resolution region."""
        cfg = self.cfg

        if len(self.deltas) == 1:
            if (
                cfg.approx_highres_width != cfg.approx_surface_width
                or cfg.approx_highres_depth != cfg.approx_surface_depth
            ):
                raise ValueError(
                    "With no transition layers, highres and surface dimensions must match"
                )

        # Innermost high-res region. Round up to nearest delta of low res region
        widths = [
            np.ceil(cfg.approx_highres_width / (2 * self.deltas[-1]))
            * (2 * self.deltas[-1])
        ]
        depths = [
            np.ceil(cfg.approx_highres_depth / (self.deltas[-1])) * (self.deltas[-1])
        ]

        # Add intermediate transition layers if needed
        if len(self.deltas) > 2:
            for i in range(1, len(self.deltas) - 1):
                # Rounding up the full region (including inner regions) to nearest delta of low res region
                width = np.ceil(
                    (
                        2 * cfg.approx_intermediate_layer_width * self.deltas[i]
                        + widths[i - 1]
                    )
                    / (2 * self.deltas[i + 1])
                ) * (
                    2 * self.deltas[i + 1]
                )  # 2 here because we add particles on both sides of the domain
                depth = (
                    np.ceil(
                        (
                            cfg.approx_intermediate_layer_width * self.deltas[i]
                            + depths[i - 1]
                        )
                        / self.deltas[i + 1]
                    )
                    * self.deltas[i + 1]
                )
                widths.append(width)
                depths.append(depth)

        # Add final region to reach surface dimensions if needed. Round up to nearest delta
        if len(self.deltas) > 1:
            width = np.ceil(cfg.approx_surface_width / (2 * self.deltas[-1])) * (
                2 * self.deltas[-1]
            )
            depth = (
                np.ceil(cfg.approx_surface_depth / (self.deltas[-1]))
                * (self.deltas[-1])
            )
            widths.append(width)
            depths.append(depth)

        return widths, depths

    def _generate_region(
        self,
        delta: float,
        Lx: float,
        Ly: float,
        Lz: float,
        density: float,
        mat_id: int,
        u: float,
    ) -> dict:
        """Fill a box with lattice particles."""
        pos = self.lattice.positions(delta, Lx, Ly, Lz)
        N = len(pos)
        particle_mass = density * (Lx * Ly * Lz) / N
        return {
            "pos": pos,
            "vel": np.zeros_like(pos),
            "m": np.full(N, particle_mass),
            "rho": np.full(N, density),
            "u": np.full(N, u),
            "mat": np.full(N, mat_id),
        }

    def _mask_inner_region(
        self, region: dict, inner_width: float, inner_depth: float
    ) -> dict:
        """Remove particles that lie inside the next-finer resolution region."""
        pos = region["pos"]
        inside = (
            (np.abs(pos[:, 0]) <= 0.5 * inner_width)
            & (np.abs(pos[:, 1]) <= 0.5 * inner_width)
            & (pos[:, 2] >= -inner_depth)
        )
        return {k: v[~inside] for k, v in region.items()}

    def make_surface(self):
        """Build the multi-resolution target surface."""
        cfg = self.cfg

        # Initialise empty surface dict
        self.surface = {
            "pos": np.empty((0, 3)),
            "vel": np.empty((0, 3)),
            "m": np.empty(0),
            "rho": np.empty(0),
            "u": np.empty(0),
            "mat": np.empty(0),
        }

        # Generate each resolution region and concatenate to surface, masking out inner regions to avoid overlap
        for i, delta in enumerate(self.deltas):
            w, d = self.widths[i], self.depths[i]
            region = self._generate_region(
                delta, w, w, d, cfg.density_surface, cfg.mat_id_surface, cfg.u_surface
            )

            if i > 0:
                region = self._mask_inner_region(
                    region, self.widths[i - 1], self.depths[i - 1]
                )

            for key in self.surface:
                self.surface[key] = np.concatenate([self.surface[key], region[key]])

    def make_impactor(self) -> None:
        """Build the spherical impactor above the surface impact point."""
        cfg = self.cfg
        delta = self.deltas[0]
        r = cfg.impactor_radius

        # Fill a cube then carve out a sphere
        cube_width = np.ceil(3 * r / delta) * delta
        cube = self._generate_region(
            delta,
            cube_width,
            cube_width,
            cube_width,
            cfg.density_impactor,
            cfg.mat_id_impactor,
            cfg.u_impactor,
        )

        # Shift cube so it is centred at the origin in z (already is in x, y)
        cube["pos"][:, 2] += 1.5 * r

        # Make impactor symmetric if BCC
        if cfg.lattice_type == "bcc":
            cube["pos"] += 0.25 * delta

        inside_sphere = np.linalg.norm(cube["pos"], axis=1) <= r
        self.impactor = {k: v[inside_sphere] for k, v in cube.items()}

        # Place impactor with offset
        self.impactor["pos"][:, 0] += cfg.impactor_x_offset
        self.impactor["pos"][:, 2] += r + cfg.impactor_z_offset

        # Apply impact velocity (angle measured from horizontal)
        self.impactor["vel"][:, 0] = -cfg.impactor_speed * np.cos(
            cfg.impactor_angle_rad
        )
        self.impactor["vel"][:, 2] = -cfg.impactor_speed * np.sin(
            cfg.impactor_angle_rad
        )

    def make_ics(self) -> None:
        """Generate both surface and impactor."""
        print("Generating surface...")
        self.make_surface()
        print("Generating impactor...")
        self.make_impactor()

    def save(self, filename: str) -> None:
        """Write initial conditions to an HDF5 file in SWIFT/GADGET format."""
        if self.surface is None or self.impactor is None:
            raise RuntimeError("Call make_ics() before save()")

        surface_width = self.widths[-1]
        surface_depth = self.depths[-1]

        # Combine surface and impactor particle data into single arrays for output
        pos = np.concatenate([self.surface["pos"], self.impactor["pos"]])
        vel = np.concatenate([self.surface["vel"], self.impactor["vel"]])
        m = np.concatenate([self.surface["m"], self.impactor["m"]])
        rho = np.concatenate([self.surface["rho"], self.impactor["rho"]])
        u = np.concatenate([self.surface["u"], self.impactor["u"]])
        mat = np.concatenate([self.surface["mat"], self.impactor["mat"]])
        h = np.cbrt(m / rho)

        # Shift coordinates so positive
        pos[:, 0] += 0.5 * surface_width
        pos[:, 1] += 0.5 * surface_width
        pos[:, 2] += surface_depth

        # Particle IDs: boundary particles get the lowest IDs for SWIFT
        N_surface = len(self.surface["pos"])
        surface_pos = pos[:N_surface]
        boundary_thickness = 0.05 * surface_depth

        is_boundary = (
            (surface_pos[:, 2] <= surface_pos[:, 2].min() + boundary_thickness)
            | (surface_pos[:, 0] <= surface_pos[:, 0].min() + boundary_thickness)
            | (surface_pos[:, 0] >= surface_pos[:, 0].max() - boundary_thickness)
            | (surface_pos[:, 1] <= surface_pos[:, 1].min() + boundary_thickness)
            | (surface_pos[:, 1] >= surface_pos[:, 1].max() - boundary_thickness)
        )
        N_boundary = is_boundary.sum()
        N_interior = N_surface - N_boundary

        ids = np.empty(len(pos))
        ids[:N_surface][is_boundary] = np.arange(1, N_boundary + 1)
        ids[:N_surface][~is_boundary] = np.arange(
            N_boundary + 1, N_boundary + N_interior + 1
        )
        ids[N_surface:] = np.arange(N_surface + 1, len(pos) + 1)

        N_part = len(pos)
        print(f"\nTotal particles : {len(pos):,}")
        print(f"Boundary particles : {N_boundary:,}")

        with h5py.File(filename, "w") as f:
            # Header
            grp = f.create_group("/Header")
            grp.attrs["BoxSize"] = [
                surface_width,
                surface_width,
                surface_depth + cfg.box_height,
            ]
            grp.attrs["NumPart_Total"] = [N_part, 0, 0, 0, 0, 0]
            grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
            grp.attrs["NumPart_ThisFile"] = [N_part, 0, 0, 0, 0, 0]
            grp.attrs["Time"] = 0.0
            grp.attrs["NumFilesPerSnapshot"] = 1
            grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            grp.attrs["Flag_Entropy_ICs"] = 0
            grp.attrs["Dimension"] = 3

            # Units
            grp = f.create_group("/Units")
            grp.attrs["Unit length in cgs (U_L)"] = 100.0
            grp.attrs["Unit mass in cgs (U_M)"] = 1000.0
            grp.attrs["Unit time in cgs (U_t)"] = 1.0
            grp.attrs["Unit current in cgs (U_I)"] = 1.0
            grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

            # Particle group
            grp = f.create_group("/PartType0")
            grp.create_dataset("Coordinates", data=pos, dtype="d")
            grp.create_dataset("Velocities", data=vel, dtype="f")
            grp.create_dataset("Masses", data=m, dtype="f")
            grp.create_dataset("Density", data=rho, dtype="f")
            grp.create_dataset("SmoothingLength", data=h, dtype="f")
            grp.create_dataset("InternalEnergy", data=u, dtype="f")
            grp.create_dataset("ParticleIDs", data=ids, dtype="L")
            grp.create_dataset("MaterialIDs", data=mat, dtype="i")

        print(f"Saved to {filename}")


if __name__ == "__main__":
    cfg = Config()
    ics = ImpactICs(cfg)
    ics.make_ics()
    ics.save("impact.hdf5")
