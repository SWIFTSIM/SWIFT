import numpy as np


class Config:
    # Impact geometry
    impactor_radius = 1e-3  # m
    impactor_angle = 90.0  # degrees from horizontal (i.e. 90 is vertical impact)
    impactor_speed = 1e2  # m/s
    impactor_x_offset = 0  # m
    impactor_z_offset = 1e-3  # m (0 corresponts to touching surface)

    # Resolution
    particles_per_projectile_radius = 5
    lattice_type = "cubic"  # "cubic" or "bcc"

    # Domain
    approx_surface_width = 20e-3  # m
    approx_surface_depth = 4e-3  # m
    approx_highres_width = 10e-3  # m
    approx_highres_depth = 2e-3  # m
    box_height = 4e-3  # m (height of box above surface)

    # Transitions
    N_transitions = 1

    # Material: surface
    mat_id_surface = 501
    density_surface = 1000  # kg/m^3
    u_surface = 0.0  # J/kg

    # Material: impactor
    mat_id_impactor = 500
    density_impactor = 1000  # kg/m^3
    u_impactor = 0.0  # J/kg

    @property
    def impactor_angle_rad(self) -> float:
        return np.deg2rad(self.impactor_angle)
