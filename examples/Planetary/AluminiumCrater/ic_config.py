import numpy as np


class Config:
    # Impact geometry
    impactor_radius = 0.5e-3  # m
    impactor_angle = 50.0  # degrees from horizontal (i.e. 90 is vertical impact)
    impactor_speed = 5e3  # m/s
    impactor_x_offset = 0.5e-3  # m
    impactor_z_offset = 0.25e-3  # m (0 corresponts to touching surface)

    # Resolution
    particles_per_projectile_radius = 5
    lattice_type = "bcc"  # "cubic" or "bcc"

    # Domain
    approx_surface_width = 20e-3  # m
    approx_surface_depth = 10e-3  # m
    approx_highres_width = 10e-3  # m
    approx_highres_depth = 2e-3  # m
    box_height = 20e-3  # m (height of box above surface)

    # Transitions
    N_transitions = 3
    approx_intermediate_layer_width = 2  # units of particle separation

    # Material: surface
    mat_id_surface = 1000
    density_surface = 2700.0  # kg/m^3
    u_surface = 1e5  # J/kg

    # Material: impactor
    mat_id_impactor = 1000
    density_impactor = 2700.0  # kg/m^3
    u_impactor = 1e5  # J/kg

    @property
    def impactor_angle_rad(self) -> float:
        return np.deg2rad(self.impactor_angle)
