from swiftsimio import load
from swiftsimio.visualisation.projection import project_pixel_grid
from swiftsimio.visualisation.smoothing_length import generate_smoothing_lengths
import numpy as np

data = load("snapshots_e6/snap_128_0031.hdf5")

# Generate smoothing lengths for the dark matter
data.dark_matter.smoothing_length = generate_smoothing_lengths(
    data.dark_matter.coordinates,
    data.metadata.boxsize,
    kernel_gamma=1.8,
    neighbours=57,
    speedup_fac=2,
    dimension=3,
)

rho_est = data.dark_matter.masses/(data.dark_matter.smoothing_length**3)

print(data.dark_matter.smoothing_length.max(), data.dark_matter.smoothing_length.min())
print(rho_est.max())
print(rho_est.mean())
print(np.percentile(rho_est, [50, 90, 99]))