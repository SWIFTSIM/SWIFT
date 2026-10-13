import numpy as np

import fsps


fiducial_norm_constant = 1.65  # value used in Schaye+15
fiducial_normalisation = 1.515e-4 * fiducial_norm_constant

ages = np.linspace(0, 10, 100)

sp = fsps.StellarPopulation(
    imf_type=0,  # Salpeter+55
    zcontinuous=1,
    sfh=3,
    logzsol=0.0,
    add_agb_dust_model=False,
    use_wr_spectra=False,
    imf_upper_limit=100,
    # imf_lower_limit=0.1,
)

sp.set_tabular_sfh(
    age=ages,
    sfr=np.ones(ages.shape),
)

salpeter_mag = sp.get_mags(tage=0.1, bands=['galex_fuv'])

sp = fsps.StellarPopulation(
    imf_type=1,  # Chabrier+03
    zcontinuous=1,
    sfh=3,
    logzsol=0.0,
    add_agb_dust_model=False,
    use_wr_spectra=False,
    imf_upper_limit=100,
    # imf_lower_limit=0.1,
)

sp.set_tabular_sfh(
    age=ages,
    sfr=np.ones(ages.shape),
)

compare_mag = sp.get_mags(tage=0.1, bands=['galex_fuv'])

new_norm_constant = salpeter_mag - compare_mag
new_normalisation = fiducial_normalisation * new_norm_constant 

print("Ratio (Custom / Salpeter):", 1 / new_norm_constant)
print("New normalisation:", new_normalisation) 


