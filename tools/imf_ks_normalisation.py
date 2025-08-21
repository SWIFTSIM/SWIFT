import numpy as np

import fsps

sp = fsps.StellarPopulation(
    imf_type=0,  # Salpeter+55
    compute_vega_mags=False,
    zcontinuous=1,
    sfh=3,
    logzsol=0.0,
    add_agb_dust_model=False,
    use_wr_spectra=False,
    imf_upper_limit=100,
    imf_lower_limit=0.1,
)

sp.set_tabular_sfh(
    age=np.array([0,1,2]),
    sfr=np.array([1,1,1]),
)

salpeter_mag = sp.get_mags(tage=0.1, bands=['galex_fuv'])

sp = fsps.StellarPopulation(
    imf_type=1,  # Chabrier+03
    compute_vega_mags=False,
    zcontinuous=1,
    sfh=3,
    logzsol=0.0,
    add_agb_dust_model=False,
    use_wr_spectra=False,
    imf_upper_limit=100,
    imf_lower_limit=0.1,
)

sp.set_tabular_sfh(
    age=np.array([0,1,2]),
    sfr=np.array([1,1,1]),
)

compare_mag = sp.get_mags(tage=0.1, bands=['galex_fuv'])

print("Ratio (Custom / Salpeter):", 1 / (salpeter_mag - compare_mag))


