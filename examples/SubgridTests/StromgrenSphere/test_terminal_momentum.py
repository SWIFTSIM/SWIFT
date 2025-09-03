#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 16:54:57 2025

@author: roduit
"""

import numpy as np
import astropy.units as u

def feedback_get_SN_terminal_momentum(gas_number_density, gas_metallicity, star_ejected_energy):
    """
    Compute the terminal momentum of a supernova (SN) blastwave.
    
    Parameters:
    - sp: Object containing feedback data (energy_ejected, weighted_gas_metallicity, weighted_gas_density)
    - phys_const: Object containing physical constants (const_solar_mass, const_proton_mass)
    - us: Object providing unit conversion factors (units_cgs_conversion_factor function)
    
    Returns:
    - Terminal momentum in internal units.
    """
    # Terminal momentum 0 (in internal units)
    p_terminal_0 = (2.5e5 * u.Msun * u.km/u.s)
    
    # Energy ejected in erg
    E_ej = star_ejected_energy.to(u.erg)
    ten_to_51 = 1e51*u.erg
    
    # Velocity factor (currently 1, see PAPER 2024)
    velocity_factor = 1.0
    
    # Metallicity factor
    Z_sun = 0.0134  # Solar metallicity
    
    if gas_metallicity / Z_sun < 0.01:
        metallicity_factor = 2.0
    elif 0.01 <= gas_metallicity / Z_sun <= 1:
        metallicity_factor = (gas_metallicity / Z_sun) ** -0.18
    else:
        metallicity_factor = (gas_metallicity / Z_sun) ** -0.14
    
    
    n_bar = gas_number_density/u.cm**(-3)
    if n_bar < 0.001:
        density_factor = 2.63
    else:
        density_factor = n_bar ** -0.143
    
    # Compute terminal momentum in internal units
    p_terminal = (p_terminal_0 * E_ej / ten_to_51 *
                  density_factor * metallicity_factor * velocity_factor)
    
    
    p_terminal = p_terminal.to(u.Msun*u.km/u.s)
    return p_terminal

gas_number_density = 1*u.cm**(-3)
gas_metallicity = 0.02
star_ejected_energy = 1e51*u.erg

p_term = feedback_get_SN_terminal_momentum(gas_number_density, gas_metallicity, star_ejected_energy)
