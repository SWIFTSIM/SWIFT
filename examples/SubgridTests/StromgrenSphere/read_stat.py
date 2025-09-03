#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 21:59:21 2025

@author: roduit
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Remplace par le chemin de ton fichier
file_path = "./statistics.txt"

# Lire le fichier en ignorant les lignes de commentaire
data = np.loadtxt(file_path)

# Liste des colonnes pertinentes
columns = [
    "Step", "Time", "Scale Factor", "Redshift", "Total Mass", "Gas Mass",
    "Dark Matter Mass", "Sink Mass", "Stellar Mass", "Black Hole Mass",
    "Gas Metal Mass", "Star Metal Mass", "BH Metal Mass", "Kinetic Energy",
    "Thermal Energy", "Potential Energy", "Radiated Energy", "Gas Entropy",
    "CoM_x", "CoM_y", "CoM_z", "Momentum_x", "Momentum_y", "Momentum_z",
    "AngMom_x", "AngMom_y", "AngMom_z", "BH Accretion Rate", "BH Accreted Mass",
    "BH Subgrid Mass", "Total H Mass", "Molecular H Mass", "Atomic H Mass",
    "Total He Mass", "Magnetic Energy", "DivB Error", "Cross Helicity",
    "Magnetic Helicity", "BH Luminosity", "BH Jet Power"
]

# Convertir en DataFrame
df = pd.DataFrame(data, columns=columns)

# Norme du moment linéaire total
df["Momentum_Total"] = np.sqrt(df["Momentum_x"]**2 + df["Momentum_y"]**2 + df["Momentum_z"]**2)

# Plot de l'évolution du moment linéaire
plt.figure(figsize=(10, 6))
# plt.plot(df["Time"], df["Momentum_Total"], label="Total Momentum", color="blue")
plt.plot(df["Time"], df["Momentum_Total"], label="Total Momentum", color="blue")
# plt.plot(df["Step"], df["Momentum_x"], label="p_x")
# plt.plot(df["Step"], df["Momentum_y"], label="p_y")
# plt.plot(df["Step"], df["Momentum_z"], label="p_z")
plt.xlabel("Time (yr)")
plt.ylabel("Momentum (kg m/s)")
plt.yscale("symlog", linthresh=1e-13)  # Échelle logarithmique si nécessaire
plt.legend()
plt.title("Évolution du moment linéaire total")
plt.grid()
plt.show()
