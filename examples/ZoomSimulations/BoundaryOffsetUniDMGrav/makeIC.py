"""Generate boundary-straddling ICs for zoom shift-path testing.

This script reuses UniformDMGravity IC generation, then applies a half-box
coordinate shift directly to both PartType1 and PartType2 so the high-resolution
region straddles all periodic corners in the IC file.
"""

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
SOURCE_DIR = HERE.parent / "UniformDMGravity"
SOURCE_SCRIPT = SOURCE_DIR / "makeIC.py"
SOURCE_IC = SOURCE_DIR / "zoom_uniform_dm_gravity.hdf5"
TARGET_IC = HERE / "boundary_offset_uni_dm_grav.hdf5"


def shift_coordinates(path: Path, delta: float = 0.5) -> None:
    """Apply periodic coordinate shift to PartType1 and PartType2."""
    with h5py.File(path, "r+") as handle:
        for ptype in ("PartType1", "PartType2"):
            dset: Any = handle[f"{ptype}/Coordinates"]
            coords = dset[:]
            coords = np.mod(coords + delta, 1.0)
            dset[:] = coords


def main() -> None:
    """Build source ICs, copy locally, and apply boundary-straddling shift."""
    if not SOURCE_SCRIPT.exists():
        raise FileNotFoundError(f"Missing source IC generator: {SOURCE_SCRIPT}")

    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=SOURCE_DIR)
    shutil.copy2(SOURCE_IC, TARGET_IC)
    shift_coordinates(TARGET_IC, delta=0.5)


if __name__ == "__main__":
    main()
