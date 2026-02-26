"""Generate the DM-only IC used by UniformDMEAGLE.

This wrapper reuses the UniformDMGravity IC generator so the EAGLE example uses
exactly the same deterministic DM layout, then copies the result to the local
IC file name expected by the EAGLE parameter file.

Why this wrapper exists:
- keeps the DM IC construction logic in one place (UniformDMGravity),
- guarantees the EAGLE and gravity examples start from identical DM ICs,
- avoids copy/pasting and drifting IC-generation implementations.

Output:
- ``zoom_uniform_dm_eagle.hdf5`` in this directory.
"""

import shutil
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE_DIR = HERE.parent / "UniformDMGravity"
SOURCE_SCRIPT = SOURCE_DIR / "makeIC.py"
SOURCE_IC = SOURCE_DIR / "zoom_uniform_dm_gravity.hdf5"
TARGET_IC = HERE / "zoom_uniform_dm_eagle.hdf5"


def main() -> None:
    """Generate base ICs with UniformDMGravity and copy to local filename."""

    # Guard against missing source example files.
    if not SOURCE_SCRIPT.exists():
        raise FileNotFoundError(f"Missing source IC generator: {SOURCE_SCRIPT}")

    # Run the canonical IC generator in its own directory.
    subprocess.run([sys.executable, str(SOURCE_SCRIPT)], check=True, cwd=SOURCE_DIR)

    # Copy to the EAGLE-specific file name expected by the YAML file.
    shutil.copy2(SOURCE_IC, TARGET_IC)
    print(f"Copied {SOURCE_IC.name} -> {TARGET_IC.name}")


if __name__ == "__main__":
    main()
