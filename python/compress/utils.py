import subprocess

import numpy as np

GIT_ROOT = subprocess.check_output(["git", "root"]).decode().strip()

def read_pgm(pgmf_path: str) -> np.ndarray:
    """Return a raster of integers from a PGM as a list of lists."""
    pgmf = open(pgmf_path, 'rb')
    assert pgmf.readline().decode() == 'P5\n'
    (width, height) = [int(i) for i in pgmf.readline().split()]
    depth = int(pgmf.readline())
    assert depth <= 255

    raster = []
    for y in range(height):
        row = []
        for y in range(width):
            row.append(ord(pgmf.read(1)))
        raster.append(row)
    return np.array(raster)
