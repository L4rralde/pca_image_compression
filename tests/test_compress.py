import sys, glob
sys.path.append("../python/compress")

import utils
import pytest
import numpy as np

@pytest.mark.parametrize(
    "fpath", 
    list(glob.glob(f"{utils.GIT_ROOT}/images/*.pgm"))
)
def test_read_pgm(fpath: str) -> None:
    img = utils.read_pgm(fpath)
    w, h = img.shape
    assert w != 0
    assert h != 0
