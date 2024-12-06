import re

import pdb2pqr
from pdb2pqr import __version__


def test_version_exists():
    assert hasattr(pdb2pqr, "__version__")


def test_version():
    assert re.match(r"[0-9]+\.[0-9]+\.[0-9]+", __version__)
