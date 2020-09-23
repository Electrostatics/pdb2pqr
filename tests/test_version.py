import pdb2pqr
import re
from pdb2pqr import config


def test_version_exists():
    assert hasattr(pdb2pqr, "__version__")


def test_version():
    assert re.match(
        r"[0-9]+\.[0-9]+", pdb2pqr.__version__
    ) or pdb2pqr.__version__.startswith("0+untagged")
    assert re.match(
        r"[0-9]+\.[0-9]+", config.VERSION
    ) or config.VERSION.startswith("0+untagged")
