import re
import pdb2pqr
from pdb2pqr import __version__, config


def test_version_exists():
    assert hasattr(pdb2pqr, "__version__")


def test_version():
    assert re.match(r"[0-9]+\.[0-9]+", __version__) or __version__.startswith(
        "0+untagged"
    )
    assert re.match(
        r"[0-9]+\.[0-9]+", config.VERSION
    ) or config.VERSION.startswith("0+untagged")
