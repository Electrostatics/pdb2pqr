"""Test basic testing functionality."""
import pytest


class TestTest:
    """Test basic testing functionality."""

    def test_assert(self):
        assert True

    def test_raises(self):
        with pytest.raises(ValueError):
            raise ValueError