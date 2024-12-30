"""Test Sphinx compilation of documentation."""

from pathlib import Path
from sphinx.application import Sphinx


SOURCE_DIR = Path("docs/source")
CONFIG_DIR = SOURCE_DIR


def test_html_documentation(tmp_path):
    """Test Sphinx build of HTML documentation."""
    source_dir = SOURCE_DIR
    config_dir = CONFIG_DIR
    output_dir = tmp_path
    doctree_dir = tmp_path / "doctrees"
    app = Sphinx(
        source_dir,
        config_dir,
        output_dir,
        doctree_dir,
        buildername="html",
        warningiserror=True,
    )
    app.build(force_all=True)


def test_text_documentation(tmp_path):
    """Test Sphinx build of text documentation."""
    source_dir = SOURCE_DIR
    config_dir = CONFIG_DIR
    output_dir = tmp_path
    doctree_dir = tmp_path / "doctrees"
    app = Sphinx(
        source_dir,
        config_dir,
        output_dir,
        doctree_dir,
        buildername="text",
        warningiserror=True,
    )
    app.build(force_all=True)
