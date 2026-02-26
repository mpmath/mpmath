#!/usr/bin/env python3
"""Regression test: Simulates frozen app behavior where package metadata is unavailable."""
import sys
import importlib.metadata
import importlib

def test_issue_metadata_fallback(monkeypatch):
    # Import mpmath normally first to get the expected version
    import mpmath
    expected_version = mpmath.__version__

    # Save original version function
    original_version = importlib.metadata.version

    def mock_version(package_name):
        if package_name == "mpmath":
            raise importlib.metadata.PackageNotFoundError(
                f"No package metadata was found for {package_name}"
            )
        return original_version(package_name)

    monkeypatch.setattr(importlib.metadata, "version", mock_version)

    # Remove mpmath from sys.modules to force reload with mocked metadata
    sys.modules.pop("mpmath", None)

    # Reload mpmath - should use fallback version
    #   also ensures the hardcoded __version__ attribute is updated correctly
    import mpmath
    assert mpmath.__version__ == expected_version
