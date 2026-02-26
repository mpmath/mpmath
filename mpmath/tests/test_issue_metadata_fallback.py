#!/usr/bin/env python3
"""Regression test: Simulates frozen app behavior where package metadata is unavailable."""
import sys
import importlib.metadata
import importlib

def test_issue_metadata_fallback(monkeypatch):
    # Mock metadata.version to simulate frozen app where metadata is unavailable
    original_version = importlib.metadata.version

    def mock_version(package_name):
        if package_name == "mpmath":
            raise importlib.metadata.PackageNotFoundError(
                f"No package metadata was found for {package_name}"
            )
        return original_version(package_name)

    monkeypatch.setattr(importlib.metadata, "version", mock_version)

    # Remove mpmath from sys.modules to force reload
    sys.modules.pop("mpmath", None)

    # Import mpmath - should fall back to "0.0.0"
    import mpmath
    assert mpmath.__version__ == "0.0.0"
