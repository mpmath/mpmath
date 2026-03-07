#!/bin/bash
#
# Test that the version number is provided correctly in a frozen (bundled)
# executable (see #1044).

set -e

# Get the directory of the current repository
MPMATH_DIR="$(realpath "$(dirname "$0")/../../")"
echo "Repo mpmath directory: $MPMATH_DIR"

echo "Install requirements..."
python3 -m pip install build pyinstaller

echo "Building the source distribution from the local repo..."
python3 -m build --sdist

# Find and install the generated tarball
TARBALL=$(ls -t dist/*.tar.gz | head -1)
echo "Generated tarball: $TARBALL"

echo "Installing mpmath from tarball..."
pip install dist/"$(basename $TARBALL)"

TEMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TEMP_DIR"
cd "$TEMP_DIR"

# Create version_script.py that prints the package version
cat << EOF > version_script.py
import mpmath
print(mpmath.__version__)
EOF

# Save local version for later comparison
DIRECT_VERSION="$(python3 -m mpmath --version)"

echo "Building version_script with PyInstaller..."
pyinstaller --onefile --clean version_script.py

echo "Run frozen executable and extract the version from the output..."
FROZEN_VERSION="$(./dist/version_script 2>&1)"

if [ "$DIRECT_VERSION" == "$FROZEN_VERSION" ]; then
    echo "Test passed: Version matches in frozen (bundled) executable."
else
    echo "Test failed: Version mismatch."
    echo "Direct version: $DIRECT_VERSION"
    echo "Frozen version: $FROZEN_VERSION"
    exit 1
fi
