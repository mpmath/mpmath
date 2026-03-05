#!/bin/bash
# Test that the version number is provided correctly in a frozen (bundled) executable (see #1044).
set -e

# Clean up any previous test artifacts
echo "Cleaning up previous test artifacts..."
rm -rf build dist *.spec .venv

# Get the directory of the current repository
MPMATH_DIR=$(pwd)
# # Working directory independent version: 
# MPMATH_DIR=$(realpath "$(dirname "$(realpath "$0")")/../..")
echo "Repo mpmath directory: $MPMATH_DIR"
TEMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TEMP_DIR"
cd "$TEMP_DIR"

echo "Creating virtual environment..."
python3 -m venv .venv
source .venv/bin/activate


cd "$MPMATH_DIR"
echo "Building tarball from local repo..."
# Build the source distribution
python3 -m pip install build pyinstaller
python3 -m build --sdist

# Find and install the generated tarball
TARBALL=$(ls -t dist/*.tar.gz | head -1)
echo "Generated tarball: $TARBALL"

echo "Installing mpmath from tarball..."
pip install dist/"$(basename $TARBALL)"
cd "$TEMP_DIR" # Return to temporary folder to run the version script


# Create version_script.py that imports mpmath and prints the package version
cat << EOF > version_script.py
import mpmath
print("mpmath version:", mpmath.__version__)
EOF

# Save local version for later comparison
DIRECT_VERSION=$(python -m mpmath --version)

echo "Building version_script with PyInstaller..."
pyinstaller --onefile --clean version_script.py

echo "Run frozen executable and extract the version from the output"
FROZEN_VERSION=$(./dist/version_script 2>&1)

# Return to the original repository directory
cd "$MPMATH_DIR" 
deactivate

if [ "$DIRECT_VERSION" == "$FROZEN_VERSION" ]; then
    echo "Test passed: Version matches in frozen (bundled) executable."
else
    echo "Test failed: Version mismatch."
    echo "Direct version: $DIRECT_VERSION"
    echo "Frozen version: $FROZEN_VERSION"
    exit 1
fi
