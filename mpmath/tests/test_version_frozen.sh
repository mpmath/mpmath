#!/bin/bash
set -e

# Clean up any previous test artifacts
echo "Cleaning up previous test artifacts..."
rm -rf test_local build dist *.spec .venv

mkdir -p test_local
cd test_local

echo "Creating virtual environment..."
python3 -m venv .venv
source .venv/bin/activate


cd .. # Leave the test_local directory to build the tarball from the local repo
echo "Building tarball from local repo..."
# Build the source distribution
python3 -m pip install --quiet build
python3 -m build --sdist

# Find and install the generated tarball
TARBALL=$(ls -t dist/*.tar.gz | head -1)
echo "Generated tarball: $TARBALL"

# Install mpmath from the generated tarball
echo "Installing mpmath from tarball..."
pip install --quiet dist/"$(basename $TARBALL)"
cd test_local # Return to test_local to run the version script


# Create version_script.py that imports mpmath and prints the package version
cat << EOF > version_script.py
import mpmath
print("mpmath version:", mpmath.__version__)
EOF

# Save version_script.py output to DIRECT_VERSION for later comparison
DIRECT_VERSION=$(python version_script.py | grep "mpmath version:" | awk '{print $3}')

# Build with PyInstaller
echo "Building version_script with PyInstaller..."
pip install --quiet pyinstaller
pyinstaller --onefile --clean version_script.py

# Run the frozen executable and capture output
echo "Running frozen executable..."
./dist/version_script 2>&1 | tee output.txt

# Extract the version from the output
FROZEN_VERSION=$(grep "mpmath version:" output.txt | awk '{print $3}')

cd ..
deactivate


# Assert that the version from the frozen/bundled executable matches the direct version
if [ "$DIRECT_VERSION" == "$FROZEN_VERSION" ]; then
    echo "Test passed: Version matches in frozen (bundled) executable."
else
    echo "Test failed: Version mismatch."
    echo "Direct version: $DIRECT_VERSION"
    echo "Frozen version: $FROZEN_VERSION"
    exit 1
fi

