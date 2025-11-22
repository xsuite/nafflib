#!/bin/bash
set -euo pipefail; IFS=$'\n\t'

# Choose Python interpreter
if command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    echo "❌ Neither python3 nor python found in PATH."
    exit 1
fi

echo "Using interpreter: $($PYTHON --version)"

NAME=$($PYTHON setup.py --name)
VER=$($PYTHON setup.py --version)


echo "========================================================================"
echo "Tagging $NAME v$VER"
echo "========================================================================"

git tag v$VER
git push origin v$VER

echo "========================================================================"
echo "Releasing $NAME v$VER on PyPI"
echo "========================================================================"

$PYTHON setup.py sdist
twine upload dist/* --verbose --repository nafflib
rm -r dist/ *.egg-info
