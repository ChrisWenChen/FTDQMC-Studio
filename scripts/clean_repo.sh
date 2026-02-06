#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

find . -name '.DS_Store' -type f -delete
find . -name '._*' -type f -delete

# Remove build artifacts
rm -rf build build-* build_* out

echo "Cleaned .DS_Store, ._* files, build*/out directories."
