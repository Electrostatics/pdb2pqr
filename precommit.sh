#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
cd "$SCRIPT_DIR"

if command -v pre-commit &> /dev/null; then
  PRE_COMMIT=pre-commit
elif [[ -x "$SCRIPT_DIR/.venv/bin/pre-commit" ]]; then
  PRE_COMMIT="$SCRIPT_DIR/.venv/bin/pre-commit"
else
  echo "pre-commit is not installed; install the project's dev dependencies first." >&2
  exit 1
fi

exec "$PRE_COMMIT" run --all-files
