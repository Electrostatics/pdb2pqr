#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
tmpfile=$(mktemp)

echo "ruff check --select I --fix '$SCRIPT_DIR'" > "$tmpfile"
echo "ruff format '$SCRIPT_DIR'" >> "$tmpfile"
cat .github/workflows/python-package.yml | grep 'ruff check' | sed -e "s/^\s+//" >> "$tmpfile"

echo "Run these commands:"
cat "$tmpfile"

while IFS= read -r command
do
  echo "Command: #${command}#"
  eval "$command"
done < "$tmpfile"

rm -f "$tmpfile"

exit 0
