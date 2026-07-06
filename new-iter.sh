#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") <iter-number> [description]

Creates a new iteration folder in .planning/ with an empty a-requirements.md file.

Arguments:
  iter-number   Required, zero-padded to two digits (e.g. 3 -> 03)
  description   Optional suffix appended to folder name (e.g. "my-analysis")

Examples:
  $(basename "$0") 3              -> .planning/iter03/a-requirements.md
  $(basename "$0") 3 my-analysis  -> .planning/iter03_my-analysis/a-requirements.md
EOF
  exit 1
}

[[ $# -lt 1 ]] && usage

iter_num=$(printf "%02d" "$1")
desc="${2:+_$2}"
folder=".planning/iter${iter_num}${desc}"

if [[ -d "$folder" ]]; then
  echo "Error: $folder already exists." >&2
  exit 1
fi

mkdir -p "$folder"
touch "$folder/a-requirements.md"

echo "Created $folder/a-requirements.md"

file="$folder/a-requirements.md"

if command -v positron &>/dev/null; then
  positron "$file" &
elif command -v vim &>/dev/null; then
  vim "$file"
else
  echo "No editor found. Open $file manually."
fi
