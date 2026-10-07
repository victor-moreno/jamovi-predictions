# Sourced by every tools/ script: find the module and read its metadata, so the
# scripts can be copied into any module unchanged.
#
# The module is the repository root when the root holds DESCRIPTION, otherwise
# the one subdirectory holding both DESCRIPTION and jamovi/ (e.g. jmvplus/).
# Set MODULE_DIR=/path to override, e.g. when there is more than one.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [ -z "${MODULE_DIR:-}" ]; then
  if [ -f "$ROOT/DESCRIPTION" ]; then
    MODULE_DIR="$ROOT"
  else
    for d in "$ROOT"/*/; do
      if [ -f "$d/DESCRIPTION" ] && [ -d "$d/jamovi" ]; then
        [ -z "${MODULE_DIR:-}" ] || { echo "error: several modules under $ROOT; set MODULE_DIR" >&2; exit 1; }
        MODULE_DIR="${d%/}"
      fi
    done
  fi
fi
[ -n "${MODULE_DIR:-}" ] && [ -f "$MODULE_DIR/DESCRIPTION" ] || {
  echo "error: no module (DESCRIPTION + jamovi/) found under $ROOT" >&2; exit 1; }

MODULE="$(awk -F': *' '$1 == "Package" { print $2; exit }' "$MODULE_DIR/DESCRIPTION")"
VERSION="$(awk -F': *' '$1 == "Version" { print $2; exit }' "$MODULE_DIR/DESCRIPTION")"
[ -n "$MODULE" ] && [ -n "$VERSION" ] || { echo "error: no Package/Version in DESCRIPTION" >&2; exit 1; }

# scratch space inside the project, never the system /tmp
mkdir -p "$ROOT/.tmp"
