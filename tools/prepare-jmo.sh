#!/usr/bin/env bash
# Repack a built .jmo into the two CPU artifacts one jamovi series needs.
# Metadata-only: it rewrites jamovi's compatibility stamp and rebuilds no R code.
#
#   bash tools/prepare-jmo.sh current 4.6.0
#   bash tools/prepare-jmo.sh solid 4.5.0 path/to/<module>_<version>.jmo
#
# TWO artifacts per series cover every operating system. jamovi's only
# compatibility gate is an exact string compare of the artifact's rVersion
# against the app's JAMOVI_R_VERSION, and that string carries the R version and
# the CPU but no OS component.
#
# Safe only for a module with no compiled code; the guard below enforces that.
# A module that bundles compiled code needs a real build on each CPU instead.
set -euo pipefail
source "$(dirname "$0")/_module.sh"

SERIES="${1:-}"
R_VERSION="${2:-}"
SOURCE="${3:-$MODULE_DIR/${MODULE}_${VERSION}.jmo}"

usage() { echo "usage: prepare-jmo.sh {solid|current} R_VERSION [source.jmo]" >&2; }

case "$SERIES" in
  solid|current) ;;
  *) usage; exit 1 ;;
esac

case "$R_VERSION" in
  [0-9]*.[0-9]*.[0-9]*) ;;
  *) usage; exit 1 ;;
esac

[ -f "$SOURCE" ] || { echo "error: source artifact not found: $SOURCE" >&2; exit 1; }
command -v unzip >/dev/null || { echo "error: unzip is required" >&2; exit 1; }
command -v zip >/dev/null || { echo "error: zip is required" >&2; exit 1; }

# A repack cannot fix compiled code: the .so/.dylib/.dll inside would be for the
# build machine's CPU whatever the stamp claims.
if unzip -l "$SOURCE" | grep -Eq '\.(so|dylib|dll)$'; then
  echo "error: $SOURCE contains native libraries and cannot be repacked" >&2
  exit 1
fi

# The R version in the stamp is a promise about the R the module was BUILT
# under; R warns "built under R version x.y.z" when an older R loads it. Build
# against the target series (tools/build-jmo.sh) rather than restamping.
BUILT="$(unzip -p "$SOURCE" "$MODULE/R/$MODULE/DESCRIPTION" | sed -n 's/^Built: R \([0-9.]*\);.*/\1/p')"
if [ -n "$BUILT" ] && [ "$BUILT" != "$R_VERSION" ]; then
  echo "!! warning: $(basename "$SOURCE") was built under R $BUILT but is being" >&2
  echo "!! stamped for R $R_VERSION. It will load with a warning on the target." >&2
fi

mkdir -p "$ROOT/dist"

for ARCH in x64 arm64; do
  TARGET="${R_VERSION}-${ARCH}"
  OUT="$ROOT/dist/${MODULE}_${VERSION}_${SERIES}_R${R_VERSION}_${ARCH}.jmo"
  # scratch space inside the project, not the system /tmp
  TMP="$(mktemp -d "$ROOT/.tmp/${MODULE}-jmo.XXXXXX")"

  unzip -q "$SOURCE" -d "$TMP"
  for META in "$TMP/$MODULE/jamovi.yaml" "$TMP/$MODULE/jamovi-full.yaml"; do
    [ -f "$META" ] || { echo "error: metadata missing: $META" >&2; rm -rf "$TMP"; exit 1; }
    perl -0pi -e "s/^rVersion: .*\$/rVersion: $TARGET/m" "$META"
    grep -qx "rVersion: $TARGET" "$META" || {
      echo "error: could not set rVersion in $META" >&2; rm -rf "$TMP"; exit 1; }
  done

  rm -f "$OUT"
  ( cd "$TMP" && zip -q -r "$OUT" "$MODULE" )
  rm -rf "$TMP"
  echo ">> wrote dist/$(basename "$OUT")  (rVersion: $TARGET)"
done
