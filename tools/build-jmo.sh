#!/usr/bin/env bash
# Build the module against one jamovi series, then repack it for both CPUs.
#
#   bash tools/build-jmo.sh current    # /Applications/jamovi.app
#   bash tools/build-jmo.sh solid      # /Applications/jamovi-solid.app
#
# The .jmo's rVersion stamp and the R package's Built: field both come from the
# R bundled inside the app passed to --home; on macOS jmc ignores --rpath and
# the R on PATH entirely (jamovi-compiler/index.js, the darwin branch). So the
# series alone decides the stamp -- never hand-edit it, and never read the R
# version from anywhere but that app's own env.conf, which is what this does.
#
# Uses jmc --build, not --install: a release build must not disturb the modules
# installed in either app. Use tools/install.sh for that.
set -euo pipefail
source "$(dirname "$0")/_module.sh"

SERIES="${1:-}"
case "$SERIES" in
  current) APP="${JAMOVI_CURRENT_APP:-/Applications/jamovi.app}" ;;
  solid)   APP="${JAMOVI_SOLID_APP:-/Applications/jamovi-solid.app}" ;;
  *) echo "usage: build-jmo.sh {solid|current}" >&2; exit 1 ;;
esac
[ -d "$APP" ] || { echo "error: $APP is not installed" >&2; exit 1; }

STAMP="$(sed -n 's/^JAMOVI_R_VERSION=//p' "$APP/Contents/Resources/env.conf" | tr -d '"[:space:]')"
R_VERSION="${STAMP%-*}"
[ -n "$R_VERSION" ] || { echo "error: no JAMOVI_R_VERSION in $APP" >&2; exit 1; }

SOURCE="$MODULE_DIR/${MODULE}_${VERSION}.jmo"
NODE="$(Rscript --vanilla -e 'cat(system.file("node-darwin", "bin", "node", package = "node"))')"
JMC="$(Rscript --vanilla -e 'cat(system.file("node_modules", "jamovi-compiler", "index.js", package = "jmvtools"))')"
[ -x "$NODE" ] && [ -f "$JMC" ] || { echo "error: jmvtools/node not available to Rscript" >&2; exit 1; }

echo ">> building $MODULE $VERSION against $SERIES ($(basename "$APP"), R $R_VERSION)"
rm -f "$SOURCE"
"$NODE" "$JMC" --build "$MODULE_DIR" --home "$APP" --jmo "$SOURCE"
[ -f "$SOURCE" ] || { echo "error: jmc did not produce $SOURCE" >&2; exit 1; }

# jmc copies the R library cache under <module>/build/ into the artifact
# wholesale, and that cache is reused across builds -- a package renamed in the
# past leaves its old copy behind and ships silently. List what is bundled.
echo ">> bundles: $(unzip -l "$SOURCE" | sed -n "s#.*$MODULE/R/\([^/]*\)/DESCRIPTION#\1#p" | sort -u | tr '\n' ' ')"

bash "$ROOT/tools/prepare-jmo.sh" "$SERIES" "$R_VERSION" "$SOURCE"
