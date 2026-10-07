#!/usr/bin/env bash
# Build every artifact a release needs and publish them as ONE GitHub release.
#
#   bash tools/release.sh            # publish (or replace) the v<version> release
#   bash tools/release.sh --prune    # ... and delete every other release and tag
#
# Four artifacts cover the supported matrix: two jamovi series (solid, current)
# x two CPUs. The tag is the module version, not an R version, because one
# release now carries every target -- so --prune is how the old
# v<version>-R<rversion> tags get cleared out.
set -euo pipefail

PRUNE=no
case "${1:-}" in
  "") ;;
  --prune) PRUNE=yes ;;
  *) echo "usage: release.sh [--prune]" >&2; exit 1 ;;
esac

command -v gh >/dev/null || { echo "error: gh CLI is required" >&2; exit 1; }

source "$(dirname "$0")/_module.sh"

REPO="$(cd "$ROOT" && gh repo view --json nameWithOwner -q .nameWithOwner)"
TAG="v${VERSION}"

# Exactly one file per (series, CPU); the R version in the name comes from the
# app, so it is read back from the filename rather than written down twice.
artifact() {
  local found=( "$ROOT"/dist/"${MODULE}_${VERSION}_$1"_R*_"$2".jmo )
  [ "${#found[@]}" -eq 1 ] && [ -f "${found[0]}" ] || {
    echo "error: expected one $1/$2 artifact in dist/, found ${#found[@]}" >&2; exit 1; }
  printf '%s' "${found[0]}"
}

# clinpred needs jamovi current (>= 28.4, File option): no solid build
for SERIES in current; do
  rm -f "$ROOT"/dist/"${MODULE}_${VERSION}_${SERIES}"_R*.jmo
  bash "$ROOT/tools/build-jmo.sh" "$SERIES"
done

CURR_ARM="$(artifact current arm64)"; CURR_X64="$(artifact current x64)"
ASSETS=("$CURR_ARM" "$CURR_X64")

rver() { echo "$1" | sed -n 's/.*_R\([0-9.]*\)_[a-z0-9]*\.jmo$/\1/p'; }
CURR_R="$(rver "$CURR_ARM")"

echo
echo ">> releasing $TAG on $REPO"
printf '   %s\n' "${ASSETS[@]##*/}"

NOTES="$(cat <<NOTESEOF
$MODULE $VERSION for jamovi current (28.4 or later). Two files, one per CPU:

| your jamovi | Apple silicon | Intel / AMD |
| --- | --- | --- |
| **current** (bundles R $CURR_R) | \`$(basename "$CURR_ARM")\` | \`$(basename "$CURR_X64")\` |

The same file works on macOS, Windows and Linux: jamovi checks the R version
and the CPU, not the operating system. jamovi reports its R version under
Help -> About if you are unsure which you have.

Install: **Modules -> jamovi library -> Sideload**, then select the \`.jmo\`.
NOTESEOF
)"

gh release delete "$TAG" --repo "$REPO" --cleanup-tag --yes >/dev/null 2>&1 || true
gh release create "$TAG" \
  --repo "$REPO" \
  --title "$MODULE $VERSION" \
  --notes "$NOTES" \
  "${ASSETS[@]}"

if [ "$PRUNE" = yes ]; then
  for OLD in $(gh release list --repo "$REPO" --limit 100 --json tagName -q '.[].tagName'); do
    [ "$OLD" = "$TAG" ] && continue
    echo ">> deleting superseded release $OLD"
    gh release delete "$OLD" --repo "$REPO" --cleanup-tag --yes
  done
fi

echo ">> done: https://github.com/$REPO/releases/tag/$TAG"
