#!/usr/bin/env bash
# Build and install the module into jamovi desktop and/or a running jamovi
# Docker container, then run tools/smoke.R (if present) against the install.
#
#   bash tools/install.sh              both targets, whichever are available
#   bash tools/install.sh desktop
#   bash tools/install.sh docker [container]     (default container: jamovi)
#
# The desktop target uses whichever R `Rscript` resolves to (respecting
# ~/.Rprofile, which may append jamovi.app's bundled module library). With rig,
# pick the active R with `rig default <version>` first: it must match the R that
# jamovi.app bundles, or rlang (behind jmvcore) fails to load.
set -euo pipefail
source "$(dirname "$0")/_module.sh"

TARGET="${1:-both}"
CONTAINER="${2:-jamovi}"
ARTIFACT="$MODULE_DIR/${MODULE}_${VERSION}.jmo"
SMOKE="$ROOT/tools/smoke.R"

# R preamble for the smoke test: jamovi's base library, the installed module and
# jmv's library first, as jamovi's engine sets them up (jmc does not bundle
# packages jmv already ships, e.g. car, so a module importing one needs jmv's).
smoke_preamble() {
  printf '.libPaths(c("%s", "%s", "%s", .libPaths()))\n' "$1" "$2" "$(dirname "$(dirname "$1")")/jmv/R"
  printf 'suppressPackageStartupMessages(library("%s", character.only = TRUE))\n' "$MODULE"
}

# ── desktop ──────────────────────────────────────────────────────────────────
install_desktop() {
  local APP LOG MODDIR OLD i
  APP="${JAMOVI_CURRENT_APP:-/Applications/jamovi.app}"
  [ -d "$APP" ] || { echo "!! desktop: $APP is not installed — skipping"; return 0; }

  MODDIR="$HOME/Library/Application Support/jamovi/modules/$MODULE"
  # the Built: stamp of any copy already installed; the new one must differ
  OLD="$(grep '^Built:' "$MODDIR/R/$MODULE/DESCRIPTION" 2>/dev/null || true)"

  echo ">> desktop: building $MODULE $VERSION with $(Rscript -e 'cat(R.version.string)')"
  LOG="$(mktemp "$ROOT/.tmp/install.XXXXXX")"
  ( cd "$MODULE_DIR" && Rscript -e 'jmvtools::install()' ) 2>&1 | tee "$LOG" | grep -vE '^\s*$' || true

  # jmvtools::install() reports errors on stdout and still exits 0, and it can
  # claim success after failing to drive jamovi.app. Judge by the artifacts.
  [ -f "$ARTIFACT" ] || { echo "error: jmvtools did not produce $ARTIFACT" >&2; rm -f "$LOG"; return 1; }
  if grep -qE 'SingletonLock|sandbox initialization failed|GPU process isn.t usable' "$LOG"; then
    echo
    echo "!! jamovi.app could not be driven (already open, or a sandboxed session)."
    echo "!! The .jmo was built. Install it by hand:"
    echo "!!   jamovi -> Modules -> jamovi library -> Sideload -> $ARTIFACT"
    rm -f "$LOG"; return 0
  fi
  if ! grep -q 'Module installed successfully' "$LOG"; then
    echo "error: jmvtools::install() did not install the module (see above)" >&2
    rm -f "$LOG"; return 1
  fi
  rm -f "$LOG"

  # jamovi unpacks the module shortly AFTER jmvtools reports success. Waiting for the
  # directory is not enough when an older copy is installed: the smoke test would load
  # the old copy (measured: a new analysis "is not an exported object"). Wait for a new
  # Built: stamp instead.
  local NEW=""
  for i in $(seq 1 60); do
    NEW="$(grep '^Built:' "$MODDIR/R/$MODULE/DESCRIPTION" 2>/dev/null || true)"
    [ -n "$NEW" ] && [ "$NEW" != "$OLD" ] && break
    sleep 1
  done
  if [ -z "$NEW" ] || [ "$NEW" = "$OLD" ]; then
    echo "!! desktop: jamovi has not replaced the installed copy of $MODULE after 60 s;"
    echo "!! restart jamovi or sideload $ARTIFACT by hand"
    return 1
  fi
  echo ">> desktop: installed at $MODDIR (${NEW#Built: })"

  if [ -f "$SMOKE" ]; then
    echo ">> desktop: smoke test"
    { smoke_preamble "$APP/Contents/Resources/modules/base/R" "$MODDIR/R"; cat "$SMOKE"; } \
      | Rscript --vanilla -
  fi
}

# ── docker ───────────────────────────────────────────────────────────────────
# jamovi's own route (jamovi >= 28.4): the jamovi compiler on this machine
# builds the module in a throwaway container off the running container's image
# (dependencies included, compiled against the image's R) and hands the .jmo to
# the server over its stdin, which installs it into $HOME/.jamovi/modules in the
# container. Nothing is added to the image. Needs here: docker, node, and a
# jamovi-src checkout (its jamovi-compiler) at the image's version; set
# JMC_COMPILER=/path/to/jamovi-compiler if it is not a sibling of this repo.
# The container must run upstream's docker-compose.yaml (stdin_open and
# --stdin-slave), and the Docker VM must see this repo read-write (it is mounted
# into the build container). To keep modules across 'down'/'up', mount a volume
# at /root/.jamovi (the skill's templates/docker/docker-compose.override.yaml).
find_compiler() {
  local c
  for c in "${JMC_COMPILER:-}" "$ROOT/../jamovi-src/jamovi-compiler" "$ROOT/tools/jamovi-src/jamovi-compiler"; do
    [ -n "$c" ] && [ -f "$c/index.js" ] && [ -f "$c/docker.js" ] && { (cd "$c" && pwd); return 0; }
  done
  return 1
}

install_docker() {
  # A sandbox that cannot read ~/.docker/config.json makes the CLI print a
  # warning and list nothing, which looks like "no container". Use a writable config.
  if [ -z "${DOCKER_CONFIG:-}" ] && ! [ -r "$HOME/.docker/config.json" ]; then
    export DOCKER_CONFIG="$ROOT/.tmp/dockercfg"; mkdir -p "$DOCKER_CONFIG"
  fi
  command -v docker >/dev/null || { echo "!! docker: no docker CLI — skipping"; return 0; }
  if ! docker ps --format '{{.Names}}' 2>/dev/null | grep -qx "$CONTAINER"; then
    echo "!! docker: container '$CONTAINER' is not running — skipping"
    return 0
  fi
  local JMC CHOME OLD NEW i
  JMC="$(find_compiler)" || {
    echo "!! docker: no jamovi-compiler with docker support (jamovi >= 28.4) found." >&2
    echo "!! Clone jamovi/jamovi beside this repo as jamovi-src, or set JMC_COMPILER." >&2
    return 1; }
  command -v node >/dev/null || { echo "!! docker: node is needed on this machine" >&2; return 1; }
  if [ ! -d "$JMC/node_modules" ]; then
    echo ">> docker: installing the compiler's npm dependencies (once)"
    ( cd "$JMC" && npm install --no-audit --no-fund >/dev/null ) || return 1
  fi
  if [ "$(docker inspect -f '{{.Config.OpenStdin}}' "$CONTAINER")" != true ]; then
    echo "!! docker: '$CONTAINER' was not started with stdin open (--stdin-slave);" >&2
    echo "!! start it with jamovi-src's docker-compose.yaml" >&2
    return 1
  fi

  CHOME="$(docker exec "$CONTAINER" sh -c 'echo $HOME')"
  # the build-time stamp of any copy already installed; the new one must differ
  OLD="$(docker exec "$CONTAINER" sh -c "grep '^build-time' '$CHOME/.jamovi/modules/$MODULE/jamovi.yaml' 2>/dev/null" || true)"

  echo ">> docker: building $MODULE $VERSION in $(docker inspect -f '{{.Config.Image}}' "$CONTAINER")"
  node "$JMC/index.js" --install "$MODULE_DIR" --home "docker:$CONTAINER" || return 1
  rm -f "$MODULE_DIR/.jmc-docker.jmo"   # the build container's artifact, already handed over

  # the server installs it asynchronously after reading its stdin
  NEW=""
  for i in $(seq 1 60); do
    NEW="$(docker exec "$CONTAINER" sh -c "grep '^build-time' '$CHOME/.jamovi/modules/$MODULE/jamovi.yaml' 2>/dev/null" || true)"
    [ -n "$NEW" ] && [ "$NEW" != "$OLD" ] && break
    sleep 1
  done
  if [ -z "$NEW" ] || [ "$NEW" = "$OLD" ]; then
    echo "!! docker: jamovi has not installed the new $MODULE after 60 s (docker logs $CONTAINER)" >&2
    return 1
  fi
  echo ">> docker: installed at $CHOME/.jamovi/modules/$MODULE"
  if ! docker inspect -f '{{range .Mounts}}{{println .Destination}}{{end}}' "$CONTAINER" | grep -qx "$CHOME/.jamovi"; then
    echo "!! docker: $CHOME/.jamovi is not a volume: the module is lost when the container is"
    echo "!! recreated (down/up). See the skill's templates/docker/docker-compose.override.yaml"
  fi

  if [ -f "$SMOKE" ]; then
    echo ">> docker: smoke test"
    { smoke_preamble /usr/lib/jamovi/modules/base/R "$CHOME/.jamovi/modules/$MODULE/R"; cat "$SMOKE"; } \
      | docker exec -i "$CONTAINER" bash -c 'source /usr/lib/jamovi/bin/env.conf 2>/dev/null || true; Rscript --vanilla -'
  fi
  echo ">> docker: installed $MODULE $VERSION"
}

case "$TARGET" in
  desktop) install_desktop ;;
  docker)  install_docker ;;
  both)    install_desktop || true; echo; install_docker || true ;;
  *)       echo "usage: install.sh [desktop|docker|both] [container]" >&2; exit 1 ;;
esac
