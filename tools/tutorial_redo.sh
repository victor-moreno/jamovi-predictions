#!/bin/sh
# Re-render the module's tutorial website with Quarto.
#
#   sh tools/tutorial_redo.sh
#
# Removes the computation cache (_freeze) to guarantee that changes to the
# module's R/ code and .qmd sources are fully re-evaluated, and executes
# Quarto in a sandboxed cache directory to prevent macOS permission errors.
set -eu

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TOOLS_DIR="$ROOT/tools"

if [ -f "$TOOLS_DIR/_module.sh" ]; then
    . "$TOOLS_DIR/_module.sh"
fi

TUTORIAL_DIR="$ROOT/tutorial"
if [ ! -d "$TUTORIAL_DIR" ]; then
    echo "error: no tutorial/ directory found in $ROOT" >&2
    exit 1
fi

# Locate quarto-sandbox.sh (check local tools/, skill directory, or fallback search)
QSANDBOX=""
if [ -f "$TOOLS_DIR/quarto-sandbox.sh" ]; then
    QSANDBOX="$TOOLS_DIR/quarto-sandbox.sh"
elif [ -f "$ROOT/../jamovi-skill/scripts/quarto-sandbox.sh" ]; then
    QSANDBOX="$(cd "$ROOT/../jamovi-skill/scripts" && pwd)/quarto-sandbox.sh"
fi

cd "$TUTORIAL_DIR"

echo ">> cleaning computation cache (_freeze/)..."
rm -rf _freeze

echo ">> rendering tutorial with Quarto..."
QUARTO_SCRATCH="$ROOT/.tmp/qhome"
export QUARTO_SCRATCH

if [ -n "$QSANDBOX" ] && [ -f "$QSANDBOX" ]; then
    sh "$QSANDBOX" render "$@"
elif command -v quarto >/dev/null 2>&1; then
    HOME="$QUARTO_SCRATCH" quarto render "$@"
else
    echo "error: neither quarto-sandbox.sh nor quarto on PATH found." >&2
    exit 1
fi

echo ">> tutorial rendered successfully in: $TUTORIAL_DIR/_site/index.html"
