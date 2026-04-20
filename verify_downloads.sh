#!/usr/bin/env bash
# Verify every downloaded FAERS ZIP:
#   1. Nonzero and > 1 MB
#   2. `unzip -t` passes (CRC integrity)
#   3. Contains at least one of DEMO* / DRUG* / REAC*
#
# Usage: bash verify_downloads.sh

set -u
cd "$(dirname "$0")/data/raw" 2>/dev/null || {
    echo "error: data/raw/ not found — run download_faers.sh first"
    exit 1
}

total=0; ok=0; bad=0
bad_files=()

for f in faers_ascii_*.zip; do
    [ -f "$f" ] || continue
    total=$((total+1))

    size=$(stat -f%z "$f" 2>/dev/null || stat -c%s "$f")
    if [ "$size" -lt 1048576 ]; then
        echo "[FAIL] $f: too small (${size} bytes)"
        bad=$((bad+1)); bad_files+=("$f")
        continue
    fi

    if ! unzip -tq "$f" >/dev/null 2>&1; then
        echo "[FAIL] $f: corrupt (unzip -t failed)"
        bad=$((bad+1)); bad_files+=("$f")
        continue
    fi

    if ! unzip -l "$f" 2>/dev/null | grep -qiE '(DEMO|DRUG|REAC)[0-9]+Q[0-9]\.(TXT|txt)'; then
        echo "[FAIL] $f: missing expected DEMO/DRUG/REAC files"
        bad=$((bad+1)); bad_files+=("$f")
        continue
    fi

    ok=$((ok+1))
done

echo
echo "Verified:  $ok / $total"
echo "Bad:       $bad"
if [ $bad -gt 0 ]; then
    echo "Re-download these by deleting and re-running download_faers.sh:"
    printf '  %s\n' "${bad_files[@]}"
    exit 1
fi
