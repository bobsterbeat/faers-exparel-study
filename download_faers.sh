#!/usr/bin/env bash
# Download FAERS quarterly ASCII extracts from 2012 Q4 through current.
# - Skips files already downloaded (size > 1 MB sanity check)
# - Retries transient errors (5x, 10s backoff)
# - Resumes interrupted downloads (-C -)
# - 404s / 500s on unreleased quarters are reported but non-fatal
#
# Usage: bash download_faers.sh
# Output: data/raw/faers_ascii_YYYYqN.zip

set -u
cd "$(dirname "$0")"
mkdir -p data/raw
cd data/raw

# Study window: 2012 Q4 (earliest FAERS ASCII available) through 2026 Q1.
# Adjust END_YEAR as new quarters are released.
START_YEAR=2012
START_Q=4
END_YEAR=2026
END_Q=1

ok=0; skipped=0; failed=0; already=0

for year in $(seq $START_YEAR $END_YEAR); do
    for q in 1 2 3 4; do
        if [ $year -eq $START_YEAR ] && [ $q -lt $START_Q ]; then continue; fi
        if [ $year -eq $END_YEAR ]   && [ $q -gt $END_Q ];   then continue; fi

        fname="faers_ascii_${year}q${q}.zip"
        url="https://fis.fda.gov/content/Exports/${fname}"

        # Skip if already downloaded and reasonably sized (>1MB)
        if [ -f "$fname" ] && [ "$(stat -f%z "$fname" 2>/dev/null || stat -c%s "$fname")" -gt 1048576 ]; then
            echo "[skip] $fname (already downloaded)"
            already=$((already+1))
            continue
        fi

        echo "[get ] $fname"
        # --retry-all-errors covers curl exit 18 (partial transfer), which
        # FDA triggers regularly by dropping mid-download. -C - resumes
        # from the byte offset of any partial file — so we keep partials
        # on failure (re-running the script will pick them up).
        if curl -fsSL \
                --retry 10 --retry-delay 15 --retry-all-errors \
                --retry-connrefused \
                -C - \
                -o "$fname" \
                "$url"; then
            ok=$((ok+1))
            sleep 3  # be polite to FDA (they appear to rate-limit)
        else
            rc=$?
            echo "       FAILED (curl exit $rc) — keeping partial for resume"
            failed=$((failed+1))
        fi
    done
done

echo
echo "Downloaded:      $ok"
echo "Already present: $already"
echo "Failed/missing:  $failed"
echo "Files in $(pwd):"
ls -lh faers_ascii_*.zip 2>/dev/null | awk '{print "  " $9, $5}'
