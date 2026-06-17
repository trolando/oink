#!/usr/bin/env bash
#
# Smoke benchmark for Oink: run the `oink` CLI over a corpus of parity games,
# reporting the read / parse / solve times and verifying each solution. It is
# meant to catch obvious correctness and performance regressions, not to be a
# scientific benchmark.
#
# Usage:
#   scripts/benchmark-smoke.sh [corpus-dir-or-file] [solver]
#
# Environment:
#   OINK     path to the oink CLI (default: build-release/oink, then build-debug/oink)
#   TIMEOUT  per-game timeout in seconds, passed to oink -z (default: 60)
#
# Exit status is non-zero if any game fails to verify (or times out).

set -u

CORPUS="${1:-tests}"
SOLVER="${2:-tl}"
TIMEOUT="${TIMEOUT:-60}"

# Locate the oink CLI.
OINK="${OINK:-}"
if [ -z "$OINK" ]; then
    for candidate in build-release/oink build-debug/oink; do
        if [ -x "$candidate" ]; then OINK="$candidate"; break; fi
    done
fi
if [ -z "$OINK" ] || [ ! -x "$OINK" ]; then
    echo "error: oink CLI not found; build it or set OINK=/path/to/oink" >&2
    exit 2
fi

if [ ! -e "$CORPUS" ]; then
    echo "error: corpus not found: $CORPUS" >&2
    exit 2
fi

echo "oink:   $OINK"
echo "solver: --$SOLVER   timeout: ${TIMEOUT}s   corpus: $CORPUS"
echo
printf '%-30s %9s %9s %9s  %s\n' "file" "read(s)" "parse(s)" "solve(s)" "result"

games=0; failed=0; skipped=0
tot_read=0; tot_parse=0; tot_solve=0

extract() { sed -n "s/.*$1 \([0-9.]*\) sec.*/\1/p" | head -1; }

while IFS= read -r f; do
    out=$("$OINK" -i "$f" "--$SOLVER" -v -z "$TIMEOUT" 2>&1)
    if printf '%s' "$out" | grep -q "parsing error"; then
        skipped=$((skipped + 1))
        continue
    fi
    read_t=$(printf '%s' "$out" | extract "reading input took");  read_t=${read_t:-0}
    parse_t=$(printf '%s' "$out" | extract "parsing took");       parse_t=${parse_t:-0}
    solve_t=$(printf '%s' "$out" | extract "total solving time:"); solve_t=${solve_t:-0}
    if printf '%s' "$out" | grep -q "solution verified"; then
        result="verified"
    else
        result="FAILED"
        failed=$((failed + 1))
    fi
    games=$((games + 1))
    tot_read=$(awk "BEGIN{print $tot_read + $read_t}")
    tot_parse=$(awk "BEGIN{print $tot_parse + $parse_t}")
    tot_solve=$(awk "BEGIN{print $tot_solve + $solve_t}")
    printf '%-30s %9s %9s %9s  %s\n' "$(basename "$f")" "$read_t" "$parse_t" "$solve_t" "$result"
done < <(if [ -d "$CORPUS" ]; then find "$CORPUS" -maxdepth 1 -type f | sort; else echo "$CORPUS"; fi)

echo
printf '%-30s %9.3f %9.3f %9.3f  %s\n' "TOTAL ($games games)" \
    "$tot_read" "$tot_parse" "$tot_solve" \
    "$([ "$failed" -eq 0 ] && echo 'all verified' || echo "$failed FAILED")"
[ "$skipped" -gt 0 ] && echo "($skipped non-parity-game files skipped)"

[ "$failed" -eq 0 ]
