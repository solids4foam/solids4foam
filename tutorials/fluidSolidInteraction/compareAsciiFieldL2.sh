#!/bin/bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <referenceField> <candidateField>" >&2
    exit 1
fi

ref_file=$1
cand_file=$2

ref_tmp=$(mktemp)
cand_tmp=$(mktemp)
trap 'rm -f "$ref_tmp" "$cand_tmp"' EXIT

extract_internal_field()
{
    awk 'function emitValues(line,    n, parts, i) { gsub(/[()]/, " ", line); gsub(/^[[:space:]]+|[[:space:]]+$/, "", line); if (line == "") return; n = split(line, parts, /[[:space:]]+/); for (i = 1; i <= n; ++i) if (parts[i] != "") print parts[i] } /^[[:space:]]*internalField[[:space:]]+uniform[[:space:]]+/ { sub(/.*internalField[[:space:]]+uniform[[:space:]]+/, "", $0); sub(/[[:space:]]*;[[:space:]]*$/, "", $0); emitValues($0); exit } /^[[:space:]]*internalField[[:space:]]+nonuniform/ { state = "count"; next } state == "count" && /^[[:space:]]*[0-9]+[[:space:]]*$/ { state = "open"; next } state == "open" && /^[[:space:]]*\([[:space:]]*$/ { state = "data"; next } state == "data" && /^[[:space:]]*\)[[:space:]]*;?[[:space:]]*$/ { exit } state == "data" { emitValues($0) }' "$1"
}

extract_internal_field "$ref_file" > "$ref_tmp"
extract_internal_field "$cand_file" > "$cand_tmp"

ref_count=$(wc -l < "$ref_tmp" | tr -d ' ')
cand_count=$(wc -l < "$cand_tmp" | tr -d ' ')

if [[ "$ref_count" != "$cand_count" ]]; then
    echo \
        "Field size mismatch: ${ref_file} has ${ref_count} numeric entries, "\
        "${cand_file} has ${cand_count}" >&2
    exit 2
fi

paste "$ref_tmp" "$cand_tmp" | awk -v ref="$ref_file" -v cand="$cand_file" 'NR>=1 { diff = $2 - $1; sumDiff2 += diff * diff; sumRef2 += $1 * $1; if (diff < 0) diff = -diff; if (diff > maxAbs) maxAbs = diff } END { if (sumRef2 > 0) rel = sqrt(sumDiff2 / sumRef2); else rel = sqrt(sumDiff2); printf "reference=%s candidate=%s components=%d relative_l2=%.12g max_abs=%.12g\n", ref, cand, NR, rel, maxAbs }'
