#!/bin/bash
###############################################################################
# conjfdr_step3_check.sh
#
# Check all conjFDR results and summarize loci counts.
# Run interactively: bash conjfdr_step3_check.sh
###############################################################################

BASEDIR=/home/lchang24/links/projects/def-gsarah/lchang24/github/pleiofdr
RESULTSDIR=${BASEDIR}/conjfdr_clean/results

echo "========================================"
echo "conjFDR Results Summary"
echo "========================================"
echo ""

TOTAL_PAIRS=0
COMPLETE=0
FAILED=0
TOTAL_LOCI_05=0
TOTAL_LOCI_01=0

printf "%-45s %8s %8s %s\n" "Pair" "FDR<0.05" "FDR<0.01" "Status"
printf "%-45s %8s %8s %s\n" "----" "--------" "--------" "------"

for d in ${RESULTSDIR}/*/; do
    [ ! -d "$d" ] && continue
    TOTAL_PAIRS=$((TOTAL_PAIRS + 1))
    pair=$(basename "$d")

    loci_file=$(ls "$d"/*_conjfdr_0.05_loci.csv 2>/dev/null | head -1)

    if [ -n "$loci_file" ] && [ -s "$loci_file" ]; then
        n05=$(awk 'END{print NR-1}' "$loci_file")

        # Count loci at FDR < 0.01
        n01=$(awk -F',' 'NR>1 && $NF<0.01' "$loci_file" | wc -l)
        # If comma doesn't work, try tab
        if [ "$n01" -eq 0 ] && [ "$n05" -gt 0 ]; then
            n01=$(awk -F'\t' 'NR>1 && $NF<0.01' "$loci_file" | wc -l)
        fi

        TOTAL_LOCI_05=$((TOTAL_LOCI_05 + n05))
        TOTAL_LOCI_01=$((TOTAL_LOCI_01 + n01))
        COMPLETE=$((COMPLETE + 1))

        printf "%-45s %8d %8d %s\n" "$pair" "$n05" "$n01" "OK"
    else
        FAILED=$((FAILED + 1))
        printf "%-45s %8s %8s %s\n" "$pair" "-" "-" "MISSING"
    fi
done

echo ""
echo "========================================"
printf "Complete: %d / %d pairs\n" $COMPLETE $TOTAL_PAIRS
printf "Failed:   %d\n" $FAILED
printf "Total loci at conjFDR < 0.05: %d\n" $TOTAL_LOCI_05
printf "Total loci at conjFDR < 0.01: %d\n" $TOTAL_LOCI_01
echo "========================================"

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "Check failed pair logs:"
    echo "  ls ${BASEDIR}/conjfdr_clean/logs/conjfdr_*.err"
fi
