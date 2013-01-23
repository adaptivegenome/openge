#!/bin/bash
source $(dirname $0)/../common.sh
rm test.bam test.out test2.bam test2.out

$OGE mergesort $DATA/simple.bam -o test.bam

[ ! -f test.bam ] && err "Failed to find simple.sorted.bam"

$OGE stats test.bam > test.out

grep -q "Sorted: *Yes" test.out || err "Failed to find expected command output (Sorted: Yes)"

$OGE mergesort $DATA/simple.bam $DATA/simple.bam $DATA/simple.sam -o test2.bam

[ ! -f test2.bam ] && err "Failed to find test2.out"

$OGE stats test2.bam > test2.out

grep -q "Sorted: *Yes" test2.out || err "Failed to find expected command output (Sorted: Yes)"

true
