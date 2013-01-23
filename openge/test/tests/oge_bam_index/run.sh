#!/bin/bash
source $(dirname $0)/../common.sh
rm -f test.bam.bai oge.bam oge.bam.bai oge2.bam oge2.bam.bai

$OGE view $DATA/simple.bam -o oge1.bam
[ ! -f oge1.bam.bai ] && err "Failed to find oge1.bam.bai"
mv oge1.bam.bai test.bai
samtools index oge1.bam

[ ! -f test.bai ] && err "Failed to find test.bai"

cmp test.bai oge1.bam.bai || err "BAM index differs (test 1)."

#test with bigger, unsorted input
$OGE mergesort $DATA/208.yhet.bam -o oge2.bam
[ ! -f oge2.bam.bai ] && err "Failed to find oge2.bam.bai"
mv oge2.bam.bai test.bai
samtools index oge2.bam

[ ! -f test.bai ] && err "Failed to find test.bai"

cmp oge2.bam.bai test.bai || err "BAM index differs (test 2)."
true
