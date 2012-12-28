#!/bin/bash
source $(dirname $0)/../common.sh
rm -f test.bai oge.bam oge.bai oge2.bam oge2.bai

$OGE view $DATA/simple.bam -o oge1.bam
[ ! -f oge2.bai ] && err "Failed to find oge1.bai"
mv oge1.bai test.bai
samtools index oge1.bam

[ ! -f test.bai ] && err "Failed to find test.bai"

cmp correct.bai oge.bai || err "BAM index differs (test 1)."

$OGE view $DATA/208.yhet.bam -o oge2.bam
[ ! -f oge2.bai ] && err "Failed to find oge2.bai"
mv oge2.bai test.bai
samtools index oge2.bam

[ ! -f test.bai ] && err "Failed to find test.bai"

cmp correct.bai oge.bai || err "BAM index differs (test 2)."
true
