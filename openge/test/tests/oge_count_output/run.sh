source $(dirname $0)/../common.sh
rm test.out test2.out

$OGE count $DATA/simple.bam > test.out

[ ! -f test.out ] && err "Failed to find test.out"

grep -q "10" test.out || err "Failed to find expected command output (10)"

$OGE count $DATA/simple.bam $DATA/simple.bam $DATA/simple.sam > test2.out

[ ! -f test2.out ] && err "Failed to find test2.out"

grep -q "30" test2.out || err "Failed to find expected command output (30)"

true
