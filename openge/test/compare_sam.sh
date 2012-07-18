#!/bin/bash
eval $1 > sam_compare_temp.txt
cat sam_compare_temp.txt | grep -vP ^@PG - > sam_compare_a.txt
cat $2 | grep -vP ^@PG - > sam_compare_b.txt
diff sam_compare_a.txt sam_compare_b.txt
ret_val=$?
echo "RET"
echo $ret_val
rm sam_compare_a.txt sam_compare_b.txt sam_compare_temp.txt
exit $ret_val
