
NB=$(grep -c '^>' $TEST1)
for i in $(seq 0 $((NB - 1))); do
    ufasta dsort $TEST1 | ufasta extract -n "read$i" | tail -n +2 | sort -C
    EXPECT_SUCCESS "Sorted entry $i"
done
