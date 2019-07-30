
ES=$(ufasta n50 -H -S $TEST1)
S=$(ufasta sizes $TEST1 | perl -ane '$c += $F[0]; END { print $c, "\n" }')
EXPECT_EQ "$ES" "$S" "Sum of sizes"

offset=0
ufasta sizes -i $TEST1 | while read s o eo; do
    test $offset -lt $o || fail "Offset not increasing enough"
    offset=$((o + s))
done
EXPECT_SUCCESS "Offset increment"

i=0
ufasta sizes -i $TEST1 | while read s o eo; do
    diff <(ufasta hgrep "^read$i\b" $TEST1 | tail -n +2) \
        <(tail -c +$((o + 1)) $TEST1 | head -c $((eo - o)))
    EXPECT_SUCCESS "Compare sequence $i"
    i=$((i+1))
done
