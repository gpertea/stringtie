
ufasta n50 -H -N10 -N25 -N50 -N75 -N90 $TEST1 | sort -c -r
EXPECT_SUCCESS "Ns statistics sorted"

S=$(ufasta n50 -H -S $TEST1)
ES=$(grep -v '^>' $TEST1 | perl -ne 'chomp; print' | wc -c)
EXPECT_EQ "$ES" "$S" "Correct sum size"
C=$(grep -c '^>' $TEST1)

EALL=$(cat <<EOF
N10 516
N25 443
N50 314
N90 112
S 9317
A 186.34
E 314.704
C $C
EOF
)
ALL=$(ufasta stats -a -N10 -N25 -N50 -N90 $TEST1)
EXPECT_EQ "$EALL" "$ALL" "All statistics"

FS=$(ufasta sizes $TEST1 | ufasta n50 -f -N10 -N25 -N50 -N90 -E -S -A -C /dev/fd/0)
EXPECT_EQ "$EALL" "$FS" "All statistics from sizes"
