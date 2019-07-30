
ori=$(cat $TEST1 | wc)
rc=$(ufasta rc $TEST1 | wc)
EXPECT_EQ "$ori" "$rc" "Output same length"

ori=$(cat $TEST1)
id=$(ufasta rc $TEST1 | ufasta rc /dev/fd/0)
EXPECT_EQ "$ori" "$id" "Involution"

paste <(ufasta rc $TEST1 | ufasta one /dev/fd/0 | grep -v '^>') <(ufasta one $TEST1 | grep -v '^>') | \
    perl -ane '$F[0] = uc($F[0]); $F[1] =~ tr/ACGTacgt/TGCATGCA/; $F[1] = reverse($F[1]); exit(1) unless $F[0] eq $F[1]'
EXPECT_SUCCESS "Reverse complement"

ori=$(cat $TEST1 | wc)
can=$(ufasta rc -C $TEST1 | wc)
EXPECT_EQ "$ori" "$rc" "Output -C same length"

first=$(ufasta rc -C $TEST1)
second=$(echo "$first" | ufasta rc -C /dev/fd/0)
EXPECT_EQ "$first" "$second" "Projection"

paste <(ufasta rc -C $TEST1 | ufasta one /dev/fd/0 | grep -v '^>') <(ufasta one $TEST1 | grep -v '^>') | \
    perl -ane '$F[0] = uc($F[0]); $F[1] = uc($F[1]); $rc = reverse($F[1]); $rc =~ tr/ACGT/TGCA/; unless(($F[0] eq $rc || $F[0] eq $F[1]) && $F[0] cmp $F[1] <= 0) { print("$.\n0 $F[0]\n1 $F[1]\n2 $rc"); exit(1) }'
EXPECT_SUCCESS "Canonical"
