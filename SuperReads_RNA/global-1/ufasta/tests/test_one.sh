
diff <(grep '^>' $TEST1) <(ufasta one $TEST1 | grep '^>')
EXPECT_SUCCESS "Identical headers"

diff <(grep -v '^>' $TEST1 | perl -ne 'chomp; print') <(ufasta one $TEST1 | grep -v '^>' | perl -ne 'chomp; print')
EXPECT_SUCCESS "Identical sequence"

first=$((echo -e 'toto\ntata'; cat $TEST1) | ufasta one /dev/fd/0 | head -n 2)
expect=$(echo -e 'toto\ntata')
EXPECT_EQ "$expect" "$first" "Keep pre-header"
