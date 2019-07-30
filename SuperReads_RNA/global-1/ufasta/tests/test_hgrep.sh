
diff <(ufasta one $TEST1 | grep -A 1 '^>read1' | grep -v -- --) <(ufasta hgrep '^read1' $TEST1 | ufasta one /dev/fd/0)
EXPECT_SUCCESS "Grep for ^read1"
