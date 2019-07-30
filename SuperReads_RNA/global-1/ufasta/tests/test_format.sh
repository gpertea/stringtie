
diff <(ufasta one $TEST1) <(ufasta format -l 15 $TEST1 | ufasta one $TEST1)
EXPECT_SUCCESS "Format then one"

ufasta format -l 15 $TEST1 | grep -v '^>' | perl -ne 'chomp; length($_) <= 15 || exit(1)'
EXPECT_SUCCESS "Short lines"

diff <(ufasta one $TEST1 | grep -v '^>' | tr '[:upper:]' '[:lower:]') <(ufasta format -L $TEST1 | ufasta one /dev/fd/0 | grep -v '^>')
EXPECT_SUCCESS "Lower case"

diff <(ufasta one $TEST1 | grep -v '^>' | tr '[:lower:]' '[:upper:]') <(ufasta format -U $TEST1 | ufasta one /dev/fd/0 | grep -v '^>')
EXPECT_SUCCESS "Uppr case"
