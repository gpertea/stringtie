
EXPECT_EQ "" "$(ufasta extract $TEST1 | head -n 1)" "Extract empty"
diff $TEST1 <(ufasta extract -v -f /dev/null $TEST1)
EXPECT_SUCCESS "Extract all"

diff <(ufasta hgrep '^read(4|5|7|10|32)\b' $TEST1) \
    <(ufasta extract -n read10 -n read5 -n toto -f <(echo read4 read5 blah read7) -f <(echo read32) -f /dev/null $TEST1)
EXPECT_SUCCESS "Extract 5 reads"

diff <(ufasta hgrep -v '^read(4|5|7|10|32)\b' $TEST1) \
    <(ufasta extract -v -n read10 -n read5 -n toto -f <(echo read4 read5 blah read7) -f <(echo read32) -f /dev/null $TEST1)
EXPECT_SUCCESS "Extract complement 5 reads"

diff <(ufasta hgrep -v -m 4 '^read(4|5|7|10|32)\b' $TEST1) \
    <(ufasta extract -v -m 4 -n read10 -n read5 -n toto -f <(echo read4 read5 blah read7) -f <(echo read32) -f /dev/null $TEST1)
EXPECT_SUCCESS "Extract complement 4 reads"

NBT=$(grep -c '^>' $TEST1)
NBE=$(ufasta extract -p 0.5 $TEST1 | grep -c '^>')
# Small probability to fail
perl -e 'exit(!($ARGV[1] >= $ARGV[0] / 4 && $ARGV[1] <= 3 * $ARGV[0] / 4))' $NBT $NBE
EXPECT_SUCCESS "Extract random number of reads"
