
RES=$(ufasta dgrep 'CacGACgcggatGaGTcGcagAgccG.*GcGatcGcCtCacAGCAaGTtC' $TEST1)
EXPECT_EQ ">read6" $(echo "$RES" | head -n 1) "Got read6"

RES=$(ufasta dgrep 'acgt' $TEST1 | grep -v '^>')
GR=$(grep 'acgt' $TEST1)
EXPECT_EQ "$GR" "$RES" "Identical to grep"

RES=$(ufasta dgrep -n 'act' $TEST1 | grep -v '>')
GR=$(grep -n 'act' $TEST1)
EXPECT_EQ "$GR" "$RES" "Identical to grep, line numbers"
